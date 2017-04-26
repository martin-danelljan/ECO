function [merged_cluster, new_cluster, merged_cluster_id, new_cluster_id, score_matrix, prior_weights, num_training_samples] = ...
    merge_clusters(samplesf, xlf_proj_perm, dist_vector, score_matrix, prior_weights,...
    num_training_samples,num_feature_blocks,max_train_samples,minimum_sample_weight,params)

merged_cluster = [];
new_cluster = [];
merged_cluster_id = -1;
new_cluster_id = -1;

% Check if we have filled the memory
if num_training_samples == max_train_samples && params.use_sample_merge
    
    % Check if any sample weight is too low
    [min_sample_weight, min_sample_id] = min(prior_weights);
    
    if min_sample_weight < minimum_sample_weight
        score_matrix(:,min_sample_id) = dist_vector;
        score_matrix(min_sample_id,:) = dist_vector;
        score_matrix(min_sample_id,min_sample_id) = inf;
        
        prior_weights(min_sample_id) = 0;
        prior_weights = prior_weights*(1 - params.learning_rate)/sum(prior_weights);
        prior_weights(min_sample_id) = params.learning_rate;
        
        new_cluster_id = min_sample_id;
        
        new_cluster = xlf_proj_perm;
        
        % cluster_frame_ids{replace_ind} = frame;
        
    else
        % Find closest match with the new sample
        [new_sample_min_dist, new_sample_closest_sample] = min(dist_vector);
        
        % Find the closest pair amongst existing samples
        [min_samples_dist, closest_samples_id] = min(score_matrix(:));
        [sample_id1,sample_id2] = ind2sub(size(score_matrix),closest_samples_id);
        
        if sample_id1 == sample_id2
            error('Score matrix diagonal filled wrongly');
        end
        
        if new_sample_min_dist < min_samples_dist
            % Renormalize sample weights
            prior_weights = prior_weights*(1 - params.learning_rate);
            
            % Set replace index
            replace_ind = new_sample_closest_sample;
            
            replace_sample = cell(1, 1, num_feature_blocks);
            
            % Extract the old sample
            for k = 1:num_feature_blocks
                replace_sample{k} = samplesf{k}(replace_ind,:,:,:);
            end
            
            % Update scores
            [score_matrix, merged_cluster] = update_score_matrix(replace_sample, xlf_proj_perm, replace_ind, -1, prior_weights(replace_ind),params.learning_rate,...
                params.weight_update_criteria, params.sample_update_criteria, score_matrix, dist_vector, samplesf,num_feature_blocks);
            
            prior_weights(new_sample_closest_sample) = prior_weights(new_sample_closest_sample) + params.learning_rate;
            
            merged_cluster_id = replace_ind;
            
            % cluster_frame_ids{replace_ind} = [cluster_frame_ids{replace_ind}; frame];
        else
            prior_weights = prior_weights*(1 - params.learning_rate);
            
            % Decide which sample to keep by checking the weights!!
            if prior_weights(sample_id2) > prior_weights(sample_id1)
                temp = sample_id1;
                sample_id1 = sample_id2;
                sample_id2 = temp;
            end
                       
            replace_ind = sample_id2;
            
            sample1 = cell(1, 1, num_feature_blocks);
            sample2 = cell(1, 1, num_feature_blocks);
            
            % Extract the old sample
            for k = 1:num_feature_blocks
                sample1{k} = samplesf{k}(sample_id1,:,:,:);
                sample2{k} = samplesf{k}(sample_id2,:,:,:);
            end
            
            % Update scores
            [score_matrix, merged_cluster] = update_score_matrix(sample1, sample2, sample_id1,sample_id2, prior_weights(sample_id1),prior_weights(sample_id2),...
                params.weight_update_criteria, params.sample_update_criteria, score_matrix, dist_vector, samplesf,num_feature_blocks);
            
            % Enter new score
            score_matrix(:,sample_id2) = dist_vector;
            score_matrix(sample_id2, :) = dist_vector;
            score_matrix(sample_id2,sample_id2) = inf;
            
            % Add sample weights
            prior_weights(sample_id1) = prior_weights(sample_id1) + prior_weights(sample_id2);
            prior_weights(sample_id2) = params.learning_rate;
            
            merged_cluster_id = sample_id1;
            new_cluster_id = replace_ind;
            
            new_cluster = xlf_proj_perm;
            % cluster_frame_ids{sample_id1} = [cluster_frame_ids{sample_id1}; cluster_frame_ids{sample_id2}];
            % cluster_frame_ids{replace_ind} = frame;
            
        end
    end
else
    % Enter new score
    num_training_samples = num_training_samples + 1;
    score_matrix(:,num_training_samples) = dist_vector;
    score_matrix(num_training_samples, :) = dist_vector;
    score_matrix(num_training_samples, num_training_samples) = inf;
    
    % Enter new weight
    if num_training_samples == 1
        prior_weights(num_training_samples) = 1;
    else
        prior_weights = prior_weights*(1 - params.learning_rate);
        prior_weights(num_training_samples) = params.learning_rate;
    end
    replace_ind = num_training_samples;
    new_cluster_id = replace_ind;
    
    new_cluster = xlf_proj_perm;
    
    % cluster_frame_ids{replace_ind} = frame;
end

if (abs(1 - sum(prior_weights)) > 1e-5)
    error('Weights not properly updated');
end
end