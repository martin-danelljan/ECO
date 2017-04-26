function [score_matrix, new_sample] = update_score_matrix(sample1, sample2, id1, id2, w1, w2, ...
    weight_update_criteria,sample_update_criteria, score_matrix, dist_vector, samplesf, num_feature_blocks)

% Update the weight matrix and find the new sample to be added

alpha1 = w1/(w1+w2);
alpha2 = 1 - alpha1;

% Build the new sample
if strcmp(sample_update_criteria, 'Replace')
    new_sample = sample2;
elseif strcmp(sample_update_criteria, 'Merge')
    new_sample = sample1;
    for k = 1:num_feature_blocks
        new_sample{k} = alpha1*sample1{k} + alpha2*sample2{k};
    end
else
    error('Invalid sample update criteria');
end

% Check whether to use dist_vector
if id2 > 0
    dist_vector = score_matrix(:,id2);
end

if strcmp(weight_update_criteria, 'Add')
    score_matrix(:,id1) = score_matrix(:,id1) + dist_vector;
    score_matrix(id1,:) = score_matrix(:,id1) ;
    
    score_matrix(id1,id1) = inf;
elseif strcmp(weight_update_criteria, 'WeightedAdd')
    score_matrix(:,id1) = alpha1*score_matrix(:,id1) + alpha2*dist_vector;
    score_matrix(id1,:) = score_matrix(:,id1) ;
    
    score_matrix(id1,id1) = inf;
elseif strcmp(weight_update_criteria, 'Replace')
    score_matrix(:,id1) = dist_vector;
    score_matrix(id1,:) = score_matrix(:,id1) ;
    
    score_matrix(id1,id1) = inf;
elseif strcmp(weight_update_criteria, 'Recompute')
    
    sample_distance = 0;
    for k = 1:num_feature_blocks
        diff_matrix = abs(bsxfun(@minus, new_sample{k}, samplesf{k})).^2;
        sample_distance = sample_distance + sum(sum(sum(diff_matrix,2),3),4);
    end
    
    if(numel(sample_distance) ~=size(score_matrix,1))
        error('Size of new distance vector incorrect');
    end
    
    score_matrix(:,id1) = sqrt(sample_distance);
    score_matrix(id1,:) = sqrt(sample_distance);
    
    score_matrix(id1,id1) = inf;
else
    error('Invalid weight update criteria');
end
