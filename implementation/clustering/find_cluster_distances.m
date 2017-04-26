function dist_vector = find_cluster_distances(samplesf, xlf_proj_perm, num_feature_blocks, num_training_samples, max_train_samples, params)

dist_vector = inf(max_train_samples,1,'single');

if params.neglect_higher_frequency
    if num_training_samples == max_train_samples
        sample_distance = 0;
        for k = 1:num_feature_blocks
            tmp_sample_sz = size(xlf_proj_perm{k});
            tmp_r1 = round(1 + tmp_sample_sz(3)*0.5*params.frequency_clip_factor);
            tmp_r2 = tmp_r1 + round(tmp_sample_sz(3)*(1-params.frequency_clip_factor)) - 1;
            
            tmp_c1 = round(1 + tmp_sample_sz(4)*params.frequency_clip_factor);
            
            diff_matrix = permute(reshape(bsxfun(@minus, xlf_proj_perm{k}(:,:,tmp_r1:tmp_r2,tmp_c1:end), samplesf{k}(:,:,tmp_r1:tmp_r2,tmp_c1:end)), num_training_samples, []), [2 3 1]);
            sample_distance = sample_distance + real(squeeze(mtimesx(diff_matrix, 'C', diff_matrix, 'speed')));
        end
        
        dist_vector(1:num_training_samples) = sqrt(sample_distance);
    elseif num_training_samples > 0
        sample_distance = 0;
        for k = 1:num_feature_blocks
            tmp_sample_sz = size(xlf_proj_perm{k});
            tmp_r1 = round(1 + tmp_sample_sz(3)*0.5*params.frequency_clip_factor);
            tmp_r2 = tmp_r1 + round(tmp_sample_sz(3)*(1-params.frequency_clip_factor)) - 1;
            
            tmp_c1 = round(1 + tmp_sample_sz(4)*params.frequency_clip_factor);
            
            diff_matrix = permute(reshape(bsxfun(@minus, xlf_proj_perm{k}(:,:,tmp_r1:tmp_r2,tmp_c1:end), samplesf{k}(1:num_training_samples,:,tmp_r1:tmp_r2,tmp_c1:end)), num_training_samples, []), [2 3 1]);
            sample_distance = sample_distance + real(squeeze(mtimesx(diff_matrix, 'C', diff_matrix, 'speed')));
        end
        
        dist_vector(1:num_training_samples) = sqrt(sample_distance);
    end
else
    if num_training_samples == max_train_samples
        % This if statement is only for speed
        sample_distance = 0;
        for k = 1:num_feature_blocks
            diff_matrix = permute(reshape(bsxfun(@minus, xlf_proj_perm{k}, samplesf{k}), num_training_samples, []), [2 3 1]);
            sample_distance = sample_distance + real(squeeze(mtimesx(diff_matrix, 'C', diff_matrix, 'speed')));
        end
        
        dist_vector(1:num_training_samples) = sqrt(sample_distance);
    elseif num_training_samples > 0
        sample_distance = 0;
        for k = 1:num_feature_blocks
            diff_matrix = permute(reshape(bsxfun(@minus, xlf_proj_perm{k}, samplesf{k}(1:num_training_samples,:,:,:)), num_training_samples, []), [2 3 1]);
            sample_distance = sample_distance + real(squeeze(mtimesx(diff_matrix, 'C', diff_matrix, 'speed')));
        end
        
        dist_vector(1:num_training_samples) = sqrt(sample_distance);
    end
end

end

