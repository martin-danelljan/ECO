function gram_vector = find_gram_vector_gpu(samplesf, new_sample, num_training_samples, params)
% Find the inner product of the new sample with the existing samples. TO be
% used for distance calculation

% Note: Since getting the 'exact' distance between the samples is not important,
% the inner product computation is done by using only the half spectrum for
% efficiency. In practice the error incurred by this is negligible. Also
% since we wannt to merge samples that are similar, and finding the best
% match is not important, small error in the distance computation doesn't
% matter

gram_vector = inf(params.nSamples,1,'single');

num_feature_blocks = numel(new_sample);

if num_training_samples > 0
    ip = zeros(1,num_training_samples,'like', params.data_type_complex);
    for k = 1:num_feature_blocks
        ip = ip + 2 * (new_sample{k}(:)' * reshape(samplesf{k}(:,:,:,1:num_training_samples), [],num_training_samples));
    end
    
    gram_vector(1:num_training_samples) = real(gather(ip));
end

end


