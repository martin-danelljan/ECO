function merged_sample = merge_samples(sample1, sample2, w1, w2, sample_merge_type)
% Merge sample1 and sample2 using weights w1 and w2. The type of merging is
% decided by sample_merge_type. 
% The sample_merge_type can be
% 1) Merge: The output is the weighted sum of the input samples
% 2) Replace: The output is the first sample. ie w2 is assumed to be 0


% Normalise the weights so that they sum to one
alpha1 = w1/(w1+w2);
alpha2 = 1 - alpha1;

% Build the merged sample
if strcmpi(sample_merge_type, 'replace')
    merged_sample = sample1;
elseif strcmpi(sample_merge_type, 'merge')
    num_feature_blocks = numel(sample1);
    
    merged_sample = cell(1, 1, num_feature_blocks);
    for k = 1:num_feature_blocks
        merged_sample{k} = alpha1*sample1{k} + alpha2*sample2{k};
    end
else
    error('Invalid sample merge type');
end
