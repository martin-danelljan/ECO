function extract_info = get_feature_extract_info(features)

% Find used image sample size
extract_info.img_sample_sizes = {};
extract_info.img_input_sizes = {};
for feat_ind = 1:length(features)
    % if not equals any previously stored size
    if ~any(cellfun(@(sz) isequal(features{feat_ind}.img_sample_sz, sz), extract_info.img_sample_sizes))
        extract_info.img_sample_sizes{end+1} = features{feat_ind}.img_sample_sz;
        extract_info.img_input_sizes{end+1} = features{feat_ind}.img_input_sz;
    end
end