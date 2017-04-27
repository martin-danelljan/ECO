function feature_map = extract_features(image, pos, scales, features, gparams, extract_info)

% Sample image patches at given position and scales. Then extract features
% from these patches.
% Requires that cell size and image sample size is set for each feature.

if ~iscell(features)
    error('Wrong input');
end;

num_features = length(features);
num_scales = length(scales);
num_sizes = length(extract_info.img_sample_sizes);

% Extract image patches
img_samples = cell(2,1);
for sz_ind = 1:num_sizes
    img_sample_sz = extract_info.img_sample_sizes{sz_ind};
    img_input_sz = extract_info.img_input_sizes{sz_ind};
    img_samples{sz_ind} = zeros(img_input_sz(1), img_input_sz(2), size(image,3), num_scales, 'uint8');
    for scale_ind = 1:num_scales
        img_samples{sz_ind}(:,:,:,scale_ind) = sample_patch(image, pos, round(img_sample_sz*scales(scale_ind)), img_input_sz, gparams);
    end
end

% Find the number of feature blocks and total dimensionality
num_feature_blocks = 0;
total_dim = 0;
for feat_ind = 1:num_features
    num_feature_blocks = num_feature_blocks + length(features{feat_ind}.fparams.nDim);
    total_dim = total_dim + sum(features{feat_ind}.fparams.nDim);
end

feature_map = cell(1, 1, num_feature_blocks);

% Extract feature maps for each feature in the list
ind = 1;
for feat_ind = 1:num_features
    feat = features{feat_ind};
    gparams.cell_size = feat.fparams.cell_size;
    
    % get the image patch index
    img_sample_ind = cellfun(@(sz) isequal(feat.img_sample_sz, sz), extract_info.img_sample_sizes);
    
    % do feature computation
    if feat.is_cell
        num_blocks = length(feat.fparams.nDim);
        feature_map(ind:ind+num_blocks-1) = feat.getFeature(img_samples{img_sample_ind}, feat.fparams, gparams);
    else
        num_blocks = 1;
        feature_map{ind} = feat.getFeature(img_samples{img_sample_ind}, feat.fparams, gparams);
    end
    
    ind = ind + num_blocks;
end
              
% Do feature normalization
if ~isempty(gparams.normalize_power) && gparams.normalize_power > 0
    if gparams.normalize_power == 2
        feature_map = cellfun(@(x) bsxfun(@times, x, ...
            sqrt((size(x,1)*size(x,2))^gparams.normalize_size * size(x,3)^gparams.normalize_dim ./ ...
            (sum(reshape(x, [], 1, 1, size(x,4)).^2, 1) + eps))), ...
            feature_map, 'uniformoutput', false);
    else
        feature_map = cellfun(@(x) bsxfun(@times, x, ...
            ((size(x,1)*size(x,2))^gparams.normalize_size * size(x,3)^gparams.normalize_dim ./ ...
            (sum(abs(reshape(x, [], 1, 1, size(x,4))).^gparams.normalize_power, 1) + eps)).^(1/gparams.normalize_power)), ...
            feature_map, 'uniformoutput', false);
    end
end
if gparams.square_root_normalization
    feature_map = cellfun(@(x) sign(x) .* sqrt(abs(x)), feature_map, 'uniformoutput', false);
end

end