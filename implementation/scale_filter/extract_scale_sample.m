function scale_sample = extract_scale_sample(im, pos, base_target_sz, scaleFactors, scale_model_sz, use_mexResize)

% Get scale filter sample.

if nargin < 6
    use_mexResize = true;
end

nScales = length(scaleFactors);

% Downsample factor.
df = floor(min(scaleFactors));
if df > 1
    im = im(1:df:end,1:df:end,:);
    pos = (pos - 1) / df + 1;
    scaleFactors = scaleFactors / df;
end

for s = 1:nScales
    patch_sz = floor(base_target_sz * scaleFactors(s));
    
    xs = floor(pos(2)) + (1:patch_sz(2)) - floor(patch_sz(2)/2);
    ys = floor(pos(1)) + (1:patch_sz(1)) - floor(patch_sz(1)/2);
    
    %check for out-of-bounds coordinates, and set them to the values at
    %the borders
    xs(xs < 1) = 1;
    ys(ys < 1) = 1;
    xs(xs > size(im,2)) = size(im,2);
    ys(ys > size(im,1)) = size(im,1);
    
    %extract image
    im_patch = im(ys, xs, :);
    
    % resize image to model size
    if use_mexResize
        im_patch_resized = mexResize(im_patch, scale_model_sz, 'auto');
    else
        im_patch_resized = imresize(im_patch, scale_model_sz, 'bilinear', 'Antialiasing',false);
    end
    
    % extract scale features
    temp_hog = fhog(single(im_patch_resized), 4);
    
    if s == 1
        dim_scale = size(temp_hog,1)*size(temp_hog,2)*31;
        scale_sample = zeros(dim_scale, nScales, 'single');
    end
    
    scale_sample(:,s) = reshape(temp_hog(:,:,1:31), dim_scale, 1);
end
