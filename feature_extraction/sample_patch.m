function resized_patch = sample_patch(im, pos, sample_sz, output_sz, gparams)

if nargin < 4
    output_sz = [];
end
if nargin < 5 || ~isfield(gparams, 'use_mexResize')
    gparams.use_mexResize = false;
end

% Pos should be integer when input, but floor in just in case.
pos = floor(pos);

% Downsample factor
resize_factor = min(sample_sz ./ output_sz);
df = max(floor(resize_factor - 0.1), 1);
if df > 1
    % pos = 1 + of + df * (npos - 1)
    
    % compute offset and new center position
    os = mod(pos - 1, df);
    pos = (pos - 1 - os) / df + 1;
    
    % new sample size
    sample_sz = sample_sz / df;
    
    % donwsample image
    im = im(1+os(1):df:end, 1+os(2):df:end, :);
end

% make sure the size is not too small and round it
sample_sz = max(round(sample_sz), 2);

xs = pos(2) + (1:sample_sz(2)) - floor((sample_sz(2)+1)/2);
ys = pos(1) + (1:sample_sz(1)) - floor((sample_sz(1)+1)/2);

%check for out-of-bounds coordinates, and set them to the values at
%the borders
xs(xs < 1) = 1;
ys(ys < 1) = 1;
xs(xs > size(im,2)) = size(im,2);
ys(ys > size(im,1)) = size(im,1);

%extract image
im_patch = im(ys, xs, :);

if isempty(output_sz) || isequal(sample_sz(:), output_sz(:))
    resized_patch = im_patch;
else
    if gparams.use_mexResize
        resized_patch = mexResize(im_patch, output_sz, 'linear');
    else
        resized_patch = imresize(im_patch, output_sz, 'bilinear', 'Antialiasing',false);
    end
end

end

