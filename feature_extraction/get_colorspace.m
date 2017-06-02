function [ feature_map ] = get_colorspace( im, fparam, gparam )

% Get a color space feature. Currently implements 'gray' and 'rgb'.

if isfield(fparam, 'cell_size')
    cell_size = fparam.cell_size;
else
    cell_size = gparam.cell_size;
end

[im_height, im_width, num_im_chan, num_images] = size(im);

single_im = single(im)/255;

if strcmpi(fparam.colorspace,'gray')
    if num_im_chan == 3
        if num_images == 1
            feature_map = rgb2gray(single_im) - 0.5;
        else
            feature_map = zeros(im_height, im_width, 1, num_images, 'single');
            for k = 1:num_images
                feature_map(:,:,:,k) = rgb2gray(single_im(:,:,:,k)) - 0.5;
            end
        end
    elseif num_im_chan == 1
        feature_map = single_im - 0.5;
    else
        except = MException('get_colorspace','Invalid input data, must have 1 or 3 dimensions');
        throw(except);
    end
elseif strcmpi(fparam.colorspace,'rgb')
    if num_im_chan == 3
        feature_map = single_im - 0.5;
    else
        except = MException('get_colorspace','Invalid input data, must have 3 dimensions for rgb');
        throw(except);
    end
end

if gparams.use_gpu
    feature_map = gpuArray(feature_map);
end

if cell_size > 1
    feature_map = average_feature_region(feature_map, cell_size);
end
end

