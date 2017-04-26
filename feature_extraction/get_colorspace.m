function [ colorspace ] = get_colorspace( im, fparam, gparam )

% Get a color space feature. Currently implements 'gray' and 'rgb'.

[im_height, im_width, num_im_chan, num_images] = size(im);

single_im = single(im)/255;

if strcmp(fparam.colorspace,'gray')
    if num_im_chan == 3
        if num_images == 1
            t_colorspace = rgb2gray(single_im) - 0.5;
        else
            t_colorspace = zeros(im_height, im_width, 1, num_images, 'single');
            for k = 1:num_images
                t_colorspace(:,:,:,k) = rgb2gray(single_im(:,:,:,k)) - 0.5;
            end
        end
    elseif num_im_chan == 1
        t_colorspace = single_im - 0.5;
    else
        except = MException('get_colorspace','Invalid input data, must have 1 or 3 dimensions');
        throw(except);
    end;
elseif strcmp(fparam.colorspace,'rgb')
    if num_im_chan == 3
        t_colorspace = single_im - 0.5;
    else
        except = MException('get_colorspace','Invalid input data, must have 3 dimensions for rgb');
        throw(except);
    end;
end;

if gparam.cell_size > 1
    colorspace = average_feature_region(t_colorspace,gparam.cell_size);
else
    colorspace = t_colorspace;
end
end

