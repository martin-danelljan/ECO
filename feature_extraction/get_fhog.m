function [ feature_image ] = get_fhog( im, fparam, gparam )

% Extract FHOG features using pdollar toolbox.

if ~isfield(fparam, 'nOrients')
    fparam.nOrients = 9;
end
if isfield(fparam, 'cell_size')
    cell_size = fparam.cell_size;
else
    cell_size = gparam.cell_size;
end

[im_height, im_width, num_im_chan, num_images] = size(im);
feature_image = zeros(floor(im_height/cell_size), floor(im_width/cell_size), fparam.nDim, num_images, 'like', gparam.data_type);

for k = 1:num_images
    hog_image = fhog(single(im(:,:,:,k)), cell_size, fparam.nOrients);
    
    %the last dimension is all 0 so we can discard it
    feature_image(:,:,:,k) = hog_image(:,:,1:end-1);
end
end