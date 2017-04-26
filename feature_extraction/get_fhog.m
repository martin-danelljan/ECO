function [ feature_image ] = get_fhog( im, fparam, gparam )
% Extract FHOG features using pdollar toolbox.
% Currently assumes 9 gradient orientations.

if ~isfield(fparam, 'nOrients')
    fparam.nOrients = 9;
end

[im_height, im_width, num_im_chan, num_images] = size(im);
feature_image = zeros(floor(im_height/gparam.cell_size), floor(im_width/gparam.cell_size), fparam.nDim, num_images, 'single');

for k = 1:num_images
    hog_image = fhog(single(im(:,:,:,k)), gparam.cell_size, fparam.nOrients);
    
    %the last dimension is all 0 so we can discard it
    feature_image(:,:,:,k) = hog_image(:,:,1:end-1);
end
end