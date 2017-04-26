function [ region_image ] = average_feature_region(im, region_size)
% compute cell wise averages, where a cell is a region_size*region_sized
% region in the image. Input can be uint8t, single or double matrices
% of arbitrary dimension

region_area = region_size.^2;

if isa(im,'double') || isa(im,'single') || isa(im,'gpuArray')
    maxval = 1.0;
elseif isa(im,'unit8')
    maxval = 255;
end

% compute the integral image
iImage = integralVecImage(im);

% region indices
i1 = (region_size:region_size:size(im,1)) + 1;
i2 = (region_size:region_size:size(im,2)) + 1;

% sum over region, divided by number of elements, and normalize to [0,1]
% range if integer image
region_image = (iImage(i1,i2,:,:) - iImage(i1,i2-region_size,:,:) - iImage(i1-region_size,i2,:,:) + iImage(i1-region_size,i2-region_size,:,:)) ./ (region_area * maxval);

end

