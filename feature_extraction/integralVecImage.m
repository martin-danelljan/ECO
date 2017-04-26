function intImage = integralVecImage(I)
%integralImage Compute integral image.
%   J = integralImage(I) computes integral image of an intensity image I.
%   The output integral image, J, is zero padded on top and left, resulting
%   in size(J) = size(I) + 1. This facilitates easy computation of pixel 
%   sums along all image boundaries. Integral image, J, is essentially 
%   a padded version of CUMSUM(CUMSUM(I),2).
%
%   Class Support
%   -------------
%   Intensity image I can be any numeric class. The class of output 
%   integral image, J, is double.
%
%   Example
%   -------
%   % Compute the integral image and use it to compute sum of pixels
%   % over a rectangular region in I.
%   I = magic(5)
%   
%   % define rectangular region as 
%   % [startingRow, startingColumn, endingRow, endingColumn]
%   [sR sC eR eC] = deal(1, 3, 2, 4);
%    
%   % compute the sum over the region using the integral image
%   J = integralImage(I);
%   regionSum = J(eR+1,eC+1) - J(eR+1,sC) - J(sR,eC+1) + J(sR,sC)
%
%   See also integralFilter, CUMSUM

%   Copyright 2010 The MathWorks, Inc.

%   References:
%      P.A. Viola and M.J. Jones. Rapid object detection using boosted
%      cascade of simple features. In CVPR (1), pages 511-518, 2001.
%
%#codegen
%#ok<*EMCLS>
%#ok<*EMCA>

% validateattributes(I, {'numeric','logical'}, {'2d', 'nonsparse', 'real'},...
%     'integralImage', 'I');

if ~isempty(I)
    if isa(I,'gpuArray')
        intImage = zeros(size(I,1)+1, size(I,2)+1, size(I,3), size(I,4), 'single','gpuArray');
    else
        intImage = zeros(size(I,1)+1, size(I,2)+1, size(I,3), size(I,4), 'single');
    end
    
    intImage(2:end, 2:end, :, :) = cumsum(cumsum(I,1),2);
else
    intImage = [];
end

