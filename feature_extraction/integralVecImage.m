function intImage = integralVecImage(I)

% Compute the integral image of I.

if ~isempty(I)
    intImage = zeros(size(I,1)+1, size(I,2)+1, size(I,3), size(I,4), 'like', I);
    intImage(2:end, 2:end, :, :) = cumsum(cumsum(I,1),2);
else
    intImage = [];
end

