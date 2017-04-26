function xf = compact_fourier_coeff(xf)

% Creates a compact fourier series representation by removing the strict right
% half plane.

if iscell(xf)
    xf = cellfun(@(xf) xf(:,1:(size(xf,2)+1)/2,:), xf, 'uniformoutput', false);
else
    xf = xf(:,1:(size(xf,2)+1)/2,:);
end