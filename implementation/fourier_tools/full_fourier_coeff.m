function xf = full_fourier_coeff(xf)

% Reconstructs the full Fourier series coefficients.

if iscell(xf)
    xf = cellfun(@(xf) cat(2, xf, conj(rot90(xf(:,1:end-1,:), 2))), xf, 'uniformoutput', false);
else
    xf = cat(2, xf, conj(rot90(xf(:,1:end-1,:), 2)));
end