function xf = interpolate_dft(xf, interp1_fs, interp2_fs)

% Performs the implicit interpolation in the Fourier domain of a sample by
% multiplying with the Fourier coefficients of the interpolation function.

xf = cellfun(@(xf, interp1_fs, interp2_fs) bsxfun(@times, bsxfun(@times, xf, interp1_fs), interp2_fs), ...
    xf, interp1_fs, interp2_fs, 'uniformoutput', false);
