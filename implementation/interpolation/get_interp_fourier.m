function [interp1_fs, interp2_fs] = get_interp_fourier(sz, params)

% Compute the Fourier series of the interpolation function (b_d in the
% paper). The interpolation method is set using params.interpolation_method:
% - 'ideal' performs ideal reconstruction, which corresponds to a periodic
%   summation of a sinc function b_d in the spatial domain.
% - 'bicubic' uses a bicubic kernel b_d in the spatial domain.

switch lower(params.interpolation_method)
    case 'none'
        % Do nothing
        interp1_fs = ones(sz(1),1);
        interp2_fs = ones(1,sz(2));
    case 'ideal'
        % Ideal reconstruction (flat in the frequency domain)
        interp1_fs = ones(sz(1),1) / sz(1);
        interp2_fs = ones(1,sz(2)) / sz(2);
    case 'bicubic'
        % Take the truncated fourier series from the cubic spline
        a = params.interpolation_bicubic_a;
        interp1_fs = real(1/sz(1) * cubic_spline_fourier((-(sz(1)-1)/2:(sz(1)-1)/2)'/sz(1), a));
        interp2_fs = real(1/sz(2) * cubic_spline_fourier((-(sz(2)-1)/2:(sz(2)-1)/2)/sz(2), a));
    otherwise
        error('Unknown dft interpolation method');
end

if params.interpolation_centering
    % Center the feature grids by shifting the interpolated features
    % Multiply Fourier coeff with e^(-i*pi*k/N)
    interp1_fs = interp1_fs .* exp(-1i*pi / sz(1) * (-(sz(1)-1)/2:(sz(1)-1)/2)');
    interp2_fs = interp2_fs .* exp(-1i*pi / sz(2) * (-(sz(2)-1)/2:(sz(2)-1)/2));
end
if params.interpolation_windowing
    % Window the Fourier series of the interpolation basis function
    win1 = hann(sz(1)+2);
    win2 = hann(sz(2)+2);
    interp1_fs = interp1_fs .* win1(2:end-1);
    interp2_fs = interp2_fs .* win2(2:end-1)';
end

interp1_fs = single(interp1_fs);
interp2_fs = single(interp2_fs);