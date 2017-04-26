function xf = shift_sample(xf, shift, kx, ky)

% Shift a sample in the Fourier domain. The shift should be normalized to
% the range [-pi, pi].

shift_exp_y = cellfun(@(ky) exp((1i * shift(1)) * ky), ky, 'uniformoutput', false);
shift_exp_x = cellfun(@(kx) exp((1i * shift(2)) * kx), kx, 'uniformoutput', false);
xf = cellfun(@(xf, sy, sx) bsxfun(@times, bsxfun(@times, xf, sy), sx), xf, shift_exp_y, shift_exp_x, 'uniformoutput', false);