function bf = cubic_spline_fourier(f, a)

% The continuous Fourier transform of a cubic spline kernel.

bf = -(- 12*a + 12*exp(-pi*f*2i) + 12*exp(pi*f*2i) + 6*a*exp(-pi*f*4i) + ...
    6*a*exp(pi*f*4i) + f.*(pi*exp(-pi*f*2i)*12i) - f.*(pi*exp(pi*f*2i)*12i) + ...
    a*f.*(pi*exp(-pi*f*2i)*16i) - a*f.*(pi*exp(pi*f*2i)*16i) + ...
    a*f.*(pi*exp(-pi*f*4i)*4i) - a*f.*(pi*exp(pi*f*4i)*4i) - 24)./(16*f.^4*pi^4);

bf(f == 0) = 1;