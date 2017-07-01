function x = cifft2(xf)

if isa(xf, 'gpuArray') || is_octave()
    x = real(ifft2(ifftshift(ifftshift(xf, 1), 2)));
else
    x = ifft2(ifftshift(ifftshift(xf, 1), 2), 'symmetric');
end