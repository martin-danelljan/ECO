function x = sample_fs(xf, grid_sz)

% Samples the Fourier series

sz = [size(xf,1) size(xf,2)];

if nargin < 2
    x = prod(sz) * cifft2(xf);
else
    if any(grid_sz < sz)
        error('The grid size must be larger than or equal to the signal size')
    end
    tot_pad = grid_sz - sz;
    pad_sz = ceil(tot_pad/2);
    xf_pad = padarray(xf, pad_sz);
    if any(mod(tot_pad,2) == 1)
        % Handle odd padding
        xf_pad = xf_pad(1:end-mod(tot_pad(1),2), 1:end-mod(tot_pad(2),2), :, :);
    end
    x = prod(grid_sz) * cifft2(xf_pad);
end


