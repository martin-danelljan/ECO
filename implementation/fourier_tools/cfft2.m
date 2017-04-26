function xf = cfft2(x)

% calculate output size
in_sz = size(x);

% if both dimensions are odd
if all(mod(in_sz(1:2), 2) == 1)
    xf = fftshift(fftshift(fft2(x), 1), 2);
else
    out_sz = in_sz;
    out_sz(1:2) = out_sz(1:2) + mod(out_sz(1:2)+1,2);
    
    % allocate
    xf = complex(zeros(out_sz, 'single'));
    
    xf(1:in_sz(1),1:in_sz(2),:,:) = fftshift(fftshift(fft2(x), 1), 2);
    
    if out_sz(1) ~= in_sz(1)
        xf(end,:,:,:) = conj(fliplr(xf(1,:,:,:)));
    end
    if out_sz(2) ~= in_sz(2)
        xf(:,end,:,:) = conj(flipud(xf(:,1,:,:)));
    end
end