function xf = cfft2(x)

% Find the data type
data_type_complex = complex(zeros(1, 'like', x));

% calculate output size
in_sz = size(x);

% if both dimensions are odd
if all(mod(in_sz(1:2), 2) == 1)
    xf = fftshift(fftshift(fft2(x), 1), 2);
else
    out_sz = in_sz;
    out_sz(1:2) = out_sz(1:2) + mod(out_sz(1:2)+1,2);
    
    % allocate
    xf = zeros(out_sz, 'like', data_type_complex);
    
    xf(1:in_sz(1),1:in_sz(2),:,:) = fftshift(fftshift(fft2(x), 1), 2);
    
    if out_sz(1) ~= in_sz(1)
        xf(end,:,:,:) = conj(xf(1,end:-1:1,:,:));
    end
    if out_sz(2) ~= in_sz(2)
        xf(:,end,:,:) = conj(xf(end:-1:1,1,:,:));
    end
end