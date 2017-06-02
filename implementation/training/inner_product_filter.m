function ip = inner_product_filter(xf, yf)

% Computes the inner product between two filters.

ip = 0;
for k = 1:length(xf)
    ip = ip + 2*gather(xf{k}(:)' * yf{k}(:)) - gather(reshape(xf{k}(:,end,:), [], 1, 1)' * reshape(yf{k}(:,end,:), [], 1, 1));
end
ip = real(ip);
