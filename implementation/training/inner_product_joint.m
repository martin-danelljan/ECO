function ip = inner_product_joint(xf, yf)

% Computes the joint inner product between two filters and projection matrices.

ip = 0;
for k = 1:size(xf,3)
    % Filter part
    ip = ip + 2*gather(xf{1,1,k}(:)' * yf{1,1,k}(:)) - gather(reshape(xf{1,1,k}(:,end,:), [], 1, 1)' * reshape(yf{1,1,k}(:,end,:), [], 1, 1));
    
    % Projection matrix part
    ip = ip + gather(xf{2,1,k}(:)' * yf{2,1,k}(:));
end
ip = real(ip);