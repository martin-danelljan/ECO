function ip = inner_product_cdcf(xf, yf)

% Computes the inner product between two filters.

if size(xf,1) == 1
    % Only the filter
    ip_cell = cellfun(@(xf, yf) real(2*(xf(:)' * yf(:)) - reshape(xf(:,end,:), [], 1, 1)' * reshape(yf(:,end,:), [], 1, 1)), xf, yf, 'uniformoutput', false');
    ip = sum(cell2mat(ip_cell));
else
    % Joint inner product with the filter and projection matrix
    ip_cell = cellfun(@(xf, yf) real(2*(xf(:)' * yf(:)) - reshape(xf(:,end,:), [], 1, 1)' * reshape(yf(:,end,:), [], 1, 1)), xf(1,1,:), yf(1,1,:), 'uniformoutput', false);
    ip = sum(cell2mat(ip_cell)) + sum(cellfun(@(xf, yf) real(xf(:)' * yf(:)), xf(2,1,:), yf(2,1,:)));
end