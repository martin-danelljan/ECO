function nx = norm_cdcf(xf)

% Computes the norm of the filter.

nx = sqrt(inner_product_cdcf(xf,xf));