function projection_matrix = init_projection_matrix(init_sample, compressed_dim, params)

% Initialize the projection matrix.

% Convert the compressed dimensions to a cell
compressed_dim_cell = permute(num2cell(compressed_dim), [2 3 1]);

% Reshape the sample
x = cellfun(@(x) reshape(x, [], size(x,3)), init_sample, 'uniformoutput', false);
x = cellfun(@(x) bsxfun(@minus, x, mean(x, 1)), x, 'uniformoutput', false);

if strcmpi(params.proj_init_method, 'pca')
    [projection_matrix, ~, ~] = cellfun(@(x) svd(x' * x), x, 'uniformoutput', false);
    projection_matrix = cellfun(@(P, dim) cast(P(:,1:dim), 'like', params.data_type), projection_matrix, compressed_dim_cell, 'uniformoutput', false);
elseif strcmpi(params.proj_init_method, 'rand_uni')
    projection_matrix = cellfun(@(x, dim) randn(size(x,2), dim, 'like', params.data_type), x, compressed_dim_cell, 'uniformoutput', false);
    projection_matrix = cellfun(@(P) bsxfun(@rdivide, P, sqrt(sum(P.^2,1))), projection_matrix, 'uniformoutput', false);
elseif strcmpi(params.proj_init_method, 'none')
    projection_matrix = [];
else
    error('Unknown initialization method for the projection matrix: %s', params.proj_init_method);
end