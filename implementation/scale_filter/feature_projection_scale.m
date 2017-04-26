function z = feature_projection_scale(x, projection_matrix, cos_window)

% do the windowing of the output
z = bsxfun(@times, cos_window, projection_matrix * x);
end