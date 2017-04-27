function scale_filter = scale_filter_update(im, pos, base_target_sz, currentScaleFactor, scale_filter, params)

% Update the scale filter.

% Get scale filter features
scales = currentScaleFactor*scale_filter.scaleSizeFactors;
xs = extract_scale_sample(im, pos, base_target_sz, scales, params.scale_model_sz, params.use_mexResize);

is_first_frame = ~isfield(scale_filter, 's_num');

if is_first_frame
    scale_filter.s_num = xs;
else
    scale_filter.s_num = (1 - params.scale_learning_rate) * scale_filter.s_num + params.scale_learning_rate * xs;
end

% Compute projection basis
bigY = scale_filter.s_num;
bigY_den = xs;
if scale_filter.max_scale_dim
    [scale_filter.basis, ~] = qr(bigY, 0);
    [scale_basis_den, ~] = qr(bigY_den, 0);
else
    [U,~,~] = svd(bigY,'econ');
    scale_filter.basis = U(:,1:params.s_num_compressed_dim);
end
scale_filter.basis = scale_filter.basis';

% Compute numerator
sf_proj = fft(feature_projection_scale(scale_filter.s_num, scale_filter.basis, scale_filter.window), [], 2);
scale_filter.sf_num = bsxfun(@times, scale_filter.yf, conj(sf_proj));

% Update denominator
xs = feature_projection_scale(xs, scale_basis_den', scale_filter.window);
xsf = fft(xs,[],2);
new_sf_den = sum(xsf .* conj(xsf),1);
if is_first_frame
    scale_filter.sf_den = new_sf_den;
else
    scale_filter.sf_den = (1 - params.scale_learning_rate) * scale_filter.sf_den + params.scale_learning_rate * new_sf_den;
end
