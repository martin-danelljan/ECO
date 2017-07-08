function scale_change_factor = scale_filter_track(im, pos, base_target_sz, currentScaleFactor, scale_filter, params)

% Track the scale using the scale filter.

% Get scale filter features
scales = currentScaleFactor*scale_filter.scaleSizeFactors;
xs = extract_scale_sample(im, pos, base_target_sz, scales, params.scale_model_sz, params.use_mexResize);

% Project
xs = feature_projection_scale(xs, scale_filter.basis, scale_filter.window);

% Get scores
xsf = fft(xs, [], 2);
scale_responsef = sum(scale_filter.sf_num .* xsf, 1) ./ (scale_filter.sf_den + params.lambda);
interp_scale_response = ifft(resizeDFT(scale_responsef, params.number_of_interp_scales), 'symmetric');

recovered_scale_index = find(interp_scale_response == max(interp_scale_response(:)), 1);

if params.do_poly_interp
    % Fit a quadratic polynomial to get a refined scale
    % estimate.
    id1 = mod(recovered_scale_index -1 -1,params.number_of_interp_scales)+1;
    id2 = mod(recovered_scale_index +1 -1,params.number_of_interp_scales)+1;
    
    poly_x = [scale_filter.interpScaleFactors(id1), scale_filter.interpScaleFactors(recovered_scale_index), scale_filter.interpScaleFactors(id2)];
    poly_y = [interp_scale_response(id1), interp_scale_response(recovered_scale_index), interp_scale_response(id2)];
    
    poly_A_mat = [poly_x(1)^2, poly_x(1), 1;...
        poly_x(2)^2, poly_x(2), 1;...
        poly_x(3)^2, poly_x(3), 1 ];
    
    poly = poly_A_mat\poly_y';
    
    scale_change_factor = -poly(2)/(2*poly(1));
else
    scale_change_factor = scale_filter.interpScaleFactors(recovered_scale_index);
end
