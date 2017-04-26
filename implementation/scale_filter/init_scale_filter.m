function [nScales, scale_step, scaleFactors, scale_filter, params] = init_scale_filter(params)

% Initialize the scale filter parameters. Uses the fDSST scale filter.

init_target_sz = params.init_sz(:)';

nScales = params.number_of_scales_filter;
scale_step = params.scale_step_filter;

scale_sigma = params.number_of_interp_scales * params.scale_sigma_factor;

scale_exp = (-floor((nScales-1)/2):ceil((nScales-1)/2)) * params.number_of_interp_scales/nScales;
scale_exp_shift = circshift(scale_exp, [0 -floor((nScales-1)/2)]);

interp_scale_exp = -floor((params.number_of_interp_scales-1)/2):ceil((params.number_of_interp_scales-1)/2);
interp_scale_exp_shift = circshift(interp_scale_exp, [0 -floor((params.number_of_interp_scales-1)/2)]);

scale_filter.scaleSizeFactors = scale_step .^ scale_exp;
scale_filter.interpScaleFactors = scale_step .^ interp_scale_exp_shift;

ys = exp(-0.5 * (scale_exp_shift.^2) /scale_sigma^2);
scale_filter.yf = single(fft(ys));
scale_filter.window = single(hann(size(ys,2)))';

%make sure the scale model is not to large, to save computation time
if params.scale_model_factor^2 * prod(init_target_sz) > params.scale_model_max_area
    params.scale_model_factor = sqrt(params.scale_model_max_area/prod(init_target_sz));
end

%set the scale model size
params.scale_model_sz = max(floor(init_target_sz * params.scale_model_factor), [8 8]);

scale_filter.max_scale_dim = strcmp(params.s_num_compressed_dim,'MAX');
if scale_filter.max_scale_dim
    params.s_num_compressed_dim = length(scale_filter.scaleSizeFactors);
end

% Scale factor for the translation filter
scaleFactors = 1;
