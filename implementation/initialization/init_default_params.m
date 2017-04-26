function params = init_default_params(params)

% Initialize default parameters

if ~isfield(params, 'interpolation_method')
    params.interpolation_method = 'none';
end
if ~isfield(params, 'interpolation_centering')
    params.interpolation_centering = false;
end
if ~isfield(params, 'interpolation_windowing')
    params.interpolation_windowing = false;
end
if ~isfield(params, 'clamp_position')
    params.clamp_position = false;
end
if ~isfield(params, 'update_projection_matrix')
    params.update_projection_matrix = true;
end
if ~isfield(params, 'proj_init_method')
    params.proj_init_method = 'pca';
end
if ~isfield(params, 'use_detection_sample')
    params.use_detection_sample = true;
end
if ~isfield(params, 'use_projection_matrix')
    params.use_projection_matrix = true;
end
if ~isfield(params, 'use_sample_merge')
    params.use_sample_merge = true;
end
if ~isfield(params, 'CG_use_FR')
    params.CG_use_FR = true;
end
if ~isfield(params, 'CG_standard_alpha')
    params.CG_standard_alpha = true;
end
if ~isfield(params, 'neglect_higher_frequency')
    params.neglect_higher_frequency = false;
end

if params.use_projection_matrix == false
    params.proj_init_method = 'none';
    params.update_projection_matrix = false;
end