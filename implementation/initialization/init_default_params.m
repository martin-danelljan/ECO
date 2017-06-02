function params = init_default_params(params)

% Initialize default parameters

default_params.use_gpu = false;
default_params.gpu_id = [];
default_params.interpolation_method = 'none';
default_params.interpolation_centering = false;
default_params.interpolation_windowing = false;
default_params.clamp_position = false;
default_params.update_projection_matrix = true;
default_params.proj_init_method = 'pca';
default_params.use_detection_sample = true;
default_params.use_projection_matrix = true;
default_params.use_sample_merge = true;
default_params.CG_use_FR = true;
default_params.CG_standard_alpha = true;

def_param_names = fieldnames(default_params);
for k = 1:numel(def_param_names)
    param_name = def_param_names{k};
    if ~isfield(params, param_name)
        params.(param_name) = default_params.(param_name);
    end
end

if params.use_projection_matrix == false
    params.proj_init_method = 'none';
    params.update_projection_matrix = false;
end

params.visualization = params.visualization || params.debug;