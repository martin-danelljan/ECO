function results = tracker(params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_train_samples = params.nSamples;
features = params.t_features;

% Set some default parameters
params = init_default_params(params);
if isfield(params, 't_global')
    global_fparams = params.t_global;
else
    global_fparams = [];
end

% Init sequence data
s_frames = params.s_frames;
pos = params.init_pos(:)';
target_sz = params.init_sz(:)';

% Debug and visualization
debug = params.debug;
visualization = params.visualization || debug;

num_frames = numel(s_frames);

max_train_samples = min(max_train_samples, num_frames);

init_target_sz = target_sz;

% Check if color image
im = imread(s_frames{1});
if size(im,3) == 3
    if all(all(im(:,:,1) == im(:,:,2)))
        is_color_image = false;
    else
        is_color_image = true;
    end
else
    is_color_image = false;
end

if size(im,3) > 1 && is_color_image == false
    im = im(:,:,1);
end

% Calculate search area and initial scale factor
search_area = prod(init_target_sz * params.search_area_scale);
if search_area > params.max_image_sample_size
    currentScaleFactor = sqrt(search_area / params.max_image_sample_size);
elseif search_area < params.min_image_sample_size
    currentScaleFactor = sqrt(search_area / params.min_image_sample_size);
else
    currentScaleFactor = 1.0;
end

% target size at the initial scale
base_target_sz = target_sz / currentScaleFactor;

% window size, taking padding into account
switch params.search_area_shape
    case 'proportional'
        img_sample_sz = floor( base_target_sz * params.search_area_scale);     % proportional area, same aspect ratio as the target
    case 'square'
        img_sample_sz = repmat(sqrt(prod(base_target_sz * params.search_area_scale)), 1, 2); % square area, ignores the target aspect ratio
    case 'fix_padding'
        img_sample_sz = base_target_sz + sqrt(prod(base_target_sz * params.search_area_scale) + (base_target_sz(1) - base_target_sz(2))/4) - sum(base_target_sz)/2; % const padding
    case 'custom'
        img_sample_sz = [base_target_sz(1)*2 base_target_sz(2)*2]; % for testing
end

[features, global_fparams, feature_info] = init_features(features, global_fparams, is_color_image, img_sample_sz, 'odd_cells');

% Set feature info
img_support_sz = feature_info.img_support_sz;
feature_sz = feature_info.data_sz;
feature_dim = feature_info.dim;
num_feature_blocks = length(feature_dim);
feature_reg = permute(num2cell(feature_info.penalty), [2 3 1]);

% Get feature specific parameters
feature_params = init_feature_params(features, feature_info);
feature_extract_info = get_feature_extract_info(features);

if params.use_projection_matrix
    compressed_dim = feature_params.compressed_dim;
else
    compressed_dim = feature_dim;
end
compressed_dim_cell = permute(num2cell(compressed_dim), [2 3 1]);

% Size of the extracted feature maps
feature_sz_cell = permute(mat2cell(feature_sz, ones(1,num_feature_blocks), 2), [2 3 1]);

% Number of Fourier coefficients to save for each filter layer. This will
% be an odd number.
filter_sz = feature_sz + mod(feature_sz+1, 2);
filter_sz_cell = permute(mat2cell(filter_sz, ones(1,num_feature_blocks), 2), [2 3 1]);

% The size of the label function DFT. Equal to the maximum filter size.
output_sz = max(filter_sz, [], 1);

% How much each feature block has to be padded to the obtain output_sz
pad_sz = cellfun(@(filter_sz) (output_sz - filter_sz) / 2, filter_sz_cell, 'uniformoutput', false);

% Compute the Fourier series indices and their transposes
ky = cellfun(@(sz) (-ceil((sz(1) - 1)/2) : floor((sz(1) - 1)/2))', filter_sz_cell, 'uniformoutput', false);
kx = cellfun(@(sz) -ceil((sz(2) - 1)/2) : 0, filter_sz_cell, 'uniformoutput', false);

% construct the Gaussian label function using Poisson formula
sig_y = sqrt(prod(floor(base_target_sz))) * params.output_sigma_factor * (output_sz ./ img_support_sz);
yf_y = cellfun(@(ky) single(sqrt(2*pi) * sig_y(1) / output_sz(1) * exp(-2 * (pi * sig_y(1) * ky / output_sz(1)).^2)), ky, 'uniformoutput', false);
yf_x = cellfun(@(kx) single(sqrt(2*pi) * sig_y(2) / output_sz(2) * exp(-2 * (pi * sig_y(2) * kx / output_sz(2)).^2)), kx, 'uniformoutput', false);
yf = cellfun(@(yf_y, yf_x) yf_y * yf_x, yf_y, yf_x, 'uniformoutput', false);

% construct cosine window
cos_window = cellfun(@(sz) single(hann(sz(1)+2)*hann(sz(2)+2)'), feature_sz_cell, 'uniformoutput', false);
cos_window = cellfun(@(cos_window) cos_window(2:end-1,2:end-1), cos_window, 'uniformoutput', false);

% Compute Fourier series of interpolation function
[interp1_fs, interp2_fs] = cellfun(@(sz) get_interp_fourier(sz, params), filter_sz_cell, 'uniformoutput', false);

% Get the reg_window_edge parameter
reg_window_edge = {};
for k = 1:length(features)
    if isfield(features{k}.fparams, 'reg_window_edge')
        reg_window_edge = cat(3, reg_window_edge, permute(num2cell(features{k}.fparams.reg_window_edge(:)), [2 3 1]));
    else
        reg_window_edge = cat(3, reg_window_edge, cell(1, 1, length(features{k}.fparams.nDim)));
    end
end

% Construct spatial regularization filter
reg_filter = cellfun(@(reg_window_edge) get_reg_filter(img_support_sz, base_target_sz, params, reg_window_edge), reg_window_edge, 'uniformoutput', false);

% Compute the energy of the filter (used for preconditioner)
reg_energy = cellfun(@(reg_filter) real(reg_filter(:)' * reg_filter(:)), reg_filter, 'uniformoutput', false);

if params.use_scale_filter
    [nScales, scale_step, scaleFactors, scale_filter, params] = init_scale_filter(params);
else
    % Use the translation filter to estimate the scale.
    nScales = params.number_of_scales;
    scale_step = params.scale_step;
    scale_exp = (-floor((nScales-1)/2):ceil((nScales-1)/2));
    scaleFactors = scale_step .^ scale_exp;
end

if nScales > 0
    %force reasonable scale changes
    min_scale_factor = scale_step ^ ceil(log(max(5 ./ img_support_sz)) / log(scale_step));
    max_scale_factor = scale_step ^ floor(log(min([size(im,1) size(im,2)] ./ base_target_sz)) / log(scale_step));
end

% Set conjugate gradient uptions
init_CG_opts.CG_use_FR = true;
init_CG_opts.tol = 1e-6;
init_CG_opts.CG_standard_alpha = true;
init_CG_opts.debug = params.debug;
CG_opts.CG_use_FR = params.CG_use_FR;
CG_opts.tol = 1e-6;
CG_opts.CG_standard_alpha = params.CG_standard_alpha;
CG_opts.debug = params.debug;

% initialize the projection matrix
rect_position = zeros(num_frames, 4);

time = 0;

% Initialize and allocate
prior_weights = zeros(max_train_samples,1, 'single');
sample_weights = prior_weights;
samplesf = cell(1, 1, num_feature_blocks);
for k = 1:num_feature_blocks
    samplesf{k} = complex(zeros(max_train_samples,compressed_dim(k),filter_sz(k,1),(filter_sz(k,2)+1)/2,'single'));
end
score_matrix = inf(max_train_samples, 'single');

latest_ind = [];
frames_since_last_train = inf;
num_training_samples = 0;
minimum_sample_weight = params.learning_rate*(1-params.learning_rate)^(2*max_train_samples);

res_norms = [];
residuals_pcg = [];

for frame = 1:num_frames
    %load image
    im = imread(s_frames{frame});
    if size(im,3) > 1 && is_color_image == false
        im = im(:,:,1);
    end

    tic();
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Target localization step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Do not estimate translation and scaling on the first frame, since we 
    % just want to initialize the tracker there
    if frame > 1
        old_pos = inf(size(pos));
        iter = 1;
        
        %translation search
        while iter <= params.refinement_iterations && any(old_pos ~= pos)
            % Extract features at multiple resolutions
            sample_pos = round(pos);
            det_sample_pos = sample_pos;
            sample_scale = currentScaleFactor*scaleFactors;
            xt = extract_features(im, sample_pos, sample_scale, features, global_fparams, feature_extract_info);
                        
            % Project sample
            xt_proj = project_sample(xt, projection_matrix);
            
            % Do windowing of features
            xt_proj = cellfun(@(feat_map, cos_window) bsxfun(@times, feat_map, cos_window), xt_proj, cos_window, 'uniformoutput', false);
            
            % Compute the fourier series
            xtf_proj = cellfun(@cfft2, xt_proj, 'uniformoutput', false);
            
            % Interpolate features to the continuous domain
            xtf_proj = interpolate_dft(xtf_proj, interp1_fs, interp2_fs);
            
            % Compute convolution for each feature block in the Fourier domain
            scores_fs_feat = cellfun(@(hf, xf, pad_sz) padarray(sum(bsxfun(@times, hf, xf), 3), pad_sz), hf_full, xtf_proj, pad_sz, 'uniformoutput', false);
            
            % Also sum over all feature blocks.
            % Gives the fourier coefficients of the convolution response.
            scores_fs = permute(sum(cell2mat(scores_fs_feat), 3), [1 2 4 3]);
            
            % Optimize the continuous score function with Newton's method.
            [trans_row, trans_col, scale_ind] = optimize_scores(scores_fs, params.newton_iterations);
            
            % Compute the translation vector in pixel-coordinates and round
            % to the closest integer pixel.
            translation_vec = [trans_row, trans_col] .* (img_support_sz./output_sz) * currentScaleFactor * scaleFactors(scale_ind);
            scale_change_factor = scaleFactors(scale_ind);
            
            % update position
            old_pos = pos;
            pos = sample_pos + translation_vec;
            
            if params.clamp_position
                pos = max([1 1], min([size(im,1) size(im,2)], pos));
            end
            
            % Do scale tracking with the scale filter
            if nScales > 0 && params.use_scale_filter
                scale_change_factor = scale_filter_track(im, pos, base_target_sz, currentScaleFactor, scale_filter, params);
            end 
            
            % Update the scale
            currentScaleFactor = currentScaleFactor * scale_change_factor;
            
            % Adjust to make sure we are not to large or to small
            if currentScaleFactor < min_scale_factor
                currentScaleFactor = min_scale_factor;
            elseif currentScaleFactor > max_scale_factor
                currentScaleFactor = max_scale_factor;
            end
            
            iter = iter + 1;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Model update step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % Extract sample and init projection matrix
    if frame == 1
        % Extract image region for training sample
        sample_pos = round(pos);
        sample_scale = currentScaleFactor;
        xl = extract_features(im, sample_pos, currentScaleFactor, features, global_fparams, feature_extract_info);
        
        % Do windowing of features
        xlw = cellfun(@(feat_map, cos_window) bsxfun(@times, feat_map, cos_window), xl, cos_window, 'uniformoutput', false);
        
        % Compute the fourier series
        xlf = cellfun(@cfft2, xlw, 'uniformoutput', false);
        
        % Interpolate features to the continuous domain
        xlf = interpolate_dft(xlf, interp1_fs, interp2_fs);
        
        % New sample to be added
        xlf = compact_fourier_coeff(xlf);
        
        % Initialize projection matrix
        xl1 = cellfun(@(x) reshape(x, [], size(x,3)), xl, 'uniformoutput', false);
        xl1 = cellfun(@(x) bsxfun(@minus, x, mean(x, 1)), xl1, 'uniformoutput', false);
        
        if strcmpi(params.proj_init_method, 'pca')
            [projection_matrix, ~, ~] = cellfun(@(x) svd(x' * x), xl1, 'uniformoutput', false);
            projection_matrix = cellfun(@(P, dim) single(P(:,1:dim)), projection_matrix, compressed_dim_cell, 'uniformoutput', false);
        elseif strcmpi(params.proj_init_method, 'rand_uni')
            projection_matrix = cellfun(@(x, dim) single(randn(size(x,2), dim)), xl1, compressed_dim_cell, 'uniformoutput', false);
            projection_matrix = cellfun(@(P) bsxfun(@rdivide, P, sqrt(sum(P.^2,1))), projection_matrix, 'uniformoutput', false);
        elseif strcmpi(params.proj_init_method, 'none')
            projection_matrix = [];
        else
            error('Unknown initialization method for the projection matrix: %s', params.proj_init_method);
        end
        clear xl1 xlw
        
        % Shift sample
        shift_samp = 2*pi * (pos - sample_pos) ./ (sample_scale * img_support_sz);
        xlf = shift_sample(xlf, shift_samp, kx, ky);
        
        % Project sample
        xlf_proj = project_sample(xlf, projection_matrix);
    elseif params.learning_rate > 0
        if ~params.use_detection_sample
            % Extract image region for training sample
            sample_pos = round(pos);
            sample_scale = currentScaleFactor;
            xl = extract_features(im, sample_pos, currentScaleFactor, features, global_fparams, feature_extract_info);
            
            % Project sample
            xl_proj = project_sample(xl, projection_matrix);
            
            % Do windowing of features
            xl_proj = cellfun(@(feat_map, cos_window) bsxfun(@times, feat_map, cos_window), xl_proj, cos_window, 'uniformoutput', false);
            
            % Compute the fourier series
            xlf1_proj = cellfun(@cfft2, xl_proj, 'uniformoutput', false);
            
            % Interpolate features to the continuous domain
            xlf1_proj = interpolate_dft(xlf1_proj, interp1_fs, interp2_fs);
            
            % New sample to be added
            xlf_proj = compact_fourier_coeff(xlf1_proj);
        else
            if debug
                % Only for visualization
                xl = cellfun(@(xt) xt(:,:,:,scale_ind), xt, 'uniformoutput', false);
            end
            
            % Use the sample that was used for detection
            sample_scale = sample_scale(scale_ind);
            xlf_proj = cellfun(@(xf) xf(:,1:(size(xf,2)+1)/2,:,scale_ind), xtf_proj, 'uniformoutput', false);
        end
        
        % Shift the sample so that the target is centered
        shift_samp = 2*pi * (pos - sample_pos) ./ (sample_scale * img_support_sz);
        xlf_proj = shift_sample(xlf_proj, shift_samp, kx, ky);
    end
    
    xlf_proj_perm = cellfun(@(xf) permute(xf, [4 3 1 2]), xlf_proj, 'uniformoutput', false);
        
    if params.use_sample_merge
        % Find the distances with existing samples
        dist_vector = find_cluster_distances(samplesf, xlf_proj_perm, num_feature_blocks, num_training_samples, max_train_samples, params);
        
        [merged_sample, new_cluster, merged_cluster_id, new_cluster_id, score_matrix, prior_weights,num_training_samples] = ...
            merge_clusters(samplesf, xlf_proj_perm, dist_vector, score_matrix, prior_weights,...
                           num_training_samples,num_feature_blocks,max_train_samples,minimum_sample_weight,params);
    else
        % Do the traditional adding of a training sample and weight update
        % of C-COT
        [prior_weights, replace_ind] = update_prior_weights(prior_weights, sample_weights, latest_ind, frame, params);
        latest_ind = replace_ind;
        
        merged_cluster_id = 0;
        new_cluster = xlf_proj_perm;
        new_cluster_id = replace_ind;
    end
    
    if frame > 1 && params.learning_rate > 0 || frame == 1 && ~params.update_projection_matrix
        % Insert the new training sample
        for k = 1:num_feature_blocks
            if merged_cluster_id > 0
                samplesf{k}(merged_cluster_id,:,:,:) = merged_sample{k};
            end
            
            if new_cluster_id > 0
                samplesf{k}(new_cluster_id,:,:,:) = new_cluster{k};
            end
        end
    end

    sample_weights = prior_weights;
           
    train_tracker = (frame < params.skip_after_frame) || (frames_since_last_train >= params.train_gap);
    
    if train_tracker     
        % Used for preconditioning
        new_sample_energy = cellfun(@(xlf) abs(xlf .* conj(xlf)), xlf_proj, 'uniformoutput', false);
        
        if frame == 1
            if params.update_projection_matrix
                hf = cell(2,1,num_feature_blocks);
                lf_ind = cellfun(@(sz) sz(1) * (sz(2)-1)/2 + 1, filter_sz_cell, 'uniformoutput', false);
                proj_energy = cellfun(@(P, yf) 2*sum(abs(yf(:)).^2) / sum(feature_dim) * ones(size(P), 'single'), projection_matrix, yf, 'uniformoutput', false);
            else
                hf = cell(1,1,num_feature_blocks);
            end
            % Initialize the filter
            for k = 1:num_feature_blocks
                hf{1,1,k} = complex(zeros([filter_sz(k,1) (filter_sz(k,2)+1)/2 compressed_dim(k)], 'single'));
            end
            
            % Initialize Conjugate Gradient parameters
            CG_opts.maxit = params.init_CG_iter; % Number of initial iterations if projection matrix is not updated
            init_CG_opts.maxit = ceil(params.init_CG_iter / params.init_GN_iter);
            sample_energy = new_sample_energy;
            rhs_samplef = cell(size(hf));
            diag_M = cell(size(hf));
            p = []; rho = []; r_old = [];
        else
            CG_opts.maxit = params.CG_iter;
            
            if params.CG_forgetting_rate == inf || params.learning_rate >= 1
                % CG will be reset
                p = []; rho = []; r_old = [];
            else
                rho = rho / (1-params.learning_rate)^params.CG_forgetting_rate;
            end
            % Update the approximate average sample energy using the learning
            % rate. This is only used to construct the preconditioner.
            sample_energy = cellfun(@(se, nse) (1 - params.learning_rate) * se + params.learning_rate * nse, sample_energy, new_sample_energy, 'uniformoutput', false);
        end
        
        % Do training
        if frame == 1 && params.update_projection_matrix
            % Initial Gauss-Newton optimization of the filter and
            % projection matrix.
            
            % Construct stuff for the proj matrix part
            init_samplef = cellfun(@(x) permute(x, [4 3 1 2]), xlf, 'uniformoutput', false);
            init_samplef_H = cellfun(@(X) conj(reshape(X, size(X,2), [])), init_samplef, 'uniformoutput', false);
           
            % Construct preconditioner
            diag_M(1,1,:) = cellfun(@(m, reg_energy) (1-params.precond_reg_param) * bsxfun(@plus, params.precond_data_param * m, (1-params.precond_data_param) * mean(m,3)) + params.precond_reg_param*reg_energy, sample_energy, reg_energy, 'uniformoutput',false);
            diag_M(2,1,:) = cellfun(@(m) params.precond_proj_param * (m + params.projection_reg), proj_energy, 'uniformoutput',false);
            
            projection_matrix_init = projection_matrix;
            
            for iter = 1:params.init_GN_iter
                % Project sample with new matrix
                init_samplef_proj = cellfun(@(x,P) mtimesx(x, P, 'speed'), init_samplef, projection_matrix, 'uniformoutput', false);
                init_hf = cellfun(@(x) permute(x, [3 4 1 2]), hf(1,1,:), 'uniformoutput', false);
                
                % Construct the right hand side vector for the filter part
                rhs_samplef(1,1,:) = cellfun(@(xf, yf) bsxfun(@times, conj(permute(xf, [3 4 2 1])), yf), init_samplef_proj, yf, 'uniformoutput', false);
                
                % Construct the right hand side vector for the projection matrix part
                fyf = cellfun(@(f, yf) reshape(bsxfun(@times, conj(f), yf), [], size(f,3)), hf(1,1,:), yf, 'uniformoutput', false);
                rhs_samplef(2,1,:) = cellfun(@(P, XH, fyf, fi) (2*real(XH * fyf - XH(:,fi:end) * fyf(fi:end,:)) - params.projection_reg * P), ...
                    projection_matrix, init_samplef_H, fyf, lf_ind, 'uniformoutput', false);
                
                % Initialize the projection matrix increment to zero
                hf(2,1,:) = cellfun(@(P) zeros(size(P), 'single'), projection_matrix, 'uniformoutput', false);
                
                % do conjugate gradient
                [hf, ~, ~, ~, res_norms_temp] = pcg_ccot(...
                    @(x) lhs_operation_joint(x, init_samplef_proj, reg_filter, feature_reg, init_samplef, init_samplef_H, init_hf, params.projection_reg),...
                    rhs_samplef, init_CG_opts, ...
                    @(x) diag_precond(x, diag_M), ...
                    [], hf);
                
                % Make the filter symmetric (avoid roundoff errors)
                hf(1,1,:) = symmetrize_filter(hf(1,1,:));
                
                % Add to the projection matrix
                projection_matrix = cellfun(@plus, projection_matrix, hf(2,1,:), 'uniformoutput', false);
                
                res_norms = [res_norms; res_norms_temp];
            end
            
            % Extract filter
            hf = hf(1,1,:);
            
            % Re-project and insert training sample
            xlf_proj = project_sample(xlf, projection_matrix);
            for k = 1:num_feature_blocks
                samplesf{k}(1,:,:,:) = permute(xlf_proj{k}, [4 3 1 2]);
            end
            
            if debug
                norm_proj_mat_init = sqrt(sum(cellfun(@(P) norm(P(:))^2, projection_matrix_init)));
                norm_proj_mat = sqrt(sum(cellfun(@(P) norm(P(:))^2, projection_matrix)));
                norm_proj_mat_change = sqrt(sum(cellfun(@(P,P2) norm(P(:) - P2(:))^2, projection_matrix_init, projection_matrix)));
                fprintf('Norm init: %f, Norm final: %f, Matrix change: %f\n', norm_proj_mat_init, norm_proj_mat, norm_proj_mat_change / norm_proj_mat_init);
            end
        else
            % Construct the right hand side vector
            rhs_samplef = cellfun(@(xf) permute(mtimesx(sample_weights, 'T', xf, 'speed'), [3 4 2 1]), samplesf, 'uniformoutput', false);
            rhs_samplef = cellfun(@(xf, yf) bsxfun(@times, conj(xf), yf), rhs_samplef, yf, 'uniformoutput', false);
            
            % Construct preconditioner
            diag_M = cellfun(@(m, reg_energy) (1-params.precond_reg_param) * bsxfun(@plus, params.precond_data_param * m, (1-params.precond_data_param) * mean(m,3)) + params.precond_reg_param*reg_energy, sample_energy, reg_energy, 'uniformoutput',false);
            
            % do conjugate gradient
            [hf, ~, ~, ~, res_norms, p, rho, r_old] = pcg_ccot(...
                @(x) lhs_operation(x, samplesf, reg_filter, sample_weights, feature_reg),...
                rhs_samplef, CG_opts, ...
                @(x) diag_precond(x, diag_M), ...
                [], hf, p, rho, r_old);
        end
        
        % Reconstruct the full Fourier series
        hf_full = full_fourier_coeff(hf);
        
        frames_since_last_train = 0;
    else
        frames_since_last_train = frames_since_last_train+1;
    end
    
    % Update the scale filter
    if nScales > 0 && params.use_scale_filter
        scale_filter = scale_filter_update(im, pos, base_target_sz, currentScaleFactor, scale_filter, params);
    end
    
    % Update the target size (only used for computing output box)
    target_sz = base_target_sz * currentScaleFactor;
    
    %save position and calculate FPS
    rect_position(frame,:) = round([pos([2,1]) - (target_sz([2,1]) - 1)/2, target_sz([2,1])]);
    
    time = time + toc();
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Visualization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % debug visualization
    if debug
        figure(20)
%         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        subplot_cols = num_feature_blocks;
        subplot_rows = 3;%ceil(feature_dim/subplot_cols);
        for disp_layer = 1:num_feature_blocks;
            subplot(subplot_rows,subplot_cols,disp_layer);
            imagesc(prod(filter_sz(disp_layer,:)) * mean(abs(ifft2(conj(ifftshift(hf_full{disp_layer})), 'symmetric')), 3)); 
            colorbar;
            axis image;
            subplot(subplot_rows,subplot_cols,disp_layer+subplot_cols);
            imagesc(mean(abs(xl{disp_layer}), 3)); 
            colorbar;
            axis image;
            if frame > 1
                subplot(subplot_rows,subplot_cols,disp_layer+2*subplot_cols);
                imagesc(prod(output_sz) * fftshift(ifft2(ifftshift(scores_fs_feat{disp_layer}(:,:,1,scale_ind)), 'symmetric')));
                colorbar;
                axis image;
            end
        end
        
        if train_tracker
            residuals_pcg = [residuals_pcg; res_norms/norm_cdcf(rhs_samplef)];
            res_start_ind = max(1, length(residuals_pcg)-300);
            figure(99);plot(res_start_ind:length(residuals_pcg), residuals_pcg(res_start_ind:end));
            axis([res_start_ind, length(residuals_pcg), 0, min(max(residuals_pcg(res_start_ind:end)), 0.2)]);
        end
    end
    
    % visualization
    if visualization == 1
        rect_position_vis = [pos([2,1]) - (target_sz([2,1]) - 1)/2, target_sz([2,1])];
        im_to_show = double(im)/255;
        if size(im_to_show,3) == 1
            im_to_show = repmat(im_to_show, [1 1 3]);
        end
        if frame == 1,  %first frame, create GUI
            fig_handle = figure('Name', 'Tracking');
%             set(fig_handle, 'Position', [100, 100, size(im,2), size(im,1)]);
            imagesc(im_to_show);
            hold on;
            rectangle('Position',rect_position_vis, 'EdgeColor','g', 'LineWidth',2);
            text(10, 10, int2str(frame), 'color', [0 1 1]);
            hold off;
            axis off;axis image;set(gca, 'Units', 'normalized', 'Position', [0 0 1 1])
            
%             output_name = 'Video_name';
%             opengl software;
%             writer = VideoWriter(output_name, 'MPEG-4');
%             writer.FrameRate = 5;
%             open(writer);
        else
            % Do visualization of the sampled confidence scores overlayed
            resp_sz = round(img_support_sz*currentScaleFactor*scaleFactors(scale_ind));
            xs = floor(det_sample_pos(2)) + (1:resp_sz(2)) - floor(resp_sz(2)/2);
            ys = floor(det_sample_pos(1)) + (1:resp_sz(1)) - floor(resp_sz(1)/2);
            
            % To visualize the continuous scores, sample them 10 times more
            % dense than output_sz. 
            sampled_scores_display = fftshift(sample_fs(scores_fs(:,:,scale_ind), 10*output_sz));
            
            figure(fig_handle);
%                 set(fig_handle, 'Position', [100, 100, 100+size(im,2), 100+size(im,1)]);
            imagesc(im_to_show);
            hold on;
            resp_handle = imagesc(xs, ys, sampled_scores_display); colormap hsv;
            alpha(resp_handle, 0.5);
            rectangle('Position',rect_position_vis, 'EdgeColor','g', 'LineWidth',2);
            text(10, 10, int2str(frame), 'color', [0 1 1]);
            hold off;
            
%                 axis off;axis image;set(gca, 'Units', 'normalized', 'Position', [0 0 1 1])
        end
        
        drawnow
%         if frame > 1
%             if frame < inf
%                 writeVideo(writer, getframe(gcf));
%             else
%                 close(writer);
%             end
%         end
%          pause
    end
end

% close(writer);

fps = numel(s_frames) / time;

disp(['fps: ' num2str(fps)])

results.type = 'rect';
results.res = rect_position;%each row is a rectangle
results.fps = fps;

