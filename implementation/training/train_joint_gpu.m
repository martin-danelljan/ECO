function [hf, projection_matrix, res_norms] = train_joint_gpu(hf, projection_matrix, init_samplef, yf, reg_filter, sample_energy, reg_energy, proj_energy, params, init_CG_opts)

% Initial Gauss-Newton optimization of the filter and
% projection matrix.

% Index for the start of the last column of frequencies
lf_ind = cellfun(@(hf) size(hf,1) * (size(hf,2)-1) + 1, hf(1,1,:), 'uniformoutput', false);

% Construct stuff for the proj matrix part
% init_samplef = cellfun(@(x) permute(x, [4 3 1 2]), init_samplef, 'uniformoutput', false);
init_samplef_H = cellfun(@(X) reshape(X, [], size(X,3))', init_samplef, 'uniformoutput', false);

% Construct preconditioner
diag_M = cell(size(hf));
diag_M(1,1,:) = cellfun(@(m, reg_energy) (1-params.precond_reg_param) * bsxfun(@plus, params.precond_data_param * m, (1-params.precond_data_param) * mean(m,3)) + params.precond_reg_param*reg_energy, sample_energy, reg_energy, 'uniformoutput',false);
diag_M(2,1,:) = cellfun(@(m) params.precond_proj_param * (m + params.projection_reg), proj_energy, 'uniformoutput',false);

% Allocate
rhs_samplef = cell(size(hf));
res_norms = [];

for iter = 1:params.init_GN_iter
    % Project sample with new matrix
    init_samplef_proj = project_sample(init_samplef, projection_matrix);
    init_hf = hf(1,1,:);
    
    % Construct the right hand side vector for the filter part
    rhs_samplef(1,1,:) = cellfun(@(xf, yf) bsxfun(@times, conj(xf), yf), init_samplef_proj, yf, 'uniformoutput', false);
    
    % Construct the right hand side vector for the projection matrix part
    fyf = cellfun(@(f, yf) reshape(bsxfun(@times, conj(f), yf), [], size(f,3)), hf(1,1,:), yf, 'uniformoutput', false);
    rhs_samplef(2,1,:) = cellfun(@(P, XH, fyf, fi) (2*real(XH * fyf - XH(:,fi:end) * fyf(fi:end,:)) - params.projection_reg * P), ...
        projection_matrix, init_samplef_H, fyf, lf_ind, 'uniformoutput', false);
    
    % Initialize the projection matrix increment to zero
    hf(2,1,:) = cellfun(@(P) zeros(size(P), 'like', params.data_type), projection_matrix, 'uniformoutput', false);
    
    % do conjugate gradient
    [hf, res_norms_temp] = pcg_ccot(...
        @(x) lhs_operation_joint_gpu(x, init_samplef_proj, reg_filter, init_samplef, init_samplef_H, init_hf, params.projection_reg),...
        rhs_samplef, init_CG_opts, ...
        @(x) diag_precond(x, diag_M), ...
        [], @inner_product_joint, hf);
    
    % Make the filter symmetric (avoid roundoff errors)
    hf(1,1,:) = symmetrize_filter(hf(1,1,:));
    
    % Add to the projection matrix
    projection_matrix = cellfun(@plus, projection_matrix, hf(2,1,:), 'uniformoutput', false);
    
    res_norms = [res_norms; res_norms_temp];
end

res_norms = res_norms/sqrt(inner_product_joint(rhs_samplef,rhs_samplef));

% Extract filter
hf = hf(1,1,:);