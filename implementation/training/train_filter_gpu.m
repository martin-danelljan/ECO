function [hf, res_norms, CG_state] = train_filter_gpu(hf, samplesf, yf, reg_filter, sample_weights, sample_energy, reg_energy, params, CG_opts, CG_state)

% Do Conjugate gradient optimization of the filter.

sample_weights = permute(sample_weights, [2 3 4 1]);

% Construct the right hand side vector
rhs_samplef = cellfun(@(xf) sum(bsxfun(@times, sample_weights, xf), 4), samplesf, 'uniformoutput', false);
rhs_samplef = cellfun(@(xf, yf) bsxfun(@times, conj(xf), yf), rhs_samplef, yf, 'uniformoutput', false);

% Construct preconditioner
diag_M = cellfun(@(m, reg_energy) (1-params.precond_reg_param) * bsxfun(@plus, params.precond_data_param * m, (1-params.precond_data_param) * mean(m,3)) + params.precond_reg_param*reg_energy, sample_energy, reg_energy, 'uniformoutput',false);

% do conjugate gradient
[hf, res_norms, CG_state] = pcg_ccot(...
    @(x) lhs_operation_gpu(x, samplesf, reg_filter, sample_weights),...
    rhs_samplef, CG_opts, ...
    @(x) diag_precond(x, diag_M), ...
    [], @inner_product_filter, ...
    hf, CG_state);

res_norms = res_norms/sqrt(inner_product_filter(rhs_samplef,rhs_samplef));