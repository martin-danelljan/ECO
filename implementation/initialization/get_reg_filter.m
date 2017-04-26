function reg_filter = get_reg_filter(sz, target_sz, params, reg_window_edge)

% Compute the spatial regularization function and derive the corresponding
% filter operation used for the optimization

if nargin < 3 || isempty(reg_window_edge)
    reg_window_edge = params.reg_window_edge;
end

if params.use_reg_window
    % create weight window
    reg_window_power = params.reg_window_power;
    
    % normalization factor
    reg_scale = 0.5 * target_sz;
    
    % construct grid
    wrg = -(sz(1)-1)/2:(sz(1)-1)/2;
    wcg = -(sz(2)-1)/2:(sz(2)-1)/2;
    [wrs, wcs] = ndgrid(wrg, wcg);
    
    % construct the regukarization window
    reg_window = (reg_window_edge - params.reg_window_min) * (abs(wrs/reg_scale(1)).^reg_window_power + abs(wcs/reg_scale(2)).^reg_window_power) + params.reg_window_min;
    
    % compute the DFT and enforce sparsity
    reg_window_dft = fft2(reg_window) / prod(sz);
    reg_window_dft(abs(reg_window_dft) < params.reg_sparsity_threshold * max(abs(reg_window_dft(:)))) = 0;
    
    % do the inverse transform, correct window minimum
    reg_window_sparse = real(ifft2(reg_window_dft));
    reg_window_dft(1,1) = reg_window_dft(1,1) - prod(sz) * min(reg_window_sparse(:)) + params.reg_window_min;
    reg_window_dft = fftshift(reg_window_dft);
    
    % find the regularization filter by removing the zeros
    reg_filter = single(real(reg_window_dft(~all(reg_window_dft==0,2), ~all(reg_window_dft==0,1))));
else
    % else use a scaled identity matrix
    reg_filter = single(params.reg_window_min);
end