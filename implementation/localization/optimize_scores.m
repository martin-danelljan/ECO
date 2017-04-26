function [disp_row, disp_col, scale_ind] = optimize_scores(scores_fs, iterations)

% Maximizes the continuous convolution response (classification scores).

% Find the size of the output.
[sz1, sz2, num_scales] = size(scores_fs);
output_sz = [sz1 sz2];

% Do the grid search step by finding the maximum in the sampled response
% for each scale.
sampled_scores = sample_fs(scores_fs);
[max_resp_row, max_row] = max(sampled_scores, [], 1);
[init_max_score, max_col] = max(max_resp_row, [], 2);
max_row_perm = permute(max_row, [2 3 1]);
col = max_col(:)';
row = max_row_perm(sub2ind(size(max_row_perm), col, 1:size(sampled_scores,3)));

% Shift and rescale the coordinate system to [-pi, pi]
trans_row = mod(row - 1 + floor((output_sz(1)-1)/2), output_sz(1)) - floor((output_sz(1)-1)/2);
trans_col = mod(col - 1 + floor((output_sz(2)-1)/2), output_sz(2)) - floor((output_sz(2)-1)/2);
init_pos_y = permute(2*pi * trans_row / output_sz(1), [1 3 2]);
init_pos_x = permute(2*pi * trans_col / output_sz(2), [1 3 2]);

% Set the current maximum to the sampled one
max_pos_y = init_pos_y;
max_pos_x = init_pos_x;

% construct grid
ky = -ceil((output_sz(1) - 1)/2) : floor((output_sz(1) - 1)/2);
kx = (-ceil((output_sz(2) - 1)/2) : floor((output_sz(2) - 1)/2))';

% pre-compute complex exponential
exp_iky = exp(bsxfun(@times, 1i * ky, max_pos_y));
exp_ikx = exp(bsxfun(@times, 1i * kx, max_pos_x));

ky2 = ky.*ky;
kx2 = kx.*kx;

iter = 1;
while iter <= iterations
    % Compute gradient
    ky_exp_ky = bsxfun(@times, ky, exp_iky);
    kx_exp_kx = bsxfun(@times, kx, exp_ikx);
    y_resp = mtimesx(exp_iky, scores_fs, 'speed');
    resp_x = mtimesx(scores_fs, exp_ikx, 'speed');
    grad_y = -imag(mtimesx(ky_exp_ky, resp_x, 'speed'));
    grad_x = -imag(mtimesx(y_resp, kx_exp_kx, 'speed'));
    
    % Compute Hessian
    ival = 1i * mtimesx(exp_iky, resp_x, 'speed');
    H_yy = real(-mtimesx(bsxfun(@times, ky2, exp_iky), resp_x, 'speed') + ival);
    H_xx = real(-mtimesx(y_resp, bsxfun(@times, kx2, exp_ikx), 'speed') + ival);
    H_xy = real(-mtimesx(ky_exp_ky, mtimesx(scores_fs, kx_exp_kx, 'speed'), 'speed'));
    det_H = H_yy .* H_xx - H_xy .* H_xy;
    
    % Compute new position using newtons method
    max_pos_y = max_pos_y - (H_xx .* grad_y - H_xy .* grad_x) ./ det_H;
    max_pos_x = max_pos_x - (H_yy .* grad_x - H_xy .* grad_y) ./ det_H;
    
    % Evaluate maximum
    exp_iky = exp(bsxfun(@times, 1i * ky, max_pos_y));
    exp_ikx = exp(bsxfun(@times, 1i * kx, max_pos_x));
    
    iter = iter + 1;
end

% Evaluate the Fourier series at the estimated locations to find the
% corresponding scores.
max_score = real(mtimesx(mtimesx(exp_iky, scores_fs, 'speed'), exp_ikx, 'speed'));

% check for scales that have not increased in score
ind = max_score < init_max_score;
max_score(ind) = init_max_score(ind);
max_pos_y(ind) = init_pos_y(ind);
max_pos_x(ind) = init_pos_x(ind);

% Find the scale with the maximum response
[max_scale_response, scale_ind] = max(max_score(:));

% Scale the coordinate system to output_sz
disp_row = (mod(max_pos_y(1,1,scale_ind) + pi, 2*pi) - pi) / (2*pi) * output_sz(1);
disp_col = (mod(max_pos_x(1,1,scale_ind) + pi, 2*pi) - pi) / (2*pi) * output_sz(2);
end