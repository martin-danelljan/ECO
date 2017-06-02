function hf_out = lhs_operation_joint(hf, samplesf, reg_filter, init_samplef, XH, init_hf, proj_reg)

% This is the left-hand-side operation in Conjugate Gradient

hf_out = cell(size(hf));

% Extract projection matrix and filter separately
P = cellfun(@real, hf(2,1,:), 'uniformoutput',false);
hf = hf(1,1,:);

% Get sizes
num_features = length(hf);
filter_sz = zeros(num_features,2);
for k = 1:num_features
    filter_sz(k,:) = [size(hf{k},1), size(hf{k},2)];
end
[~, k1] = max(filter_sz(:,1));  % Index for the feature block with the largest spatial size
block_inds = 1:num_features;
block_inds(k1) = [];
output_sz = [size(hf{k1},1), 2*size(hf{k1},2)-1];


% Compute the operation corresponding to the data term in the optimization
% (blockwise matrix multiplications)
%implements: A' diag(sample_weights) A f

% sum over all features and feature blocks
pad_sz = cell(1,1,num_features);
sh = mtimesx(samplesf{k1}, permute(hf{k1}, [3 4 1 2]), 'speed');    % assumes the feature with the highest resolution is first
for k = block_inds
    pad_sz{k} = (output_sz - [filter_sz(k,1), 2*filter_sz(k,2)-1]) / 2;
    
    sh(:,1,1+pad_sz{k}(1):end-pad_sz{k}(1), 1+pad_sz{k}(2):end) = ...
        sh(:,1,1+pad_sz{k}(1):end-pad_sz{k}(1), 1+pad_sz{k}(2):end) + mtimesx(samplesf{k}, permute(hf{k}, [3 4 1 2]), 'speed');
end

% weight all the samples
% sh = bsxfun(@times,sample_weights,sh);

% multiply with the transpose
hf_out1 = cell(1,1,num_features);
hf_out1{k1} = permute(conj(mtimesx(sh, 'C', samplesf{k1}, 'speed')), [3 4 2 1]);
for k = block_inds
    hf_out1{k} = permute(conj(mtimesx(sh(:,1,1+pad_sz{k}(1):end-pad_sz{k}(1), 1+pad_sz{k}(2):end), 'C', samplesf{k}, 'speed')), [3 4 2 1]);
end

% compute the operation corresponding to the regularization term (convolve
% each feature dimension with the DFT of w, and the tramsposed operation)
% add the regularization part
% hf_conv = cell(1,1,num_features);
for k = 1:num_features
    reg_pad = min(size(reg_filter{k},2)-1, size(hf{k},2)-1);
    
    % add part needed for convolution
    hf_conv = cat(2, hf{k}, conj(rot90(hf{k}(:, end-reg_pad:end-1, :), 2)));
    
    % do first convolution
    hf_conv = convn(hf_conv, reg_filter{k});
    
    % do final convolution and put toghether result
    hf_out1{k} = hf_out1{k} + convn(hf_conv(:,1:end-reg_pad,:), reg_filter{k}, 'valid');
end

% Stuff related to the projection matrix

% B * P
BP_cell = cell(1,1,num_features);
for k = 1:num_features
    BP_cell{k} = mtimesx(mtimesx(init_samplef{k}, P{k}, 'speed'), init_hf{k}, 'speed');
end

BP = BP_cell{k1};
for k = block_inds
    BP(1,1,1+pad_sz{k}(1):end-pad_sz{k}(1), 1+pad_sz{k}(2):end) = ...
        BP(1,1,1+pad_sz{k}(1):end-pad_sz{k}(1), 1+pad_sz{k}(2):end) + BP_cell{k};
end

% multiply with the transpose: A^H * BP
hf_out{1,1,k1} = hf_out1{k1} +  permute(bsxfun(@times, BP, conj(samplesf{k1})), [3 4 2 1]);

% B^H * BP
fBP = cell(1,1,num_features);
fBP{k1} = reshape(bsxfun(@times, conj(init_hf{k1}), BP), size(init_hf{k1},1), []).';

% Compute proj matrix part: B^H * A_m * f
shBP = cell(1,1,num_features);
shBP{k1} = reshape(bsxfun(@times, conj(init_hf{k1}), sh), size(init_hf{k1},1), []).';

for k = block_inds
    % multiply with the transpose: A^H * BP
    hf_out{1,1,k} = hf_out1{k} +  permute(bsxfun(@times, BP(1,1,1+pad_sz{k}(1):end-pad_sz{k}(1), 1+pad_sz{k}(2):end), conj(samplesf{k})), [3 4 2 1]);
    
    % B^H * BP
    fBP{k} = reshape(bsxfun(@times, conj(init_hf{k}), BP(1,1,1+pad_sz{k}(1):end-pad_sz{k}(1), 1+pad_sz{k}(2):end)), size(init_hf{k},1), []).';
    
    % Compute proj matrix part: B^H * A_m * f
    shBP{k} = reshape(bsxfun(@times, conj(init_hf{k}), sh(1,1,1+pad_sz{k}(1):end-pad_sz{k}(1), 1+pad_sz{k}(2):end)), size(init_hf{k},1), []).';
end

% hf_out2 = cell(1,1,num_features);
for k = 1:num_features
    fi = size(hf{k},1) * (size(hf{k},2)-1) + 1;  % index where the last frequency column starts
    
    % B^H * BP
    hf_out2 = 2*real(XH{k} * fBP{k} - XH{k}(:,fi:end) * fBP{k}(fi:end,:)) + proj_reg * P{k};
    
    % Compute proj matrix part: B^H * A_m * f
    hf_out{2,1,k} = hf_out2 + (2*real(XH{k} * shBP{k} - XH{k}(:,fi:end) * shBP{k}(fi:end,:)));
end
end