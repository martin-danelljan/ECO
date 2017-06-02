function hf_out = lhs_operation_joint(hf, samplesf, reg_filter, init_samplef, XH, init_hf, proj_reg)

% This is the left-hand-side operation in Conjugate Gradient

hf_out = cell(size(hf));

% Extract projection matrix and filter separately
P = cellfun(@real, hf(2,1,:), 'uniformoutput',false);
hf = hf(1,1,:);

% size of the padding
num_features = length(hf);
output_sz = [size(hf{1},1), 2*size(hf{1},2)-1];

% Compute the operation corresponding to the data term in the optimization
% (blockwise matrix multiplications)
%implements: A' diag(sample_weights) A f

% sum over all features and feature blocks
sh = sum(bsxfun(@times, samplesf{1}, hf{1}), 3);    % assumes the feature with the highest resolution is first
pad_sz = cell(1,1,num_features);
for k = 2:num_features
    pad_sz{k} = (output_sz - [size(hf{k},1), 2*size(hf{k},2)-1]) / 2;
    
    sh(1+pad_sz{k}(1):end-pad_sz{k}(1), 1+pad_sz{k}(2):end,1,:) = ...
        sh(1+pad_sz{k}(1):end-pad_sz{k}(1), 1+pad_sz{k}(2):end,1,:) + sum(bsxfun(@times, samplesf{k}, hf{k}), 3);
end

% weight all the samples and take conjugate
% sh = bsxfun(@times,sample_weights,sh);
sh = conj(sh);

% multiply with the transpose
hf_out1 = cell(1,1,num_features);
hf_out1{1} = conj(sum(bsxfun(@times, sh, samplesf{1}), 4));
for k = 2:num_features
    hf_out1{k} = conj(sum(bsxfun(@times, sh(1+pad_sz{k}(1):end-pad_sz{k}(1), 1+pad_sz{k}(2):end,1,1,:), samplesf{k}), 4));
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

% Set the filter to the current one instead (test)
% init_hf = cellfun(@(hf) permute(hf, [3 4 1 2]), hf, 'uniformoutput', false);

% B * P
BP_cell = cell(1,1,num_features);
for k = 1:num_features
    BP_cell{k} = sum(bsxfun(@times, reshape(reshape(init_samplef{k}, [], size(init_samplef{k},3)) * P{k}, size(init_samplef{k},1), size(init_samplef{k},2), []), init_hf{k}), 3);
end

BP = BP_cell{1};
for k = 2:num_features
    BP(1+pad_sz{k}(1):end-pad_sz{k}(1), 1+pad_sz{k}(2):end) = ...
        BP(1+pad_sz{k}(1):end-pad_sz{k}(1), 1+pad_sz{k}(2):end) + BP_cell{k};
end

% multiply with the transpose: A^H * BP
hf_out{1,1,1} = hf_out1{1} +  bsxfun(@times, BP, conj(samplesf{1}));

% B^H * BP
fBP = cell(1,1,num_features);
fBP{1} = reshape(bsxfun(@times, conj(init_hf{1}), BP), [], size(init_hf{1},3));

% Compute proj matrix part: B^H * A_m * f
shBP = cell(1,1,num_features);
shBP{1} = reshape(bsxfun(@times, conj(init_hf{1}), sh), [], size(init_hf{1},3));

for k = 2:num_features
    % multiply with the transpose: A^H * BP
    hf_out{1,1,k} = hf_out1{k} +  bsxfun(@times, BP(1+pad_sz{k}(1):end-pad_sz{k}(1), 1+pad_sz{k}(2):end), conj(samplesf{k}));
    
    % B^H * BP
    fBP{k} = reshape(bsxfun(@times, conj(init_hf{k}), BP(1+pad_sz{k}(1):end-pad_sz{k}(1), 1+pad_sz{k}(2):end)), [], size(init_hf{k},3));
    
    % Compute proj matrix part: B^H * A_m * f
    shBP{k} = reshape(bsxfun(@times, conj(init_hf{k}), sh(1+pad_sz{k}(1):end-pad_sz{k}(1), 1+pad_sz{k}(2):end)), [], size(init_hf{k},3));
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