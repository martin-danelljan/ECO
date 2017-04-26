function hf = symmetrize_filter(hf)

% Ensure hermetian symmetry.

for k = 1:length(hf)
    dc_ind = (size(hf{k},1) + 1) / 2;
    hf{k}(dc_ind+1:end,end,:) = conj(flipud(hf{k}(1:dc_ind-1,end,:)));
end