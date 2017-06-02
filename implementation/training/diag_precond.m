function hf_out = diag_precond(hf, M_diag)

% This is the preconditioner operation in Conjugate Gradient

hf_out = cellfun(@rdivide, hf, M_diag, 'uniformoutput',false);

end