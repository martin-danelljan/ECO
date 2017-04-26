function x = project_sample(x, P)

if ~isempty(P)
    x = cellfun(@(x, P) permute(mtimesx(permute(x, [4 3 1 2]), P, 'speed'), [3 4 2 1]), x, P, 'uniformoutput', false);
end