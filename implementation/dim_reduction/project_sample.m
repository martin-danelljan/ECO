function x = project_sample(x, P)

if isempty(P)
    return;
end

if isa(x{1}, 'gpuArray')
    for k = 1:length(x)
        x{k} = reshape(pagefun(@mtimes, reshape(x{k}, [], size(x{k},3), size(x{k},4)), P{k}), size(x{k},1), size(x{k},2), [], size(x{k},4));
    end
else
    for k = 1:length(x)
        x{k} = reshape(mtimesx(reshape(x{k}, [], size(x{k},3), size(x{k},4)), P{k}, 'speed'), size(x{k},1), size(x{k},2), [], size(x{k},4));
    end
end