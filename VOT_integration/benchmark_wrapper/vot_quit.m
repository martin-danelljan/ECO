function vot_quit(regions)

fid = fopen('output.txt', 'w');

for i = 1:numel(regions)
    region = regions{i};
    
    if numel(region) == 1
        fprintf(fid, '%f\n', region);
    elseif numel(region) == 4
        fprintf(fid, '%f,%f,%f,%f\n', region(1), region(2), region(3), region(4));
    elseif numel(region) >= 6 && mod(numel(region), 2) == 0
        fprintf(fid, '%f,', region(1:end-1));
        fprintf(fid, '%f\n', region(end));
    else
        error('VOT: Illegal result format');
    end;
    
end;

fclose(fid);