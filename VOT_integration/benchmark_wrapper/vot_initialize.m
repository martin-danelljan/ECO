function [images, region] = vot_initialize()

    % read the image file paths
    fid = fopen('images.txt','r'); 
    images = textscan(fid, '%s', 'delimiter', '\n');
    fclose(fid);
    
    images = images{1};

    % read the region
    region = dlmread('region.txt');
    region = region(:);
end