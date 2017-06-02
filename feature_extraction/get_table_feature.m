function [ feature_map ] = get_table_feature( im, fparam, gparam)
%get per-pixel features using a lookup table, if the gparam feature cell
%size is set to something large than one, the resulting data will be
%averaged in cells of the specified size.
%
%tables are loaded dynamically when needed from the lookup_tables folder,
%and stored in a persistent variable

persistent tables;

if isempty(tables)
    tables = {};
end

if isfield(fparam, 'cell_size')
    cell_size = fparam.cell_size;
else
    cell_size = gparam.cell_size;
end

tab_ind = 0;
for k = 1:length(tables)
    if isfield(tables{k}, fparam.tablename)
        tab_ind = k;
        break;
    end
end

if tab_ind == 0
    tables{end+1} = load(['lookup_tables/' fparam.tablename]);
    tab_ind = length(tables);
end

if strcmp(tables{tab_ind}.inputType,'color')
    if size(im,3) ~= 3
        except = MException('cannot get colorfeature from non color image');
        raise(except);
    end
elseif strcmp(tables{tab_ind}.inputType,'gray')
    if size(im,3) == 3
        im_gray = zeros(size(im,1),size(im,2),1,size(im,4),'uint8');
        for k = 1:size(im,4)
            im_gray(:,:,:,k) = rgb2gray(im(:,:,:,k));
        end
        im = im_gray;
    end
end

% Extract features from table
if gparam.use_gpu
    feature_map = gpuArray(table_lookup(im,tables{tab_ind}.(fparam.tablename)));
else
    feature_map = table_lookup(im,tables{tab_ind}.(fparam.tablename));
end

if cell_size > 1
    feature_map = average_feature_region(feature_map, cell_size);
end

end

