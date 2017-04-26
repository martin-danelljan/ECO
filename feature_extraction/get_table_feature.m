function [ features ] = get_table_feature( im, fparam, gparam)
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

%else extract the feature and return it
temp_features = table_lookup(im,tables{tab_ind}.(fparam.tablename));

if gparam.cell_size > 1
    features = average_feature_region(temp_features,gparam.cell_size);
else
    features = temp_features;
end

end

