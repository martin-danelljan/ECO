function [seq, ground_truth] = load_video_info(video_path)

ground_truth = dlmread([video_path '/groundtruth_rect.txt']);

seq.format = 'otb';
seq.len = size(ground_truth, 1);
seq.init_rect = ground_truth(1,:);

img_path = [video_path '/img/'];

if exist([img_path num2str(1, '%04i.png')], 'file'),
    img_files = num2str((1:seq.len)', [img_path '%04i.png']);
elseif exist([img_path num2str(1, '%04i.jpg')], 'file'),
    img_files = num2str((1:seq.len)', [img_path '%04i.jpg']);
elseif exist([img_path num2str(1, '%04i.bmp')], 'file'),
    img_files = num2str((1:seq.len)', [img_path '%04i.bmp']);
else
    error('No image files to load.')
end

seq.s_frames = cellstr(img_files);

end

