function benchmark_tracker_wrapper(tracker_name, runfile_name, do_cleanup)

if nargin < 3
    do_cleanup = true;
end

% *************************************************************
% VOT: Always call exit command at the end to terminate Matlab!
% *************************************************************
if do_cleanup
    cleanup = onCleanup(@() exit() );
else
    [pathstr, ~, ~] = fileparts(mfilename('fullpath'));
    cd_ind = strfind(pathstr, filesep());
    pathstr = pathstr(1:cd_ind(end)-1);
    cleanup = onCleanup(@() cd(pathstr));
end

try

% *************************************************************
% VOT: Set random seed to a different value every time.
% *************************************************************
RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', sum(clock)));

% **********************************
% VOT: Get initialization data
% **********************************
[images, region] = vot_initialize();

results = cell(length(images), 1);

bb_scale = 1;

% If the provided region is a polygon ...
if numel(region) > 4
    % Init with an axis aligned bounding box with correct area and center
    % coordinate
    cx = mean(region(1:2:end));
    cy = mean(region(2:2:end));
    x1 = min(region(1:2:end));
    x2 = max(region(1:2:end));
    y1 = min(region(2:2:end));
    y2 = max(region(2:2:end));
    A1 = norm(region(1:2) - region(3:4)) * norm(region(3:4) - region(5:6));
    A2 = (x2 - x1) * (y2 - y1);
    s = sqrt(A1/A2);
    w = s * (x2 - x1) + 1;
    h = s * (y2 - y1) + 1;
else
    cx = region(1) + (region(3) - 1)/2;
    cy = region(2) + (region(4) - 1)/2;
    w = region(3);
    h = region(4);
end

init_c = [cx cy];
init_sz = bb_scale * [w h];

im_size = size(imread(images{1}));
im_size = im_size([2 1]);

init_pos = min(max(round(init_c - (init_sz - 1)/2), [1 1]), im_size);
init_sz = min(max(round(init_sz), [1 1]), im_size - init_pos + 1);

seq.s_frames = images;
seq.init_rect = [init_pos, init_sz];

[file_path, file_name_start, file_ext] = fileparts(seq.s_frames{1});
[~, file_name_end, ~] = fileparts(seq.s_frames{end});
seq.path = file_path;

cd_ind = strfind(file_path, filesep());
seq.name = [file_path(cd_ind(end)+1:end), '00'];

% seq.name = 'vot_seq';
seq.ext = file_ext(2:end);
seq.len = length(seq.s_frames);
seq.nz = length(file_name_start);
seq.startFrame = str2num(file_name_start);
seq.endFrame = str2num(file_name_end);

setup_tracker_paths(tracker_name);

disp([runfile_name '(seq, [], []);']);
otb_res = eval([runfile_name '(seq, [], []);']);

num_frames = numel(images);

for frame = 1:num_frames
    bb = otb_res.res(frame,:);
    sz = bb(3:4);
    c = bb(1:2) + (sz - 1)/2;
    new_sz = sz / bb_scale;
    new_tl = c - (new_sz - 1)/2;
    results{frame} = round([new_tl, new_sz]);
end

% **********************************
% VOT: Output the results
% **********************************
vot_quit(results);

catch err
    [wrapper_pathstr, ~, ~] = fileparts(mfilename('fullpath'));
    cd_ind = strfind(wrapper_pathstr, filesep());
    VOT_path = wrapper_pathstr(1:cd_ind(end));
    
    error_report_path = [VOT_path '/error_reports/'];
    if ~exist(error_report_path, 'dir')
        mkdir(error_report_path);
    end
    
    report_file_name = [error_report_path tracker_name '_' runfile_name datestr(now,'_yymmdd_HHMM') '.mat'];
    
    save(report_file_name, 'err')
    
    rethrow(err);
end
