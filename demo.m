
% This demo script runs the ECO tracker on the included "Crossing" video.

% Add paths
setup_paths();

% Load video information
video_path = 'sequences/Crossing';
[seq, ground_truth] = load_video_info(video_path);

% Run ECO
results = testing(seq);