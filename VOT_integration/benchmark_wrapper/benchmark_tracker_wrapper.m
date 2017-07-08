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

seq.format = 'vot';

setup_tracker_paths(tracker_name);

disp([runfile_name '(seq, [], []);']);
res = eval([runfile_name '(seq, [], []);']);

catch err
    [wrapper_pathstr, ~, ~] = fileparts(mfilename('fullpath'));
    cd_ind = strfind(wrapper_pathstr, filesep());
    VOT_path = wrapper_pathstr(1:cd_ind(end));
    
    error_report_path = [VOT_path 'error_reports/'];
    if ~exist(error_report_path, 'dir')
        mkdir(error_report_path);
    end
    
    report_file_name = [error_report_path tracker_name '_' runfile_name datestr(now,'_yymmdd_HHMM') '.mat'];
    
    save(report_file_name, 'err')
    
    rethrow(err);
end
