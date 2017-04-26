% Copy this template configuration file to your VOT workspace.
% Enter the full path to the ECO repository root folder.

ECO_repo_path = ########

tracker_label = 'ECO';
tracker_command = generate_matlab_command('benchmark_tracker_wrapper(''ECO'', ''VOT2016_DEEP_settings'', true)', {[ECO_repo_path '/VOT_integration/benchmark_wrapper']});
tracker_interpreter = 'matlab';
tracker_trax = false;