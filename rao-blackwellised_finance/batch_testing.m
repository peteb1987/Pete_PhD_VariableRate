% Batch testing for variable rate filtering and smoothing algorithms

%% Set Up

clup

% Add toolbox folders to path
addpath('ekfukf/','arraylab/','lightspeed/');

% Get test flag
test_flag = str2double(getenv('SGE_TASK_ID'));
% test_flag = 1;

% Set random seed
s = RandStream('mt19937ar', 'seed', test_flag);
RandStream.setDefaultStream(s);

%% Parameters and Data

set_parameters;
set_parameters_gen;

% Generate some data
[tau, type, observ, times, interp_state] = generate_data( params );

%% Filtering

% Call filtering algorithm
[ filt_part_sets, filt_weight_sets ] = rb_vr_filter( flags, params, times, observ );

%% Smoothing

% Call smoothing algorithm
[ smooth_pts] = rb_vr_smoother( flags, params, times, observ, filt_part_sets, filt_weight_sets);

%% Kitigawa smoothing

% Create an array of filtering particles
[ filt_pts ] = create_filter_set( params.Np, filt_part_sets, filt_weight_sets );
kiti_pts = kalman_smooth_pts(flags, params, times, filt_part_sets{end});

%% Evaluation

% Errors
[filt.mNs, filt.mospa, filt.mean_rmse, filt.MAP_rmse ] = performance_measures(params, filt_pts, times, observ, tau, interp_state);
[kiti.mNs, kiti.mospa, kiti.mean_rmse, kiti.MAP_rmse ] = performance_measures(params, kiti_pts, times, observ, tau, interp_state);
[VRPS.mNs, VRPS.mospa, VRPS.mean_rmse, VRPS.MAP_rmse ] = performance_measures(params, smooth_pts, times, observ, tau, interp_state);

% Unique particles
[~, kiti.unique_sequences, kiti.unique_times] = count_unique_particles(times, filt_part_sets{end});
[~, VRPS.unique_sequences, VRPS.unique_times] = count_unique_particles(times, smooth_pts);

%% Save

filename  = ['RBVRPS_results_' num2str(test_flag) '.mat'];
save(filename, 'filt', 'kiti', 'VRPS');