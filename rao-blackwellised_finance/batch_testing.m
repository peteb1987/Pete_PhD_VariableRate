% Batch testing for variable rate filtering and smoothing algorithms

%% Set Up

clup

% Add toolbox folders to path
addpath('ekfukf/','arraylab/','lightspeed/');

% Get test flag
% test_flag = str2double(getenv('SGE_TASK_ID'));
test_flag = 1;

% How many random seed to try?
num_seeds = 10;

for r = 1:num_seeds;
    
    r
    
    % Set random seed
    s = RandStream('mt19937ar', 'seed', r);
    RandStream.setDefaultStream(s);
    
    % Parameters
    set_parameters;
    set_parameters_gen;
    
    % Generate some data
    [tau, type, observ, times, interp_state] = generate_data( params );
    
    % Set random seed
    s = RandStream('mt19937ar', 'seed', r);
    RandStream.setDefaultStream(s);
    
    % Call filtering algorithm
    [ filt_part_sets, filt_weight_sets ] = rb_vr_filter( flags, params, times, observ );
    
    % Call smoothing algorithm
    [ smooth_pts] = rb_vr_smoother( flags, params, times, observ, filt_part_sets, filt_weight_sets);
    
    % Create an array of filtering particles
    [ filt_pts ] = create_filter_set( params.Np, filt_part_sets, filt_weight_sets );
    kita_pts = kalman_smooth_pts(flags, params, times, filt_part_sets{end});
    
    % Errors
    [filt.mNs, filt.mospa, filt.mean_rmse, filt.MAP_rmse ] = performance_measures(params, filt_pts, times, observ, tau, interp_state);
    [kita.mNs, kita.mospa, kita.mean_rmse, kita.MAP_rmse ] = performance_measures(params, kita_pts, times, observ, tau, interp_state);
    [VRPS.mNs, VRPS.mospa, VRPS.mean_rmse, VRPS.MAP_rmse ] = performance_measures(params, smooth_pts, times, observ, tau, interp_state);
    
    % Unique particles
    [kita.unique_over_time, kita.unique_sequences, kita.unique_times] = count_unique_particles(times, filt_part_sets{end});
    [VRPS.unique_over_time, VRPS.unique_sequences, VRPS.unique_times] = count_unique_particles(times, smooth_pts);
    
    % Store
    results.filt(r) = filt;
    results.kita(r) = kita;
    results.VRPS(r) = VRPS;
    
end

%% Save

filename  = ['results_' num2str(test_flag) '.mat'];
save(filename, 'results');