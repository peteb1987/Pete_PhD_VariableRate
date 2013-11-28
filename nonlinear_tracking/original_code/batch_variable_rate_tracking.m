% Batch test for variable rate filtering and smoothing algorithms

%% Set Up

clup

% Add toolbox folders to path
addpath('ekfukf/','arraylab/','lightspeed/');

% How many random seed to try?
num_seeds = 10;

% Results arrays
results_filter_pts = cell(num_seeds,1);
results_filter_weights = cell(num_seeds,1);
results_smoother_pts = cell(num_seeds,1);
results_filt_MMSE_rmse = cell(num_seeds,1);
results_filt_MAP_rmse = cell(num_seeds,1);

% Get test flag
% test_flag = str2double(getenv('SGE_TASK_ID'));
test_flag = 1;

% Set random seed to determine data generation
s = RandStream('mt19937ar', 'seed', test_flag);
RandStream.setDefaultStream(s);

%% Flags and Parameters

% Tracking model flags
flags.space_dim = 3;            % 2 for 2D, 3 for 3D
flags.dyn_mod = 2;              % 1 = tangential and normal accelarations only
                                % 2 = plus linear velocities
flags.obs_mod = 2;              % 1 = linear gaussian
                                % 2 = radar with gaussian noise

% Observed velocity?
flags.obs_vel = false;           % true = observations of position and velocity. false = only position observations

% Use resample-move?
flags.resam_move = true;

% Set tracking parameters
set_parameters;

%% Data

% Generate some data
[true_x, true_tau, observs, times, true_intx, true_w] = generate_data( flags, params );

% Loop through seeds
for rand_seed = 1:10;
    
    % Set random seed
    s = RandStream('mt19937ar', 'seed', rand_seed);
    RandStream.setDefaultStream(s);
    
    %% Filtering
    
    % Call filtering algorithm
    [ filt_part_sets, filt_weight_sets ] = vr_filter( flags, params, times, observs );
    results_filter_pts{rand_seed} = filt_part_sets{end};
    results_filter_weights{rand_seed} = filt_weight_sets{end};
    
    %% Smoothing
    
    % Call MCMC smoothing algorithm
    [ results_smoother_pts{rand_seed} ] = mcmc_vr_smoother( flags, params, filt_part_sets, times, observs );
    
    %% Errors
    [ filt_rmse, filt_MAP_rmse ] = filter_errors( flags, params, filt_part_sets, times, true_tau, true_intx );
    results_filt_MMSE_rmse{rand_seed} = filt_rmse;
    results_filt_MAP_rmse{rand_seed} = filt_MAP_rmse;
    
end

save(['VRPFS_scen' num2str(test_flag) '.mat'], ...
     'flags', 'params', ...
     'true_x', 'true_tau', 'observs', 'times', 'true_intx', 'true_w', ...
     'results_filt_MMSE_rmse', 'results_filt_MAP_rmse', ...
     'results_filter_pts', ...
     'results_filter_weights', ...
     'results_smoother_pts' );