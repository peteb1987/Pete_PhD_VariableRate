% Batch script for testing models on benchmark trajectories

%% Set Up

clup

% Add toolbox folders to path
addpath('ekfukf/','arraylab/','lightspeed/');

% How many random seed to try?
num_seeds = 10;

% Get test flag
% test_flag = str2double(getenv('SGE_TASK_ID'));
test_flag = 1;

%% Flags and Parameters

% Tracking model flags
flags.space_dim = 3;            % 2 for 2D, 3 for 3D
flags.obs_mod = 2;              % 1 = linear gaussian
                                % 2 = radar with gaussian noise

% Use resample-move?
flags.resam_move = true;

for r = 1:num_seeds;
    
    r
    
    %% Model 1, no velocities
    
    % Set random seed
    s = RandStream('mt19937ar', 'seed', r);
    RandStream.setDefaultStream(s);
    
    flags.dyn_mod = 1;
    flags.obs_vel = false;
    
    % Set parameters
    set_parameters_bench;
    
    [ params, times, observs, true_intx ] = load_benchmark(flags, params);
    true_tau = [];
    true_x = [];
    true_w = [];
    
    % Filter
    [ filt_part_sets, filt_weight_sets ] = vr_filter( flags, params, times, observs );
    [~, parents] = systematic_resample(exp(filt_weight_sets{end}), params.Np); kita_pts = filt_part_sets{end}(parents);
    [ filt_rmse, ~, ~] = filter_errors( flags, params, filt_part_sets, filt_weight_sets, times, true_tau, true_w, true_intx );
    [ kita_mNs, ~, kita_rmse, ~, ~] = performance_measures( flags, params, kita_pts, times, true_tau, true_w, true_intx );

    % Store
    results.M1withoutVel.filt_rmse(r) = filt_rmse;
    results.M1withoutVel.kita_rmse(r) = kita_rmse;
    results.M1withoutVel.num_cps(r) = kita_mNs;

%     %% Model 1, with velocities
%     
%     % Set random seed
%     s = RandStream('mt19937ar', 'seed', r);
%     RandStream.setDefaultStream(s);
%     
%     flags.dyn_mod = 1;
%     flags.obs_vel = true;
%     
%     % Set parameters
%     set_parameters_bench;
%     
%     [ params, times, observs, true_intx ] = load_benchmark(flags, params);
%     true_tau = [];
%     true_x = [];
%     true_w = [];
%     
%     % Filter
%     [ filt_part_sets, filt_weight_sets ] = vr_filter( flags, params, times, observs );
%     [~, parents] = systematic_resample(exp(filt_weight_sets{end}), params.Np); kita_pts = filt_part_sets{end}(parents);
%     [ filt_rmse, ~, ~] = filter_errors( flags, params, filt_part_sets, filt_weight_sets, times, true_tau, true_intx );
%     [ kita_mNs, ~, kita_rmse, ~, ~] = performance_measures( flags, params, kita_pts, times, true_tau, true_intx );% 
%     % Store
%     results.M1withVel.filt_rmse(r) = filt_rmse;
%     results.M1withVel.kita_rmse(r) = kita_rmse;
%     results.M1withVel.num_cps(r) = kita_mNs;

%     %% Model 2, no velocities
%     
%     % Set random seed
%     s = RandStream('mt19937ar', 'seed', r);
%     RandStream.setDefaultStream(s);
%     
%     flags.dyn_mod = 2;
%     flags.obs_vel = false;
%     
%     % Set parameters
%     set_parameters_bench;
%     
%     [ params, times, observs, true_intx ] = load_benchmark(flags, params);
%     true_tau = [];
%     true_x = [];
%     true_w = [];
%     
%     % Filter
%     [ filt_part_sets, filt_weight_sets ] = vr_filter( flags, params, times, observs );
%     [~, parents] = systematic_resample(exp(filt_weight_sets{end}), params.Np); kita_pts = filt_part_sets{end}(parents);
%     [ ~, ~, filt_rmse] = filter_errors( flags, params, filt_part_sets, filt_weight_sets, times, true_tau, true_intx );
%     [ kita_mNs, ~, ~, ~, kita_rmse] = performance_measures( flags, params, kita_pts, times, true_tau, true_intx );% 
%     % Store
%     results.M2withoutVel.filt_rmse(r) = filt_rmse;
%     results.M2withoutVel.kita_rmse(r) = kita_rmse;
%     results.M2withoutVel.num_cps(r) = kita_mNs;

%     %% Model 2, with velocities
%     
%     % Set random seed
%     s = RandStream('mt19937ar', 'seed', r);
%     RandStream.setDefaultStream(s);
%     
%     flags.dyn_mod = 2;
%     flags.obs_vel = true;
%     
%     % Set parameters
%     set_parameters_bench;
%     
%     [ params, times, observs, true_intx ] = load_benchmark(flags, params);
%     true_tau = [];
%     true_x = [];
%     true_w = [];
%     
%     % Filter
%     [ filt_part_sets, filt_weight_sets ] = vr_filter( flags, params, times, observs );
%     [~, parents] = systematic_resample(exp(filt_weight_sets{end}), params.Np); kita_pts = filt_part_sets{end}(parents);
%     [ ~, ~, filt_rmse] = filter_errors( flags, params, filt_part_sets, filt_weight_sets, times, true_tau, true_intx );
%     [ kita_mNs, ~, ~, ~, kita_rmse] = performance_measures( flags, params, kita_pts, times, true_tau, true_intx );% 
%     % Store
%     results.M2withVel.filt_rmse(r) = filt_rmse;
%     results.M2withVel.kita_rmse(r) = kita_rmse;
%     results.M2withVel.num_cps(r) = kita_mNs;

    %% Model 3, no velocities
    
    % Set random seed
    s = RandStream('mt19937ar', 'seed', r);
    RandStream.setDefaultStream(s);
    
    flags.dyn_mod = 3;
    flags.obs_vel = false;
    
    % Set parameters
    set_parameters_bench;
    
    [ params, times, observs, true_intx ] = load_benchmark(flags, params);
    true_tau = [];
    true_x = [];
    true_w = [];
    
    % Filter
    [ filt_part_sets, filt_weight_sets ] = vr_filter( flags, params, times, observs );
    [~, parents] = systematic_resample(exp(filt_weight_sets{end}), params.Np); kita_pts = filt_part_sets{end}(parents);
    [ filt_rmse, ~, ~] = filter_errors( flags, params, filt_part_sets, filt_weight_sets, times, true_tau, true_w, true_intx );
    [ kita_mNs, ~, kita_rmse, ~, ~] = performance_measures( flags, params, kita_pts, times, true_tau, true_w, true_intx );
    
    % Store
    results.M3withoutVel.filt_rmse(r) = filt_rmse;
    results.M3withoutVel.kita_rmse(r) = kita_rmse;
    results.M3withoutVel.num_cps(r) = kita_mNs;
    
end

save('current_3Dmodel_comparison_results.mat', 'results', 'flags', 'params');
