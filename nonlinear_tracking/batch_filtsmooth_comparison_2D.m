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
flags.space_dim = 2;            % 2 for 2D, 3 for 3D
flags.obs_mod = 2;              % 1 = linear gaussian
                                % 2 = radar with gaussian noise

% Use resample-move?
flags.resam_move = true;

for r = 1:num_seeds;
    
    r
    
    %% Run it
    
    % Set random seed
    s = RandStream('mt19937ar', 'seed', r);
    RandStream.setDefaultStream(s);
    
    flags.dyn_mod = 2;
    flags.obs_vel = false;
    
    % Set tracking parameters
    set_parameters_gen;
    
    % Generate some data
    [true_x, true_tau, observs, times, true_intx, true_w] = generate_data( flags, params );
    
    % Set random seed
    s = RandStream('mt19937ar', 'seed', r);
    RandStream.setDefaultStream(s);
    
    % Filter
    [ filt_part_sets, filt_weight_sets ] = vr_filter( flags, params, times, observs );
    [~, parents] = systematic_resample(exp(filt_weight_sets{end}), params.Np); kita_pts = filt_part_sets{end}(parents);
    [ ~, ~, filt_rmse] = filter_errors( flags, params, filt_part_sets, filt_weight_sets, times, true_tau, true_w, true_intx );
    [ kita_mNs, ~, ~, ~, kita_rmse] = performance_measures( flags, params, kita_pts, times, true_tau, true_w, true_intx );
    [ kita_Nup, kita_Nut, kita_Nujt ] = count_unique_particles( times, kita_pts );
    
    % Smoother
    [ smooth_pts ] = mcmc_vr_smoother( flags, params, filt_part_sets, filt_weight_sets, times, observs );
    [ smooth_mNs, ~, ~, ~, smooth_rmse] = performance_measures( flags, params, smooth_pts, times, true_tau, true_w, true_intx );
    [ smooth_Nup, smooth_Nut, smooth_Nujt ] = count_unique_particles( times, smooth_pts );
    
    % Store
    results.filt_rmse(r) = filt_rmse;
    results.kita_rmse(r) = kita_rmse;
    results.smooth_rmse(r) = smooth_rmse;
    results.kita_Ns_err(r) = kita_mNs - length(true_tau);
    results.smooth_Ns_err(r) = smooth_mNs - length(true_tau);
    results.kita_Nup(r,:) = kita_Nup;
    results.kita_Nujt(r,:) = kita_Nujt;
    results.smooth_Nup(r,:) = smooth_Nup;
    results.smooth_Nujt(r,:) = smooth_Nujt;
    
end

save('test_filtsmooth_comparison_results.mat', 'results');
