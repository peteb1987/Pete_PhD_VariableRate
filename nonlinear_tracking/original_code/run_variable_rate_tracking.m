% Test script for variable rate filtering and smoothing algorithms

%% Set Up

clup
dbstop if error
% dbstop if warning

% DEFINE RANDOM SEED
rand_seed = 0;

% Set random seed
s = RandStream('mt19937ar', 'seed', rand_seed);
RandStream.setDefaultStream(s);

%% Flags and Parameters

flags.gen_data = true;          % true = generate data. false = load real data

% Tracking model flags
flags.space_dim = 2;            % 2 for 2D, 3 for 3D
flags.dyn_mod = 2;              % 1 = tangential and normal accelarations only
                                % 2 = plus linear velocities
                                % 3 = Cartesian accelerations
flags.obs_mod = 2;              % 1 = linear gaussian
                                % 2 = radar with gaussian noise

% Observed velocity?
flags.obs_vel = false;           % true = observations of position and velocity. false = only position observations

% Use resample-move?
flags.resam_move = true;

%% Data

if flags.gen_data
    
    % Set tracking parameters
    set_parameters_gen;
    
    % Generate some data
    [true_x, true_tau, observs, times, true_intx, true_w] = generate_data( flags, params );
else
    
    % Set tracking parameters
    set_parameters_bench;
    
    % Load in some data
    [ params, times, observs, true_intx ] = load_benchmark(flags, params);
    true_tau = [];
    true_x = [];
    true_w = [];
end

% Plot data
f1 = figure;
plot_results( flags, params, f1, true_x, true_tau, times, true_intx, observs, [], [], [] );

% %% FUDGE
% 
% flags.dyn_mod = 1;
% set_parameters;

%% Filtering

% Set random seed
s = RandStream('mt19937ar', 'seed', rand_seed);
RandStream.setDefaultStream(s);

tic;
% Call filtering algorithm
[ filt_part_sets, filt_weight_sets ] = vr_filter( flags, params, times, observs );
toc

% Calculate jump time kernel density estimate
[filt_kd_times, filt_jump_kdest] = jump_time_kernel_density(times(params.K), filt_part_sets{params.K});

% Resample final particles
[~, parents] = systematic_resample(exp(filt_weight_sets{end}), params.Np);
kita_pts = filt_part_sets{end}(parents);

% Plot filtering results
plot_results( flags, params, f1, true_x, true_tau, times, true_intx, observs, kita_pts, filt_kd_times, filt_jump_kdest );

% Histogram number of states
figure(3), hist([filt_part_sets{params.K}.Ns])

%% Smoothing

if flags.dyn_mod == 2
    
    tic;
    % Call MCMC smoothing algorithm
    [ smooth_pts ] = mcmc_vr_smoother( flags, params, filt_part_sets, filt_weight_sets, times, observs );
    % [ smooth_pts ] = mcmc_vr_smoother_newsample( flags, params, filt_part_sets, filt_weight_sets, times, observs );
    % [ smooth_pts ] = vr_smoother( flags, params, filt_part_sets, times, observs );
    toc
    
    % Calculate jump time kernel density estimate
    [smooth_kd_times, smooth_jump_kdest] = jump_time_kernel_density(times(params.K), smooth_pts);
    
    % Plot smoothing results
    f2 = figure(2);
    plot_results( flags, params, f2, true_x, true_tau, times, true_intx, observs, smooth_pts, smooth_kd_times, smooth_jump_kdest );
    
    % Histogram number of states
    figure(4), hist([smooth_pts.Ns])
    
else
    smooth_pts = [];
    smooth_kd_times = [];
    smooth_jump_kdest = [];
end

%% Movies

% % Movie
% fmov = figure(10);
% filter_results_movie( flags, params, fmov, true_x, true_tau, times, true_intx, observs, filt_part_sets );

%% Evaluation

% Evaluate various performance measures
[ filt_rmse, filt_MAP_rmse, corr_filt_rmse ] = filter_errors( flags, params, filt_part_sets, filt_weight_sets, times, true_tau, true_w, true_intx );
[ kiti_mNs, kiti_mospa, kiti_rmse, kiti_MAP_rmse, kiti_corr_rmse] = performance_measures( flags, params, kita_pts, times, true_tau, true_w, true_intx );
[ VRPS_mNs, VRPS_mospa, VRPS_rmse, VRPS_MAP_rmse, VRPS_corr_rmse ] = performance_measures( flags, params, smooth_pts, times, true_tau, true_w, true_intx );

% ENEES
filt_ENEES = calc_filter_ENEES(flags, params, times, true_tau, true_w, true_intx, filt_part_sets, filt_weight_sets);
kita_ENEES = calc_smoother_ENEES(flags, params, times, true_tau, true_w, true_intx, kita_pts);
smooth_ENEES = calc_smoother_ENEES(flags, params, times, true_tau, true_w, true_intx, smooth_pts);

if flags.dyn_mod == 2
    filt_rmse = corr_filt_rmse;
    kiti_rmse = kiti_corr_rmse;
    VRPS_rmse = VRPS_corr_rmse;
end

fprintf(1, 'Number of states: True: %u. Kitigawa: %f. VRPS: %f.\n', length(true_tau), kiti_mNs, VRPS_mNs);
fprintf(1, 'MMSE Position errors: Filter: %f. Kitigawa: %f. VRPS: %f.\n', filt_rmse.pos, kiti_rmse.pos, VRPS_rmse.pos);
fprintf(1, 'MMSE Velocity errors: Filter: %f. Kitigawa: %f. VRPS: %f.\n', filt_rmse.vel, kiti_rmse.vel, VRPS_rmse.vel);
fprintf(1, 'MAP Position errors: Filter: %f. Kitigawa: %f. VRPS: %f.\n', filt_MAP_rmse.pos, kiti_MAP_rmse.pos, VRPS_MAP_rmse.pos);
fprintf(1, 'MAP Velocity errors: Filter: %f. Kitigawa: %f. VRPS: %f.\n', filt_MAP_rmse.vel, kiti_MAP_rmse.vel, VRPS_MAP_rmse.vel);
fprintf(1, 'Jump Time MOSPA: Kitigawa: %f. VRPS: %f.\n', kiti_mospa, VRPS_mospa);

figure, hold on
plot(times, filt_rmse.pos_over_time, 'r');
plot(times, kiti_rmse.pos_over_time, 'b');
if ~isempty(VRPS_rmse.pos_over_time)
    plot(times, VRPS_rmse.pos_over_time, 'g');
end

figure, hold on
plot(times, filt_rmse.vel_over_time, 'r');
plot(times, kiti_rmse.vel_over_time, 'b');
if ~isempty(VRPS_rmse.vel_over_time)
    plot(times, VRPS_rmse.vel_over_time, 'g');
end

% Count the number of unique particles
[ kiti_Nup, kiti_Nut, kiti_Nujt ] = count_unique_particles( times, filt_part_sets{end} );
if ~isempty(smooth_pts)
    [ VRPS_Nup, VRPS_Nut, VRPS_Nujt ] = count_unique_particles( times, smooth_pts );
end

figure, hold on
plot(times, kiti_Nup, 'r');
if ~isempty(smooth_pts)
    plot(times, VRPS_Nup, 'b');
end

figure, hold on
plot(times, filt_ENEES, 'r');
plot(times, kita_ENEES, 'b');
plot(times, smooth_ENEES, 'g');