% Test script for variable rate filtering and smoothing algorithms

%% Set Up

clup
dbstop if error
% dbstop if warning

% DEFINE RANDOM SEED
rand_seed = 1;

% Set random seed
s = RandStream('mt19937ar', 'seed', rand_seed);
RandStream.setDefaultStream(s);

%% Flags and Parameters

flags.gen_data = true;          % true = generate data. false = load real data

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

if flags.gen_data
    % Generate some data
    [true_x, true_tau, observs, times, true_intx, true_w] = generate_data( flags, params );
else
    % Load some data
    error('You haven''t got any data in the right format yet, chump');
end

% Plot data
f1 = figure;
plot_results( flags, params, f1, true_x, true_tau, times, true_intx, observs, [], [], [] );

%% Filtering

% Call filtering algorithm
[ filt_part_sets, filt_weight_sets ] = vr_filter( flags, params, times, observs );

% Calculate jump time kernel density estimate
[filt_kd_times, filt_jump_kdest] = jump_time_kernel_density(times(params.K), filt_part_sets{params.K});

% Plot filtering results
plot_results( flags, params, f1, true_x, true_tau, times, true_intx, observs, filt_part_sets{params.K}, filt_kd_times, filt_jump_kdest );

% Histogram number of states
figure(3), hist([filt_part_sets{params.K}.Ns])

%% Smoothing

% Call MCMC smoothing algorithm
[ smooth_pts ] = mcmc_vr_smoother( flags, params, filt_part_sets, times, observs );
% [ smooth_pts ] = vr_smoother( flags, params, filt_part_sets, filt_weight_sets, times, observs );

% Calculate jump time kernel density estimate
[smooth_kd_times, smooth_jump_kdest] = jump_time_kernel_density(times(params.K), smooth_pts);

% Plot smoothing results
f2 = figure(2);
plot_results( flags, params, f2, true_x, true_tau, times, true_intx, observs, smooth_pts, smooth_kd_times, smooth_jump_kdest );

% Histogram number of states
figure(4), hist([smooth_pts.Ns])

%% Movies

% % Movie
% fmov = figure(10);
% filter_results_movie( flags, params, fmov, true_x, true_tau, times, true_intx, observs, filt_part_sets );

%% Evaluation

% Evaluate various performance measures
[ kiti_mNs, kiti_mospa, kiti_rmse, kiti_MAP_rmse ] = performance_measures( flags, params, filt_part_sets{end}, times, true_tau, true_intx );
[ VRPS_mNs, VRPS_mospa, VRPS_rmse, VRPS_MAP_rmse ] = performance_measures( flags, params, smooth_pts, times, true_tau, true_intx );

fprintf(1, 'Number of states: True: %u. Kitigawa: %f. VRPS: %f.\n', length(true_tau), kiti_mNs, VRPS_mNs);
fprintf(1, 'MMSE Position errors: Kitigawa: %f. VRPS: %f.\n', kiti_rmse.pos, VRPS_rmse.pos);
fprintf(1, 'MMSE Velocity errors: Kitigawa: %f. VRPS: %f.\n', kiti_rmse.vel, VRPS_rmse.vel);
fprintf(1, 'MAP Position errors: Kitigawa: %f. VRPS: %f.\n', kiti_MAP_rmse.pos, VRPS_MAP_rmse.pos);
fprintf(1, 'MAP Velocity errors: Kitigawa: %f. VRPS: %f.\n', kiti_MAP_rmse.vel, VRPS_MAP_rmse.vel);
fprintf(1, 'Jump Time MOSPA: Kitigawa: %f. VRPS: %f.\n', kiti_mospa, VRPS_mospa);

figure, hold on
plot(times, kiti_rmse.pos_over_time, 'r');
plot(times, VRPS_rmse.pos_over_time, 'b');

figure, hold on
plot(times, kiti_rmse.vel_over_time, 'r');
plot(times, VRPS_rmse.vel_over_time, 'b');

% Count the number of unique particles
[ kiti_Nup ] = count_unique_particles( times, filt_part_sets{end} );
[ VRPS_Nup ] = count_unique_particles( times, smooth_pts );

figure, hold on
plot(times, kiti_Nup, 'r');
plot(times, VRPS_Nup, 'b');