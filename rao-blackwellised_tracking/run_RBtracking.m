% Test script for Rao-Blackwellised variable rate filtering and smoothing algorithms

%% Set Up

clup
dbstop if error
% dbstop if warning

% DEFINE RANDOM SEED
rand_seed = 2;

% Set random seed
s = RandStream('mt19937ar', 'seed', rand_seed);
RandStream.setDefaultStream(s);

%% Flags and Parameters

flags.obs_mod = 1;              % 1 = linear gaussian, 2 = radar
flags.obs_vel = false;          % true = observed velocity, false = position only

flags.resam_move = false;

% Set tracking parameters
set_parameters;

%% Data

[all_x, cp_x, cp_tau, cp_m, cp_u, times, observs] = generate_data( flags, params );

figure, hold on, plot(all_x(1,:), all_x(2,:)), plot(cp_x(1,:), cp_x(2,:), '*g'), plot(observs(1,:), observs(2,:), 'r');
figure, hold on, plot(times, magn(all_x(3:4, :))), plot(cp_tau, magn(cp_x(3:4,:)), '*g');
figure, hold on, plot(times, magn(all_x(5:6, :))), plot(cp_tau, magn(cp_x(5:6,:)), '*g');
figure, plot(times, dot(all_x(3:4,:), all_x(5:6,:)))

% % Plot data
% f1 = figure;
% plot_results( flags, params, f1, true_x, true_tau, times, true_intx, observs, [], [], [] );

%% Filtering

% Call filtering algorithm
[ filt_part_sets, filt_weight_sets ] = rb_vr_filter( flags, params, times, observs );

% % Calculate jump time kernel density estimate
% [filt_kd_times, filt_jump_kdest] = jump_time_kernel_density(times(params.K), filt_part_sets{params.K});

% Plot filtering results
f1 = figure;
plot_results( flags, params, f1, cp_x, cp_tau, cp_m, times, all_x, observs, filt_part_sets{params.K} );

% RTS smooth
kita_pts = kalman_smooth_pts(flags, params, times, filt_part_sets{end});
f2 = figure;
plot_results( flags, params, f2, cp_x, cp_tau, cp_m, times, all_x, observs, kita_pts );

% 
% % Histogram number of states
% figure(3), hist([filt_part_sets{params.K}.Ns])
% 
% %% Smoothing
% 
% % Call MCMC smoothing algorithm
% [ smooth_pts ] = mcmc_vr_smoother( flags, params, filt_part_sets, filt_weight_sets, times, observs );
% % [ smooth_pts ] = vr_smoother( flags, params, filt_part_sets, times, observs );
% 
% % Calculate jump time kernel density estimate
% [smooth_kd_times, smooth_jump_kdest] = jump_time_kernel_density(times(params.K), smooth_pts);
% 
% % Plot smoothing results
% f2 = figure(2);
% plot_results( flags, params, f2, true_x, true_tau, times, true_intx, observs, smooth_pts, smooth_kd_times, smooth_jump_kdest );
% 
% % Histogram number of states
% figure(4), hist([smooth_pts.Ns])
% 
% %% Movies
% 
% % % Movie
% % fmov = figure(10);
% % filter_results_movie( flags, params, fmov, true_x, true_tau, times, true_intx, observs, filt_part_sets );
% 
% %% Evaluation
% 
% % Evaluate various performance measures
% [ filt_rmse, filt_MAP_rmse ] = filter_errors( flags, params, filt_part_sets, times, true_tau, true_intx );
% [ kiti_mNs, kiti_mospa, kiti_rmse, kiti_MAP_rmse ] = performance_measures( flags, params, filt_part_sets{end}, times, true_tau, true_intx );
% [ VRPS_mNs, VRPS_mospa, VRPS_rmse, VRPS_MAP_rmse ] = performance_measures( flags, params, smooth_pts, times, true_tau, true_intx );
% 
% fprintf(1, 'Number of states: True: %u. Kitigawa: %f. VRPS: %f.\n', length(true_tau), kiti_mNs, VRPS_mNs);
% fprintf(1, 'MMSE Position errors: Filter: %f. Kitigawa: %f. VRPS: %f.\n', filt_rmse.pos, kiti_rmse.pos, VRPS_rmse.pos);
% fprintf(1, 'MMSE Velocity errors: Filter: %f. Kitigawa: %f. VRPS: %f.\n', filt_rmse.vel, kiti_rmse.vel, VRPS_rmse.vel);
% fprintf(1, 'MAP Position errors: Filter: %f. Kitigawa: %f. VRPS: %f.\n', filt_MAP_rmse.pos, kiti_MAP_rmse.pos, VRPS_MAP_rmse.pos);
% fprintf(1, 'MAP Velocity errors: Filter: %f. Kitigawa: %f. VRPS: %f.\n', filt_MAP_rmse.vel, kiti_MAP_rmse.vel, VRPS_MAP_rmse.vel);
% fprintf(1, 'Jump Time MOSPA: Kitigawa: %f. VRPS: %f.\n', kiti_mospa, VRPS_mospa);
% 
% figure, hold on
% plot(times, filt_rmse.pos_over_time, 'r');
% plot(times, kiti_rmse.pos_over_time, 'b');
% plot(times, VRPS_rmse.pos_over_time, 'g');
% 
% figure, hold on
% plot(times, filt_rmse.vel_over_time, 'r');
% plot(times, kiti_rmse.vel_over_time, 'b');
% plot(times, VRPS_rmse.vel_over_time, 'g');
% 
% % Count the number of unique particles
% [ kiti_Nup, kiti_Nut, kiti_Nujt ] = count_unique_particles( times, filt_part_sets{end} );
% [ VRPS_Nup, VRPS_Nut, VRPS_Nujt ] = count_unique_particles( times, smooth_pts );
% 
% figure, hold on
% plot(times, kiti_Nup, 'r');
% plot(times, VRPS_Nup, 'b');