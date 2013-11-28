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

%% Filtering

% Call filtering algorithm
[ filt_part_sets, filt_weight_sets ] = rb_vr_filter( flags, params, times, observs );

% Plot filtering results
f1 = figure;
plot_results( flags, params, f1, cp_x, cp_tau, cp_m, times, all_x, observs, filt_part_sets{params.K} );

% RTS smooth
kita_pts = kalman_smooth_pts(flags, params, times, filt_part_sets{end});
f2 = figure;
plot_results( flags, params, f2, cp_x, cp_tau, cp_m, times, all_x, observs, kita_pts );

%% Smoothing

[ smooth_pts ] = rb_vr_smoother( flags, params, times, observs, filt_part_sets, filt_weight_sets);

% Plot smoothing results
f2 = figure;
plot_results( flags, params, f2, cp_x, cp_tau, cp_m, times, all_x, observs, smooth_pts );
