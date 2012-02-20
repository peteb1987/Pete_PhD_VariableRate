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

% Set flags and parameters
flags.gen_data = false;          % true = generate data. false = load real data

% Set finance parameters
set_parameters;
set_parameters_gen;

%% Data

if flags.gen_data
    % Generate some data
    [tau, type, observ, times, interp_state] = generate_data( params );
else
    % Load some data
    convert_data_Apr08;
    MAXK = params.K;
    offset = 0*MAXK;
    times = reg_grid(offset+1:offset+MAXK)-reg_grid(offset+1);
    observ = EURUSD_uniform(offset+1:offset+MAXK);
    interp_state = [];
    tau = [];
    type = [];
end

% Plot data
fig = figure(1);
plot_results(flags, params, times, observ, tau, type, interp_state, [], [], fig);

%% Filtering

% Call filtering algorithm
[ filt_part_sets, filt_weight_sets ] = rb_vr_filter( flags, params, times, observ );

% Calculate jump time kernel density estimate
[filt_kd_est] = jump_kernel_est(length(filt_part_sets{params.K}), times(params.K), cat(2,filt_part_sets{params.K}.tau), cat(2,filt_part_sets{params.K}.type));

% Plot filtering results
f = figure(2);
plot_results(flags, params, times, observ, tau, type, interp_state, filt_part_sets{end}, filt_kd_est, f);

% Histogram number of states
figure(3), hist([filt_part_sets{params.K}.Ns])

%% Smoothing

% Call smoothing algorithm
[ smooth_pts] = rb_vr_smoother( flags, params, times, observ, filt_part_sets, filt_weight_sets);

% Calculate jump time kernel density estimate
[smooth_kd_est] = jump_kernel_est(length(smooth_pts), times(params.K), cat(2,smooth_pts.tau), cat(2,smooth_pts.type));

% Plot smoothing results
f = figure(4);
plot_results(flags, params, times, observ, tau, type, interp_state, smooth_pts, smooth_kd_est, f);

% Histogram number of states
figure(5), hist([smooth_pts.Ns])

%% Kitigawa smoothing

% Create an array of filtering particles
[ filt_pts ] = create_filter_set( params.Np, filt_part_sets, filt_weight_sets );
kiti_pts = kalman_smooth_pts(flags, params, times, filt_part_sets{end});

% Plot Kitagawa results
f = figure(5);
plot_results(flags, params, times, observ, tau, type, interp_state, kiti_pts, filt_kd_est, f);

%% Evaluation

% Errors
[filt_mNs, filt_mospa, filt_mean_rmse, filt_MAP_rmse ] = performance_measures(params, filt_pts, times, observ, tau, interp_state);
[kiti_mNs, kiti_mospa, kiti_mean_rmse, kiti_MAP_rmse ] = performance_measures(params, kiti_pts, times, observ, tau, interp_state);
[VRPS_mNs, VRPS_mospa, VRPS_mean_rmse, VRPS_MAP_rmse ] = performance_measures(params, smooth_pts, times, observ, tau, interp_state);
if ~isempty(interp_state)
    figure(6), hold on
    plot(times, filt_mean_rmse.value_over_time, 'r')
    plot(times, kiti_mean_rmse.value_over_time, 'g')
    plot(times, VRPS_mean_rmse.value_over_time, 'b')
    figure(7), bar([filt_mean_rmse.value, kiti_mean_rmse.value, VRPS_mean_rmse.value]);
    figure(8), hold on
    plot(times, filt_mean_rmse.trend_over_time, 'r')
    plot(times, kiti_mean_rmse.trend_over_time, 'g')
    plot(times, VRPS_mean_rmse.trend_over_time, 'b')
    figure(9), bar([filt_mean_rmse.trend, kiti_mean_rmse.trend, VRPS_mean_rmse.trend]);
end
if ~isempty(interp_state)
    figure(10), hold on
    plot(times, filt_MAP_rmse.value_over_time, 'r')
    plot(times, kiti_MAP_rmse.value_over_time, 'g')
    plot(times, VRPS_MAP_rmse.value_over_time, 'b')
    figure(11), bar([filt_MAP_rmse.value, kiti_MAP_rmse.value, VRPS_MAP_rmse.value]);
    figure(12), hold on
    plot(times, filt_MAP_rmse.trend_over_time, 'r')
    plot(times, kiti_MAP_rmse.trend_over_time, 'g')
    plot(times, VRPS_MAP_rmse.trend_over_time, 'b')
    figure(13), bar([filt_MAP_rmse.trend, kiti_MAP_rmse.trend, VRPS_MAP_rmse.trend]);
end

% Unique particles
kiti_UP = count_unique_particles(times, filt_part_sets{end});
VRPS_UP = count_unique_particles(times, smooth_pts);
figure(14), hold on
plot(times, kiti_UP, 'g');
plot(times, VRPS_UP, 'b');
