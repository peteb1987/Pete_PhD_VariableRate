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

% Set flags and parameters
flags.app = 1;                  % 1 = finance (linear jump diffusion), 2 = tracking (curvilinear dynamics)
flags.gen_data = true;          % true = generate data. false = load real data

% Set finance parameters
set_parameters;
set_parameters_finance_gen;

%% Data

if flags.gen_data
    % Generate some data
    [tau, observ, times, interp_state] = generate_data( params );
else
    % Load some data
    convert_data_Apr08;
    MAXK = params.K;
    offset = 0*MAXK;
    times = reg_grid(offset+1:offset+MAXK)-reg_grid(offset+1);
    observ = EURUSD_uniform(offset+1:offset+MAXK); interp_state = observ;
end

% Plot data
figure(1), plot(times, interp_state(1,:), 'b', times, observ, 'r'); drawnow;

%% Filtering

% Call filtering algorithm
[ filt_part_sets, filt_weight_sets ] = rb_vr_filter( flags, params, times, observ );

% Calculate jump time kernel density estimate
[filt_kd_est] = jump_kernel_est(length(filt_part_sets{params.K}), times(params.K), cat(2,filt_part_sets{params.K}.tau), cat(2,filt_part_sets{params.K}.type));

% Plot filtering results
f = figure(2);
plot_results(flags, params, times, observ, tau, interp_state, filt_part_sets{end}, filt_kd_est, f);

% Histogram number of states
figure(3), hist([filt_part_sets{params.K}.Ns])

%% Smoothing

% Call smoothing algorithm
[ smooth_pts] = rb_vr_smoother( flags, params, times, observ, filt_part_sets, filt_weight_sets);

% Calculate jump time kernel density estimate
[smooth_kd_est] = jump_kernel_est(length(smooth_pts), times(params.K), cat(2,smooth_pts.tau), cat(2,smooth_pts.type));

% Plot filtering results
f = figure(2);
plot_results(flags, params, times, observ, tau, interp_state, smooth_pts, smooth_kd_est, f);

% Histogram number of states
figure(3), hist([smooth_pts.Ns])