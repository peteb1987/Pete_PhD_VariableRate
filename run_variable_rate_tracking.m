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
flags.dyn_mod = 1;              % 1 = tangential and perpendicular accelarations only
                                % 2 = added x and y noise
                                % 3 = added bearing and speed noise
flags.obs_mod = 2;              % 1 = linear gaussian
                                % 2 = radar with gaussian noise

% Observed velocity?
flags.obs_vel = false;           % true = observations of position and velocity. false = only position observations

% Set tracking parameters
set_parameters_tracking;

%% Data

if flags.gen_data
    % Generate some data
    [state, tau, observ, times, interp_state, ranvar] = generate_data_tracking( flags, params );
else
    % Load some data
    error('You don''t have any real tracking data, chump.');
end

% Plot data
figure(1);
plot_tracking_data;

%% Filtering

% Call filtering algorithm
[ filt_part_sets ] = vr_filter_tracking( flags, params, times, observ );

% Calculate jump time kernel density estimate
[filt_kd_times, filt_jump_kdest] = jump_kernel_est_tracking(times(params.K), filt_part_sets{params.K});

% Plot filtering results
figure(1);
plot_tracking_data;

% Histogram number of states
figure(3), hist([filt_part_sets{params.K}.Ns])

%% Smoothing

%%% WRITE THIS BIT
