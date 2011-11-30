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
flags.dyn_mod = 2;              % 1 = tangential and perpendicular accelarations only
                                % 2 = added x and y noise
                                % 3 = added bearing and speed noise
                                % 4 = added bearing and speed jumps
                                % 5 = 3D intrinsic model: random aT, aN and aN-plane
flags.obs_mod = 2;              % 1 = linear gaussian
                                % 2 = radar with gaussian noise

% Observed velocity?
flags.obs_vel = false;           % true = observations of position and velocity. false = only position observations

% Resample-move?
flags.resam_move = true;

% Set tracking parameters
set_parameters_tracking;

%% Data

if flags.gen_data
    % Generate some data
    [state, tau, observ, times, interp_state, ranvar] = generate_data_tracking( flags, params );
else
    % Load some data
    load('D:\pb404\VariableRateAlgorithms\Trajectories\TARGET2.MAT')
end

% Plot data
figure(1);
if flags.dyn_mod>=5
    plot_tracking_data_3D;
else
    plot_tracking_data;
end

%% Filtering

% Call filtering algorithm
[ filt_part_sets ] = vr_filter_tracking( flags, params, times, observ );

% Calculate jump time kernel density estimate
[filt_kd_times, filt_jump_kdest] = jump_kernel_est_tracking(times(params.K), filt_part_sets{params.K});

% Plot filtering results
figure(1);
if flags.dyn_mod>=5
    plot_tracking_data_3D;
    plot_particles_3D
else
    plot_tracking_data;
end

% Histogram number of states
figure(3), hist([filt_part_sets{params.K}.Ns])

%% Smoothing

% Call MCMC smoothing algorithm
[ smooth_pts ] = mcmc_vr_smoother_tracking( flags, params, filt_part_sets, times, observ );

% Calculate jump time kernel density estimate
[smooth_kd_times, smooth_jump_kdest] = jump_kernel_est_tracking(times(params.K), smooth_pts);

% Plot smoothing results
figure(2);
if flags.dyn_mod>=5
    plot_tracking_data_3D;
    plot_smooth_particles_3D;
else
    plot_smoothing_trajectories;
end

% Histogram number of states
figure(4), hist([smooth_pts.Ns])