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
flags.app = 1;                  % 1 = finance (linear jump diffusion), 2 = tracking (curvilinear dynamics)
flags.gen_data = false;          % true = generate data. false = load real data

if flags.app == 1
    
    % Set finance parameters
    set_parameters_finance;
    set_parameters_finance_gen;
    
elseif flags.app == 2
    
    % Set tracking parameters
    
end

%% Data

if flags.gen_data
    % Generate some data
    [state_x, state_tau, observ, times, interp_x] = generate_data( flags, params );
else
    % Load some data
    convert_data_Apr08;
    MAXK = 1000; params.K = MAXK;
    times = reg_grid(1:MAXK);
    observ = EURUSD_uniform(1:MAXK); interp_x = observ;
end

% Plot data
figure(1), plot(times, interp_x(1,:), 'b', times, observ, 'r'); drawnow;

%% Filtering

% Call filtering algorithm
[ filt_part_sets ] = vr_filter( flags, params, times, observ );

% Calculate jump time kernel density estimate
[filt_kd_times, filt_x_jump_kdest, filt_xdot_jump_kdest] = jump_kernel_est(times(params.K), filt_part_sets{params.K}.pts_tau, filt_part_sets{params.K}.pts_type);

% Plot filtering results
figure(2); clf;
subplot(3,1,1), hold on, plot(times, filt_part_sets{params.K}.pts_intmu(:,:,1)'); ylabel('x'); plot(times, interp_x(1,:), 'b', times, observ, 'r', 'LineWidth', 2);
subplot(3,1,2), hold on, plot(times, filt_part_sets{params.K}.pts_intmu(:,:,2)'); ylabel('x-dot');
    if flags.gen_data, plot(times, interp_x(2,:), 'b', 'LineWidth', 2); end
subplot(3,1,3), hold on, plot(filt_kd_times, filt_x_jump_kdest, 'b'); plot(filt_kd_times, filt_xdot_jump_kdest, 'g'); ylabel('jump kd estimate'); legend('x', 'x-dot');
    if flags.gen_data, for tt=1:length(state_tau), plot([state_tau(tt),state_tau(tt)], [0,1]','r'); end, end
drawnow;

% Histogram number of states
figure(3), hist(filt_part_sets{params.K}.pts_Ns)

%% Smoothing

% Call smoothing algorithm
[ smooth_part_sets] = rb_vr_smoother( flags, params, times, observ, filt_part_sets);

% Calculate jump time kernel density estimate
[smooth_kd_times, smooth_x_jump_kdest, smooth_xdot_jump_kdest] = jump_kernel_est(times(params.K), smooth_part_sets{2}.pts_tau, smooth_part_sets{2}.pts_type);

% Plot smoothing results
figure(4); clf;
subplot(3,1,1), hold on, plot(times, squeeze(smooth_part_sets{1}.pts_intmu(:,1,:))); ylabel('x'); plot(times, interp_x(1,:), 'b', times, observ, 'r', 'LineWidth', 2);
subplot(3,1,2), hold on, plot(times, squeeze(smooth_part_sets{1}.pts_intmu(:,2,:))); ylabel('x-dot');
    if flags.gen_data, plot(times, interp_x(2,:), 'b', 'LineWidth', 2); end
subplot(3,1,3), hold on, plot(smooth_kd_times, smooth_x_jump_kdest, 'b'); plot(smooth_kd_times, smooth_xdot_jump_kdest, 'g'); ylabel('jump kd estimate'); legend('x', 'x-dot');
    if flags.gen_data, for tt=1:length(state_tau), plot([state_tau(tt),state_tau(tt)], [0,1]','r'); end, end
drawnow;

% Histogram number of states
figure(5), hist(smooth_part_sets{1}.pts_Ns)