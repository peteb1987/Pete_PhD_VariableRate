% Test script for variable rate filtering and smoothing algorithms

clup
dbstop if error
% dbstop if warning

% DEFINE RANDOM SEED
rand_seed = 1;

% Set random seed
s = RandStream('mt19937ar', 'seed', rand_seed);
RandStream.setDefaultStream(s);

% Set flags and parameters
flags.app = 1;                  % 1 = finance (linear jump diffusion), 2 = tracking (curvilinear dynamics)

if flags.app == 1
    
    % Set finance parameters
    set_parameters_finance;
    set_parameters_finance_gen;
    
elseif flags.app == 2
    
    % Set tracking parameters
    
end

% % Generate some data
% [state_x, state_tau, observ, times, interp_x] = generate_data( flags, params );

% Load some data
convert_data_Apr08;
MAXK = 1000; params.K = MAXK;
times = reg_grid(1:MAXK);
observ = EURUSD_uniform(1:MAXK); interp_x = observ;

% Plot data
figure(1), plot(times, interp_x(1,:), 'b', times, observ, 'r'); drawnow;

% Call filtering algorithm
[ filt_part_sets ] = vr_filter( flags, params, times, observ );

% Calculate jump time kernel density estimate
[filt_kd_times, filt_jump_kdest] = jump_kernel_est(times(params.K), filt_part_sets{params.K}.pts_tau);

% Plot filtering results
figure(2);
subplot(3,1,1), hold on, plot(times, filt_part_sets{params.K}.pts_intmu(:,:,1)'); ylabel('x'); plot(times, interp_x(1,:), 'b', times, observ, 'r', 'LineWidth', 2);
subplot(3,1,2), hold on, plot(times, filt_part_sets{params.K}.pts_intmu(:,:,2)'); ylabel('x-dot'); %plot(times, interp_x(2,:), 'b', 'LineWidth', 2);
subplot(3,1,3), hold on, plot(filt_kd_times, filt_jump_kdest); ylabel('jump kd estimate'); %for tt=1:length(state_tau), plot([state_tau(tt),state_tau(tt)], [0,1]','r'); end
drawnow;

% Histogram number of states
figure(3), hist(filt_part_sets{params.K}.pts_Ns)

% Call smoothing algorithm
[ smooth_part_sets] = rb_vr_smoother( flags, params, times, observ, filt_part_sets);

% Calculate jump time kernel density estimate
[smooth_kd_times, smooth_jump_kdest] = jump_kernel_est(times(params.K), filt_part_sets{params.K}.pts_tau);

% Plot smoothing results
figure(4)
subplot(3,1,1), hold on, plot(times, smooth_part_sets{2}.intmu(:,:,1)'); ylabel('x'); plot(times, interp_x(1,:), 'b', times, observ, 'r', 'LineWidth', 2);
subplot(3,1,2), hold on, plot(times, smooth_part_sets{2}.intmu(:,:,2)'); ylabel('x-dot'); %plot(times, interp_x(2,:), 'b', 'LineWidth', 2);
subplot(3,1,3), hold on, plot(smooth_kd_times, smooth_jump_kdest); ylabel('jump kd estimate'); %for tt=1:length(state_tau), plot([state_tau(tt),state_tau(tt)], [0,1]','r'); end
drawnow;