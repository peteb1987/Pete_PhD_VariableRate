% Test script for variable rate filtering and smoothing algorithms

clup
dbstop if error

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

% Generate some data
[state_x, state_tau, observ, times, interp_x] = generate_data( flags, params );

% Plot data
figure(1), plot(times, interp_x(1,:), 'b', times, observ, 'r'); drawnow;

% Call inference algorithm
[ part_sets ] = vr_inference( flags, params, times, observ );

% Calculate jump time kernel density estimate
[kd_times, jump_kdest] = jump_kernel_est(params.T, part_sets{params.K}.pts_tau);

% Plot results
figure(2);
subplot(3,1,1), hold on, plot(times, part_sets{params.K}.pts_intmu(:,:,1)'); plot(times, interp_x(1,:), 'b', times, observ, 'r', 'LineWidth', 2);
subplot(3,1,2), hold on, plot(times, part_sets{params.K}.pts_intmu(:,:,2)'); plot(times, interp_x(2,:), 'b', 'LineWidth', 2);
subplot(3,1,3), hold on, plot(kd_times, jump_kdest), for tt=1:length(state_tau), plot([state_tau(tt),state_tau(tt)], [0,1]','r'); end
drawnow;

% Histogram number of states
figure(3), hist(part_sets{params.K}.pts_Ns)