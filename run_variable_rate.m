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
figure, plot(times, interp_x(1,:), 'b', times, observ, 'r');

% Call inference algorithm
[ part_sets ] = vr_inference( flags, params, times, observ );

% Plot results