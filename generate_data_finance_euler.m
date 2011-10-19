function [state_x, state_tau, observ, times, interp_x] = generate_data_finance_euler( params )
%GENERATE_DATA_FINANCE Generate jump diffusion finance data

% This version uses a regular time grid.

% Generate x jump times
x_jump_times = 0;
x_jump_mags = 0;
num_x_jumps = 0;
last_x_jump = 0;
while (num_x_jumps==0)||(x_jump_times(num_x_jumps)<params.T)
    num_x_jumps = num_x_jumps + 1;
    x_jump_times(num_x_jumps) = last_x_jump + rande(1/params.x_jump_rate);
    x_jump_mags(num_x_jumps) = normrnd(params.x_jump_mu, params.x_jump_sigma);
    last_x_jump = x_jump_times(num_x_jumps);
end

% Generate xdot jump times
xdot_jump_times = 0;
xdot_jump_mags = 0;
num_xdot_jumps = 0;
last_xdot_jump = 0;
while (num_xdot_jumps==0)||(xdot_jump_times(num_xdot_jumps)<params.T)
    num_xdot_jumps = num_xdot_jumps + 1;
    xdot_jump_times(num_xdot_jumps) = last_xdot_jump + rande(1/params.xdot_jump_rate);
    xdot_jump_mags(num_xdot_jumps) = normrnd(params.xdot_jump_mu, params.xdot_jump_sigma);
    last_xdot_jump = xdot_jump_times(num_xdot_jumps);
end

% Create a state vector array and observation vector array
state = zeros(2,params.K);
observ = zeros(1,params.K);

% Loop through observation times and generate states
last_state = [params.x_start; params.xdot_start];
xi=1;xdi=1; % jump counters
for k=1:params.K
    
    t = k*params.dt;
    
    % See if any jumps have occured
    u = [0; 0];
    if t > x_jump_times(xi)
        u(1) = u(1) + x_jump_mags(xi);
        xi = xi + 1;
    end
    if t > xdot_jump_times(xdi)
        u(2) = u(2) + xdot_jump_mags(xdi);
        xdi = xdi + 1;
    end
    
    % Create transition model matrices
    A = [params.x_decay*params.dt, params.dt; 0 params.xdot_decay*params.dt];
    Q = [params.dt*params.x_sigma^2, 0; 0, params.dt*params.xdot_sigma^2];
    
    % Sample new state
    state(:,k) = mvnrnd(A*last_state+u, Q);
    
    % Sample observation
    observ(1,k) = normrnd(state(1,k), params.obs_sigma);
    
    % Keep it for next time
    last_state = state(:,k);
    
end

state_x = [];
state_tau = sort([x_jump_times, xdot_jump_times], 'ascend');
interp_x = state;
times = (0:params.K-1)*params.dt;

end

