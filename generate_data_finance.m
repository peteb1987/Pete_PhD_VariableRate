function [state_x, state_tau, observ, times, interp_x] = generate_data_finance( params )
%GENERATE_DATA_FINANCE Generate jump diffusion finance data

% This version uses a regular time grid.

dt = params.dt;
F = params.F;
C = params.C;

% Generate x jump times
x_jump_times = 0;
x_jump_mags = [0 0];
num_x_jumps = 0;
last_x_jump = 0;
while (num_x_jumps==0)||(x_jump_times(num_x_jumps)<params.T)
    num_x_jumps = num_x_jumps + 1;
    x_jump_times(num_x_jumps,1) = last_x_jump + rande(1/params.x_jump_rate);
    x_jump_mags(num_x_jumps,:) = [normrnd(params.x_jump_mu, params.x_jump_sigma), 0];
    last_x_jump = x_jump_times(num_x_jumps);
end

% Generate xdot jump times
xdot_jump_times = 0;
xdot_jump_mags = [0 0];
num_xdot_jumps = 0;
last_xdot_jump = 0;
while (num_xdot_jumps==0)||(xdot_jump_times(num_xdot_jumps)<params.T)
    num_xdot_jumps = num_xdot_jumps + 1;
    xdot_jump_times(num_xdot_jumps,1) = last_xdot_jump + rande(1/params.xdot_jump_rate);
    xdot_jump_mags(num_xdot_jumps,:) = [0 normrnd(params.xdot_jump_mu, params.xdot_jump_sigma)];
    last_xdot_jump = xdot_jump_times(num_xdot_jumps);
end

% Combine the lists
temp = [x_jump_times x_jump_mags; xdot_jump_times xdot_jump_mags];
temp = sortrows(temp, 1);
jump_times = temp(:,1);
jump_mags = temp(:,2:3);

% Create a state vector array and observation vector array
state = zeros(2,params.K);
observ = zeros(1,params.K);

% Loop through observation times and generate states
last_state = [params.x_start; params.xdot_start];
last_t = 0;
ti=1;               % jump counter
for k=1:params.K
    
    t = last_t + dt;
    
    % Iteratively sample forwards to the next jump and through it
    interm_state = last_state;
    interm_t = last_t;
    while jump_times(ti) < t
        
        % Diffusion
        [A, Q] = diffusion_dynamics(F,C,(jump_times(ti)-interm_t));
        interm_state = mvnrnd(A*interm_state, Q)';
        
        % Jump
        interm_state = interm_state + jump_mags(ti,:)';
        
        % Increment jump counter
        ti = ti + 1;
        
    end
    
    %Sample up to the next time point
    [A, Q] = diffusion_dynamics(F,C,t-interm_t);
    state(:,k) = mvnrnd(A*last_state, Q)';
    
    % Sample observation
    observ(1,k) = normrnd(state(1,k), params.obs_sigma);
    
    % Keep it for next time
    last_state = state(:,k);
    last_t = t;
    
end

state_x = [];
state_tau = jump_times;
interp_x = state;
times = (0:params.K-1)*dt;

end

