function [state, tau, observ, times, interp_state, ranvar] = generate_data_tracking( flags, params )
%GENERATE_DATA_FINANCE Generate curvilinear tracking data

% This version uses a regular time grid.

dt = params.dt;
Q = params.Q;
R = params.R;
ds = params.state_dim; do = params.obs_dim; dr = params.rnd_dim;

% Generate state times
tau = 0;
Ns = 1;
last_tau = 0;
while (Ns==0)||(tau(Ns)<params.T)
    Ns = Ns + 1;
    tau(Ns,1) = last_tau + gamrnd(params.rate_shape, params.rate_scale);
    last_tau = tau(Ns,1);
end

% Create arrays
times = cumsum(dt*ones(params.K,1))-dt;     % Array of times
state = zeros(ds,Ns);                       % Array of states
ranvar = zeros(dr,Ns);                      % Array of random variables
interp_state = zeros(ds,params.K);          % Array of interpolated states
observ = zeros(do,params.K);                % Array of observations

% Initialise for loop
ti=2;                               % jump counter
w = mvnrnd(zeros(1,dr), Q)';        % starting r.v.s
x = params.start_state;             % starting state
last_x = x;
last_tau = 0;

% Store first jump values
state(:, 1) = x;
ranvar(:, 1) = w;
mu = tracking_calc_obs_mean(flags, params, x);
observ(:,1) = mvnrnd(mu', R)';
interp_state(:,1) = x;

% Loop through observation times and generate states
for k=2:params.K
    
    t = times(k);
    
    % Iteratively sample forwards to the next jump and through it
    while tau(ti) < t
        
        % New state
        x = tracking_calc_next_state(flags, last_x, tau(ti)-last_tau, w);
        w = mvnrnd(zeros(1,dr), Q)';
        
        % Store state
        state(:,ti) = x;
        ranvar(:,ti) = w;
        
        % Update last_ variables
        last_tau = tau(ti);
        last_x = x;
        
        % Increment jump counter
        ti = ti + 1;
        
    end
    
    % Interpolate state
    x = tracking_calc_next_state(flags, last_x, t-last_tau, w);
    
    % Store it
    interp_state(:,k) = x;
    
    % Sample observation
    mu = tracking_calc_obs_mean(flags, params, x);
    observ(:,k) = mvnrnd(mu', R)';
    
end

% Remove states after the end of time
tau(end)=[];
state(:,end) = [];

end