function [all_x, all_tau, observs, times, all_intx, all_w] = generate_data( flags, params )
%GENERATE_DATA_FINANCE Generate curvilinear tracking data

% This version uses a regular time grid.

dt = params.dt;
Q = params.Q;
R = params.R;
ds = params.state_dim; do = params.obs_dim; dr = params.rnd_dim;

% Generate state times
all_tau = 0;
Ns = 1;
last_tau = 0;
while (Ns==0)||(all_tau(Ns)<params.T)
    Ns = Ns + 1;
    all_tau(Ns,1) = last_tau + gamrnd(params.rate_shape, params.rate_scale);
    last_tau = all_tau(Ns,1);
end

% Create arrays
times = cumsum(dt*ones(params.K,1))-dt;         % Array of times
all_x = zeros(ds,Ns);                           % Array of states
all_w = zeros(dr,Ns);                           % Array of random variables
all_intx = zeros(ds,params.K);                  % Array of interpolated states
observs = zeros(do,params.K);                   % Array of observsations

% Initialise for loop
ti=2;                               % jump counter
w = mvnrnd(zeros(1,dr), Q)';        % starting r.v.s
x = params.start_state;             % starting state
last_x = x;
last_tau = 0;

% Store first jump values
all_x(:, 1) = x;
all_w(:, 1) = w;
mu = observation_mean(flags, params, x, w);
observs(:,1) = mvnrnd(mu', R)';
all_intx(:,1) = x;

% Loop through observsation times and generate states
for k=2:params.K
    
    t = times(k);
    
    % Iteratively sample forwards to the next jump and through it
    while all_tau(ti) < t
        
        % New state
        x = next_state(flags, params, last_x, w, all_tau(ti)-last_tau);
        w = mvnrnd(zeros(1,dr), Q)';
        
        % Store state
        all_x(:,ti) = x;
        all_w(:,ti) = w;
        
        % Update last_ variables
        last_tau = all_tau(ti);
        last_x = x;
        
        % Increment jump counter
        ti = ti + 1;
        
    end
    
    % Interpolate state
    x = next_state(flags, params, last_x, w, t-last_tau);
    
    % Store it
    all_intx(:,k) = x;
    
    % Sample observsation
    mu = observation_mean(flags, params, x, w);
    observs(:,k) = mvnrnd(mu', R)';
    
end

% Remove states after the end of time
all_tau(end)=[];
all_x(:,end) = [];

end