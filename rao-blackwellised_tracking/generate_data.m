function [all_x, cp_x, cp_tau, cp_m, cp_u, times, observs] = generate_data( flags, params )
%GENERATE_DATA Generate model-switching tracking data

% This version uses a regular time grid.

dt = params.dt;
R = params.R;
ds = params.state_dim; do = params.obs_dim;

% Generate changepoint times and marks
cp_tau = 0;
cp_u = [1; 0];
cp_x = params.start_state;
Ns = 1;
last_tau = 0;
while (Ns==0)||(cp_tau(Ns)<params.T)
    Ns = Ns + 1;
    cp_tau(1, Ns) = last_tau + gamrnd(params.rate_shape, params.rate_scale);
    cp_u(:, Ns) = [unidrnd(2); 0];
    if cp_u(1, Ns) == 1
        cp_u(2, Ns) = mvnrnd(0, params.accel_var);
    elseif cp_u(1, Ns) == 2
        cp_u(2, Ns) = (unidrnd(2)*2-3) * raylrnd(sqrt(2*params.tr_var/(4-pi)));
    end
    last_tau = cp_tau(Ns);
end

% Create arrays
times = cumsum(dt*ones(params.K,1));            % Array of times
all_x = zeros(ds,Ns);                           % Array of states
observs = zeros(do,params.K);                   % Array of observsations

% Initialise for loop
ti=1;                               % jump counter
last_x = params.start_state;             % starting state
last_t = 0;

% Loop through observsation times and generate states
for k=1:params.K
    
    t = times(k);
    
    interm_t = last_t;
    interm_x = last_x;
    
    % Iteratively sample forwards to the next jump and through it
    while cp_tau(ti+1) < t
        
        [A, Q, ~] = construct_transmats(cp_tau(ti+1)-interm_t, cp_u(1,ti), cp_u(2,ti), params.proc_var);
        
        % Increment jump counter
        ti = ti + 1;
        interm_t = cp_tau(ti);
        
        % Sample state
        x = mvnrnd((A*interm_x)', Q)';
        if magn(x(3:4))>params.max_vel
            x(3:4) = params.max_vel*x(3:4)/magn(x(3:4));
        end
        
        % Set new acceleration
        [~, ~, Ajump] = construct_transmats(0, cp_u(1,ti), cp_u(2,ti), params.proc_var);
        x = Ajump*x;
        
        % Store
        cp_x(:, ti) = x;
        
        interm_x = x;
        
    end
    
    [A, Q, ~] = construct_transmats(t-interm_t, cp_u(1,ti), cp_u(2,ti), params.proc_var);
    
    x = mvnrnd((A*interm_x)', Q)';
    if magn(x(3:4))>params.max_vel
        x(3:4) = params.max_vel*x(3:4)/magn(x(3:4));
    end

    % Store state
    all_x(:,k) = x;
	
    % Sample observsation
    mu = observation_mean(flags, params, x);
    observs(:,k) = mvnrnd(mu', R)';
    
    last_t = t;
    last_x = x;
    
end

% Remove states after the end of time
cp_tau = cp_tau(1,1:ti);
cp_m = cp_u(1,1:ti);
cp_u = cp_u(2,1:ti);

end