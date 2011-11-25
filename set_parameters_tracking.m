%%% TRACKING PARAMETERS %%%

%%% Data generation %%%
params.K = 500;
params.dt = 0.1;
params.T = params.dt*params.K;

%%% Model %%%

% Set model dimensionalities
if flags.dyn_mod >= 5
    params.state_dim = 6;               % State dimension
    params.obs_dim = 3;                 % Observation dimension
else
    params.state_dim = 4;               % State dimension
    params.obs_dim = 2;                 % Observation dimension
end
if flags.obs_vel
    params.obs_dim = params.obs_dim * 2;
end

% Dynamic Model
if flags.dyn_mod == 1
    params.Q = diag([0.01, 10]);                        % Covariance for aT, aP
    params.rnd_dim = 2;                                 % Random variables dimension
elseif flags.dyn_mod == 2
    params.Q = diag([0.01, 10, 1, 1]);                  % Covariance for aT, aP, aX, aY
    params.rnd_dim = 4;                                 % Random variables dimension
elseif flags.dyn_mod == 3
    params.Q = diag([0.01, 10, (pi/90)^2, 0.01]);       % Covariance for aT, aP, aB, aS
    params.rnd_dim = 4;                                 % Random variables dimension
elseif flags.dyn_mod == 4
    params.Q = diag([0.01, 10, (pi/90)^2, 0.1]);        % Covariance for aT, aP, aB, aS
    params.rnd_dim = 4;                                 % Random variables dimension
elseif flags.dyn_mod == 5
    params.Q = diag([0.01, 1, 1]);                      % Covariance for aT, aPv, aPh
    params.rnd_dim = 3;                                 % Random variables dimension
elseif flags.dyn_mod == 6
    params.Q = diag([0.01, 1, 1, 1, 1, 1]);             % Covariance for aT, aPv, aPh, aX, aY, aZ
    params.rnd_dim = 6;                                 % Random variables dimension
end

% Jump time model
params.rate_shape = 5;                  % State time gamma distribution shape parameter (this is "a", or "k")
params.rate_scale = 1;                  % State time gamma distribution scale parameter (this is "b", or "theta")

% Observation model
if flags.obs_mod == 1
    if flags.obs_vel
%         params.C = [1 0 0 0; 0 1 0 0];
        params.C = [eye(params.obs_dim), zeros(params.obs_dim, params.state_dim-params.obs_dim)];
        params.R = 10*eye(params.obs_dim);
    else
        params.C = eye(params.state_dim);
        params.R = diag([10, 10, 1, 1]);
    end
elseif flags.obs_mod == 2
    if ~flags.obs_vel
        if params.obs_dim == 2
            params.R = diag([(pi/90)^2, 10]);
        elseif params.obs_dim == 3
            params.R = diag([(pi/90)^2, (pi/180)^2, 10]);
        end
    else
        if params.obs_dim == 4
            params.R = diag([(pi/90)^2, 10, (pi/90)^2, 0.1]);
        elseif params.obs_dim == 6
            params.R = diag([(pi/90)^2, (pi/180)^2, 10, (pi/30)^2, (pi/30)^2, 0.1]);
        end
    end
end

% Other parameters
params.min_speed = 0.5;


% Algorithm
params.Np = 50;            % Target number of filtering particles
params.S = 50;             % Number of smoothing trajectories
params.start_bng = 0;
params.start_speed = 5;
params.ppsl_move_time_sd = 0.1*(params.rate_shape*params.rate_scale);
params.opt_ppsl_window_length = 5;

if flags.dyn_mod >= 5
    params.start_state = [50; 50; 50; 0; 5; 0];
    params.start_var = diag([10, 10, 10, 0.1, 0.1, 0.1]);
else
	params.start_state = [50; 50; 0; 5];
    params.start_var = diag([10, 10, (pi/90)^2, 0.1]);
end