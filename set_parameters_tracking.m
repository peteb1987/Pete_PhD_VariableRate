% Application
flags.rb = true;

% Data generation
params.K = 500;
params.dt = 0.1;
params.T = params.dt*params.K;

% Model
params.state_dim = 4;                   % State dimension
params.obs_dim = 2;                     % Observation dimension

params.rate_shape = 5;                  % State time gamma distribution shape parameter (this is "a", or "k")
params.rate_scale = 1;                  % State time gamma distribution scale parameter (this is "b", or "theta")

if flags.dyn_mod == 1
    params.Q = diag([0.01, sqrt(10)]);            % Covariance for aT, aP
    params.rnd_dim = 2;                     % Random variables dimension
elseif flags.dyn_mod == 2
    params.Q = diag([0.01, sqrt(10), 1, 1]);      % Covariance for aT, aP, aX, aY
    params.rnd_dim = 4;                     % Random variables dimension
elseif flags.dyn_mod == 3
    params.Q = diag([0.01, sqrt(10), (pi/90)^2, 0.01]);      % Covariance for aT, aP, aB, aS
    params.rnd_dim = 4;                     % Random variables dimension
end
if flags.obs_mod == 1
    if params.obs_dim == 2
        params.C = [1 0 0 0; 0 1 0 0];
        params.R = 10*eye(2);
    else
        params.R = diag([10, 10, 1, 1]);
    end
elseif flags.obs_mod == 2
    if params.obs_dim == 2
        params.R = diag([(pi/90)^2, 1]);
    else
        params.R = diag([(pi/90)^2, 10, (pi/90)^2, 0.1]);
    end
end

params.start_state = [50; 50; 0; 5];
params.min_speed = 0.5;


% Algorithm
params.Np = 500;            % Target number of filtering particles
params.S = 100;             % Number of smoothing trajectories
params.start_var = diag([10, 10, (pi/90)^2, 0.1]);
params.start_bng = 0;
params.start_speed = 5;