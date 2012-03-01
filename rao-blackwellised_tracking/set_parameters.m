%%% TRACKING PARAMETERS %%%

%%% Miscellaneous %%%
params.K = 500;                         % Number of time steps
params.dt = 0.1;                        % Sampling time (this should only be used for generating data. Otherwise, use sampling times provided)
params.T = params.dt*params.K;          % Time of last observation
params.min_speed = 0.5;                 % Minimum speed allowed

params.state_dim = 6;

% Starting point distribution (assumed known by algorithm)
% start_state is used by the data generation function (i.e. it's part of the model)
% start_var is used to initiialise particles (i.e. it's part of the algorithm)
params.start_state = [-50; 50; 5; 0; 0; 0];
params.start_var = diag([10, 10, 1, 1, 0.1, 0.1]);


%%% Model %%%

% Jump times
params.rate_shape = 5;                  % State time gamma distribution shape parameter (this is "a", or "k")
params.rate_scale = 1;                  % State time gamma distribution scale parameter (this is "b", or "theta")

% Process variance
params.proc_var = 0.01;

% Turn rate
params.tr_var = 0.01;

% Acceleration jump
params.accel_var = 0.001;

% Maximum velocity
params.max_vel = 10;

% Observations
range_var = 10;
bear_var = (pi/90)^2;
range_rate_var = 0.1;
bear_rate_var = (pi/30)^2;
x_var = 1;
x_rate_var = 1;
params.H = [1 0 0 0 0 0; 0 1 0 0 0 0];

%%% Algorithm 
params.Np = 50;                            % Target number of filtering particles
params.Ns = 10;                            % Number of smoothing trajectories

% params.M = 1;
% params.ppsl_move_time_sd = ...          % Standard deviation for proposal distribution for moving jump times
%     0.1*(params.rate_shape*params.rate_scale);
% params.min_num_ppsl_frames = 20;         % Minimum number of frames over which the UKF-approximated OID proposal is constructed
% params.prop_ppsl_frames = 0.1;          % Proportion of frames in a window used for acceleration proposal



%%% Set up secondary parameters - don't edit this bit

% Set state dimensionality
params.state_dim = 6;

% Set observation dimensionality
if flags.obs_vel
    params.obs_dim = 4;
else
    params.obs_dim = 2;
end

% Set covariance matrix for observations
if flags.obs_mod == 1
    if ~flags.obs_vel
        cov = [x_var, x_var];
    else
        cov = [x_var, x_var, x_rate_var, x_rate_var];
    end
elseif flags.obs_mod== 2
    if ~flags.obs_vel
        cov = [bear_var, range_var];
    else
        cov = [bear_var, range_var, bear_rate_var, range_rate_var];
    end
else
    error('unhandled option');
end
params.R = diag(cov);