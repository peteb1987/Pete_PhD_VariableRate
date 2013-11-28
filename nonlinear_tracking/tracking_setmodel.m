function model = tracking_setmodel(test)

% Model parameters

% Using a 2D intrinsic coordinate model and a
% bearing-range observation model.

%%%%%%%%%%%%%%%%

% General things
model.K = 100;                  % Number of time points
model.ds = 4;                   % Dimension of the state
model.do = 2;                   % Dimension of the observations
model.dp = 2;

% Time
model.T = 1;                    % Sampling period

% Acceleration variances
model.aT_vr = 0.01;
model.aN_vr = 1;

% Changepoint period model
model.period_trans_shape = 2;
model.period_trans_scale = 9/2;
model.period_trans_shift = 1;

% Changepoint parameter model
model.param_trans_mn = zeros(model.dp,1);
model.param_trans_vr = diag([model.aT_vr, model.aN_vr]);

% Initial state
model.x0 = [-100; 50; 10; 0];

% Minimum speed
model.vmin = 0.1;

% Observation variances
sigtheta = ( 2*(pi/180) )^2; % Bearing covariance (0.25/5)
sigr     = 10;                 % Range covariance (0.1/100)

% Observation covariance matrix
model.R = diag([sigtheta sigr]);

end
