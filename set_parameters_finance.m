% Application
flags.rb = true;

% Model

% mn=mean, sd=standard deviation, fb=feedback

params.state_dim = 2;
params.obs_dim = 1;

params.x_start = 0;
params.x_jump_rate = 10;
params.x_jump_mn = 0;
params.x_jump_sd = 0.005;
params.x_fb = 0;
params.x_sd = 0.005;

params.xdot_start = 0;
params.xdot_jump_rate = 10;
params.xdot_jump_mn = 0;
params.xdot_jump_sd = 0.0001;
params.xdot_fb = -0.05;
params.xdot_sd = 0.05;

params.H = [1 0];
params.obs_sd = 0.001;

params.F = [params.x_fb 1; 0 params.xdot_fb];
params.C = [params.x_sd^2 0; 0 params.xdot_sd^2];
params.R = params.obs_sd^2;

% Algorithm
params.Np = 100;            % Target number of filtering particles
params.S = 100;             % Number of smoothing trajectories
params.x_start_sd = 40;
params.xdot_start_sd = 40;