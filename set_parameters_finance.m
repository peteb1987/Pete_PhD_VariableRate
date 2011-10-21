% Application
flags.rb = true;

% Model

% mn=mean, sd=standard deviation, fb=feedback

params.state_dim = 2;
params.obs_dim = 1;

params.x_start = 0;
params.x_jump_rate = 0.1;
params.x_jump_mn = 0;
params.x_jump_sd = 10;
params.x_fb = 0;
params.x_sd = 1;

params.xdot_start = 0;
params.xdot_jump_rate = 0.1;
params.xdot_jump_mn = 0;
params.xdot_jump_sd = 10;
params.xdot_fb = -0.2;
params.xdot_sd = 1;

params.H = [1 0];
params.obs_sd = 10;

params.F = [params.x_fb 1; 0 params.xdot_fb];
params.C = [params.x_sd^2 0; 0 params.xdot_sd^2];
params.R = params.obs_sd^2;

% Algorithm
params.Np = 100;
params.x_start_sd = 40;
params.xdot_start_sd = 40;