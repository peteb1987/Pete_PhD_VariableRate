% Application
flags.rb = true;

% Model
params.x_start = 0;
params.x_jump_rate = 0.1;
params.x_jump_mu = 0;
params.x_jump_sigma = 10;
params.x_decay = 1;
params.x_sigma = 0.1;

params.xdot_start = 0;
params.xdot_jump_rate = 0.1;
params.xdot_jump_mu = 0;
params.xdot_jump_sigma = 10;
params.xdot_decay = 0.8;
params.xdot_sigma = 1;

params.obs_sigma = 10;

% Algorithm
params.Np = 1000;