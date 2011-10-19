% Application
flags.rb = true;

% Model
params.state_dim = 2;
params.obs_dim = 1;

params.x_start = 0;
params.x_jump_rate = 0.1;
params.x_jump_mu = 0;
params.x_jump_sigma = 10;
params.x_lambda = 0;
params.x_sigma = 0.1;

params.xdot_start = 0;
params.xdot_jump_rate = 0.1;
params.xdot_jump_mu = 0;
params.xdot_jump_sigma = 10;
params.xdot_lambda = -0.2;
params.xdot_sigma = 1;

params.obs_sigma = 10;

params.F = [params.x_lambda 1; 0 params.xdot_lambda];
params.C = [params.x_sigma^2 0; 0 params.xdot_sigma^2];

% Algorithm
params.Np = 1000;
params.x_start_sigma = 40;
params.xdot_start_var = 40;