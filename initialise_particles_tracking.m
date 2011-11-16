function [pts] = initialise_particles_tracking(flags, params, observ)
%INITIALISE_STATE Initialise particles

ds = params.state_dim;
do = params.obs_dim;
dr = params.rnd_dim;
Np = params.Np;

% Find a mean first state from the first observation
if flags.obs_mod == 1
    x_init = [observ(1:2,1)', params.start_bng, params.start_speed];
elseif flags.obs_mod == 2
    [x1, x2] = pol2cart(observ(1,1), observ(2,1));
    x_init = [x1, x2, params.start_bng, params.start_speed];
end

% Generate a random set of starting points and accelerations
start_pos = mvnrnd(repmat(x_init, Np, 1), params.start_var)';
start_pos(4,:) = max(start_pos(4,:), params.min_speed);
start_accel = mvnrnd(zeros(Np,dr), params.Q)';

% Put the values into cell array
pts = struct('x', num2cell(start_pos, [4,1])',...
             'w', num2cell(start_accel, [4,1])',...
             'tau', 0,...
             'Ns', 1,...
             'intx', num2cell(start_pos, [4,1])',...
             'lhood', 0);

end

