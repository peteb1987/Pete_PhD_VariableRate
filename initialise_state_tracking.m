function [pts_x, pts_w] = initialise_state_tracking(flags, params, Np, MNJ, observ)
%INITIALISE_STATE Initialise particles

ds = params.state_dim; do = params.obs_dim; dr = params.rnd_dim;

pts_x = zeros(Np, MNJ, ds);                    % Particle states
pts_w = zeros(Np, MNJ, dr);

if flags.obs_mod == 1
    x_init = [observ(1:2,1)', params.start_bng, params.start_speed];
elseif flags.obs_mod == 2
    [x1, x2] = pol2cart(observ(1,1), observ(2,1));
    x_init = [x1, x2, params.start_bng, params.start_speed];
end
start_pos = mvnrnd(repmat(x_init, Np, 1), params.start_var);
start_pos(:,4) = max(start_pos(:,4), params.min_speed);
pts_x(:,1,:) = start_pos;
pts_w(:,1,:) = mvnrnd(zeros(Np,dr), params.Q);

end

