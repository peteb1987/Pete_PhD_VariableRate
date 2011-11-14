function [pts_mu, pts_P] = initialise_state_finance(flags, params, Np, MNJ, observ)
%INITIALISE_STATE Initialise particles

ds = params.state_dim; do = params.obs_dim;

pts_mu = zeros(Np, MNJ, ds);                   % Particle means (for rb case)
pts_P = zeros(Np, MNJ, ds, ds);                % Particle covariances (for rb case)

x_init = [observ(1), zeros(1,ds-do)];          % Initialise value as first observation padded with zeros (i.e. partial linear observation)
pts_mu(:,1,:) = repmat(x_init, Np, 1);
pts_P(:,1,:,:) = permute(repmat([params.x_start_sd^2, 0; 0, params.xdot_start_sd^2], [1, 1, Np]), [3,1,2]);
pts_x(:,1,:) = repmat(x_init, Np, 1);

end

