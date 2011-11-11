function [ mu ] = tracking_calc_obs_mean( flags, params, x )
%TRACKING_CALC_OBS_MEAN Calculate the Gaussian mean of the observation
%distribution for various observation models

if flags.obs_mod == 1
    % Direct observation
    mu = x(1:params.obs_dim);
elseif flags.obs_mod == 2
    % Radar observation
    mu = zeros(params.obs_dim,1);
    [mu(1), mu(2)] = cart2pol(x(1), x(2));
    if params.obs_dim == 4
        % Bearing and range rate
        mu(4) = x(4)*cos(x(3)-mu(1));
        mu(3) = x(4)*sin(x(3)-mu(1))/mu(2);
    end
end

end

