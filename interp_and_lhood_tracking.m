function [lhood, intx] = interp_and_lhood_tracking(flags, params, tau, t, w, x, obs)
%INTERP_AND_LHOOD interpolate the state at a given time and calculate the
%likelihood at this point

% Interpolate state
intx = tracking_calc_next_state(flags, x, t-tau, w);
obs_mn = tracking_calc_obs_mean(flags, params, intx, w);
% Unwrap bearing rate
if (params.obs_dim==4)&&(abs(obs(3)-obs_mn(3))>pi)
    obs_mn(3) = obs_mn(3) - 2*pi*round((obs_mn(3)-obs(3))/(2*pi));
end
% Calculate likilood
lhood = log(mvnpdf(obs', obs_mn', params.R));

end

