function [lhood, intx, intmu, intP] = interp_and_lhood(flags, params, tau, t, w, x, mu, P, obs)
%INTERP_AND_LHOOD interpolate the state at a given time and calculate the
%likelihood at this point

intx = zeros(size(x));
intmu = zeros(size(mu));
intP = zeros(size(P));

if flags.app == 1
    % Update Kalman filter to current time
    [A, Q] = lti_disc(params.F,eye(2),params.C,t-tau);
    [p_mu, p_P] = kf_predict(mu, P, A, Q);
    [intmu, intP, ~, ~, ~, lhood] = kf_update(p_mu, p_P, obs', params.H, params.R);
    lhood = log(lhood);
elseif flags.app == 2
    % Interpolate state
    intx = tracking_calc_next_state(flags, x, t-tau, w);
    obs_mn = tracking_calc_obs_mean(flags, params, intx);
    % Unwrap bearing rate
    if (params.obs_dim==4)&&(abs(obs(3)-obs_mn(3))>pi)
        obs_mn(3) = obs_mn(3) - 2*pi*round((obs_mn(3)-obs(3))/(2*pi));
    end
    % Calculate likilood
    lhood = log(mvnpdf(obs', obs_mn', params.R));
end

end

