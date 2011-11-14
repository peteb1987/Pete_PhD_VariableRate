function [lhood, intmu, intP] = interp_and_lhood_finance(flags, params, tau, t, mu, P, obs)
%INTERP_AND_LHOOD interpolate the state at a given time and calculate the
%likelihood at this point

% Update Kalman filter to current time
[A, Q] = lti_disc(params.F,eye(2),params.C,t-tau);
[p_mu, p_P] = kf_predict(mu, P, A, Q);
[intmu, intP, ~, ~, ~, lhood] = kf_update(p_mu, p_P, obs', params.H, params.R);
lhood = log(lhood);

end

