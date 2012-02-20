function [lhood, new_mu, new_P] = interp_and_lhood(flags, params, old_mu, old_P, obs, dt, type)
%INTERP_AND_LHOOD interpolate the state at a given time and calculate the
%likelihood at this point

% Calculate covariance matrix
[A, Q] = lti_disc(params.F,params.L,params.C,dt);
if type == 1
    Q = Q + [params.x_jump_sd^2 0; 0 0];
elseif type == 2
    Q = Q + [0 0; 0 params.xdot_jump_sd^2];
end

% Update Kalman filter
[p_mu, p_P] = kf_predict(old_mu, old_P, A, Q);
[new_mu, new_P, ~, ~, ~, lhood] = kf_update(p_mu, p_P, obs', params.H, params.R);
lhood = log(lhood);

end

