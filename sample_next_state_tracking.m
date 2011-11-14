function [ next_tau, w, next_x ] = sample_next_state_tracking( flags, params, last_t, last_tau, last_x, last_w )
%SAMPLE_NEXT_TAU Generates the next variable jump time given the last one

% last_tau is the previous jump time. last_t is the most recent time at
% which a jump is known not to have occured.
% The remaining last_ variables describe the last state values at last_tau

% Propose a new state from the transition density (conditional on no jumps
% before last_t).

% Gamma distributed inter-jump times
lower_lim = gamcdf(last_t-last_tau, params.rate_shape, params.rate_scale);
u = unifrnd(lower_lim, 1);
next_tau = last_tau + gaminv(u, params.rate_shape, params.rate_scale);
next_x = tracking_calc_next_state(flags, last_x, next_tau-last_tau, last_w);
w = mvnrnd(zeros(1, params.rnd_dim), params.Q)';

end

