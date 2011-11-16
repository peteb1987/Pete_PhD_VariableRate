function [ lhood, interp_x ] = tracking_calc_likelihood(flags, params, x, tau, w, times, observ, start_idx)
%TRACKING_CALC_LIKELIHOOD Calculate the likelihood of all observations
%since the last state given the accelerations (random process variables)

K = length(times);

if nargin == 7
    % Work out where to start
    start_idx = find(min(times(times>tau))==times);
end

% Initialise likelihood
lhood = zeros(K,1);

% Array for interpolated states
interp_x = zeros(params.state_dim, K);

% Interpolate
u = tracking_calc_next_state_batch_time(flags, x, times(start_idx:K)-tau, w);
interp_x(:,start_idx:K) = u;

% Calculate likelihood
mu = tracking_calc_obs_mean(flags, params, u);
lhood(start_idx:K) = log(mvnpdf(observ(:,start_idx:K)', mu', params.R));

end

