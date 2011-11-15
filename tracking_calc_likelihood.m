function [ lhood, interp_x ] = tracking_calc_likelihood(flags, params, x, tau, w, times, observ)
%TRACKING_CALC_LIKELIHOOD Calculate the likelihood of all observations
%since the last state given the accelerations (random process variables)

% Work out where to start
start_idx = find(min(times(times>tau))==times);

% Initialise likelihood
lhood = zeros(length(times)-start_idx,1);

% Array for interpolated states
interp_x = zeros(1, length(times)-start_idx+1, params.state_dim);

% Loop through time
for k = start_idx:length(times)
    
    u = tracking_calc_next_state(flags, x, times(k)-tau, w);
    mu = tracking_calc_obs_mean(flags, params, u);
    lhood(k-start_idx+1) = log_mvnpdf_fast_batch(observ(:,k), mu, params.R);
    
    interp_x(1,k,:) = u;
    
end

end

