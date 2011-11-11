function [ lhood ] = tracking_calc_likelihood(flags, params, x, tau, w, times, observ)
%TRACKING_CALC_LIKELIHOOD Calculate the likelihood of all observations
%since the last state given the accelerations (random process variables)

% Initialise likelihood
lhood = 0;

% Work out where to start
start_idx = find(min(times>tau)==times);

% Loop through time
for k = start_idx:length(times)
    
    u = tracking_calc_next_state(flags, x, times(k)-tau, w);
    mu = tracking_calc_obs_mean(flags, params, u);
    lhood = lhood + log_mvnpdf_fast_batch(observ(:,k), mu, params.R);
    
end

end

