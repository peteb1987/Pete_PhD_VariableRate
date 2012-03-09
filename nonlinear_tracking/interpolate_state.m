function [ intx, lhood ] = interpolate_state( flags, params, tau, x, w, times, observs )
%INTERPOLATE_STATE Interpolate state over an array of times, given the most
%recent state and the current accelerations. Also calculate the likelihood.

if isempty(times)
    intx = zeros(4,0);
    lhood = [];
    return
end

% Interpolate
intx = next_state(flags, params, x, w, times'-tau);

% Likelihood
obs_mean = observation_mean(flags, params, intx, w);
lhood = log_mvnpdf(observs', obs_mean', params.R);

end

