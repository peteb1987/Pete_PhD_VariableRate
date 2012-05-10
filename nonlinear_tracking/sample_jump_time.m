function [ new_tau, prob ] = sample_jump_time( flags, params, old_tau, no_jump_t, new_tau )
%SAMPLE_JUMP_TIME Sample next jump time given previous jump time and last
%known time before which a jump has not occured.

if isempty(no_jump_t)
    no_jump_t = old_tau;
end

% Sample jump time transition density
if nargin == 4
    lower_lim = gamcdf(no_jump_t-old_tau, params.rate_shape, params.rate_scale);
    u = unifrnd(lower_lim, 1);
    new_tau = old_tau + gaminv(u, params.rate_shape, params.rate_scale);
end

% Calculate probability
if ~isempty(new_tau)
    prob = log(gampdf(new_tau-old_tau, params.rate_shape, params.rate_scale)) ...
          -log(1-gamcdf(no_jump_t-old_tau, params.rate_shape, params.rate_scale));
else
    prob = log(1-gamcdf(no_jump_t-old_tau, params.rate_shape, params.rate_scale));
end

end

