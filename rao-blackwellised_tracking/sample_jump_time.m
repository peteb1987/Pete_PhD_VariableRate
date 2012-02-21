function [ new_tau, new_m, new_u, prob ] = sample_jump_time( flags, params, old_tau, no_jump_t, new_tau, new_m, new_u )
%SAMPLE_JUMP_TIME Sample next jump time given previous jump time and last
%known time before which a jump has not occured.

if isempty(no_jump_t)
    no_jump_t = old_tau;
end

% Sample jump time transition density
if nargin == 4
    lower_lim = gamcdf(no_jump_t-old_tau, params.rate_shape, params.rate_scale);
    tmp = unifrnd(lower_lim, 1);
    new_tau = old_tau + gaminv(tmp, params.rate_shape, params.rate_scale);
    new_m = unidrnd(2);
    if new_m == 1
        new_u = mvnrnd(0, params.accel_var);
    elseif new_m == 2
        new_u = mvnrnd(0, params.tr_var);
    end
end

% Calculate probability
if nargout > 3
    
    % Jump time
    if ~isempty(new_tau)
        prob = log(gampdf(new_tau-old_tau, params.rate_shape, params.rate_scale)) ...
            -log(1-gamcdf(no_jump_t-old_tau, params.rate_shape, params.rate_scale));
    else
        prob = log(1-gamcdf(old_tau, params.rate_shape, params.rate_scale)) ...
            -log(1-gamcdf(no_jump_t-old_tau, params.rate_shape, params.rate_scale));
    end
    
    % Parameter
    if new_m == 1
        prob = prob + log(mvnpdf(new_u, 0, params.accel_var));
    elseif new_m == 2
        prob = prob + log(mvnpdf(new_u, 0, params.tr_var));
    end
    
    % Model
    prob = prob + log(0.5);
    
end

end

