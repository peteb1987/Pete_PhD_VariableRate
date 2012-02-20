function [ next_tau, next_type, prob ] = sample_next_state( flags, params, last_t, last_tau, next_tau, next_type )
%SAMPLE_NEXT_STATE Sample jump time transition density and/or return
%probability

% last_tau is the previous jump time.
% last_t is the most recent time at which a jump is known not to have occured.

% Sample next_tau if it is not given
if nargin<5
    next_x_tau = last_t + rande(1/params.x_jump_rate);
    next_xdot_tau = last_t + rande(1/params.xdot_jump_rate);
    next_tau = min(next_x_tau, next_xdot_tau);
    next_type = (next_tau==next_xdot_tau)+1;
end

% Calculate probability
if nargout > 2
    if next_type == 0
        prob = log(exppdf(next_tau-last_t, 1/(params.x_jump_rate+params.xdot_jump_rate)));
    elseif next_type == 1
        prob = log( exppdf(next_tau-last_t, 1/(params.x_jump_rate)) ) ...
            +log( params.x_jump_rate/(params.x_jump_rate+params.xdot_jump_rate) );
    elseif next_type == 2
        prob = log(exppdf(next_tau-last_t, 1/(params.xdot_jump_rate))) ...
            +log( params.xdot_jump_rate/(params.x_jump_rate+params.xdot_jump_rate) );
    end
end

end

