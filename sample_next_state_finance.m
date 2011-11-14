function [ next_tau, type, next_mu, next_P ] = sample_next_state_finance( flags, params, last_t, last_tau, last_mu, last_P )
%SAMPLE_NEXT_TAU Generates the next variable jump time given the last one

% last_tau is the previous jump time. last_t is the most recent time at
% which a jump is known not to have occured.
% The remaining last_ variables describe the last state values at last_tau

% Propose a new state from the transition density (conditional on no jumps
% before last_t).

F = params.F; C = params.C;

% Poisson process jumps - exponential inter-jump times
next_x_tau = last_t + rande(1/params.x_jump_rate);
next_xdot_tau = last_t + rande(1/params.xdot_jump_rate);
next_tau = min(next_x_tau, next_xdot_tau);
type = (next_tau==next_xdot_tau)+1;

% Diffusion
[A, Q] = lti_disc(F,eye(2),C,(next_tau-last_t));
next_mu = A*last_mu;
next_P = A*last_P*A'+Q;

% Jump
if type==1
    % x jump
    next_P = next_P + [params.x_jump_sd^2 0; 0 0];
elseif type==2
    % xdot jump
    next_P = next_P + [0 0; 0 params.xdot_jump_sd^2];
end

end

