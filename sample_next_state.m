function [ next_x, next_tau ] = sample_next_state( flags, params, last_x, last_tau )
%SAMPLE_NEXT_TAU Generates the next variable jump time given the last one

if flags.app == 1
    % Poisson process jumps - exponential inter-jump times
    nextx_x_tau = last_tau + rande(1/rate);
elseif flags.app == 1
    % Gamma distributed inter-jump times
    
    % NOT WRITTEN YET
    
end

end

