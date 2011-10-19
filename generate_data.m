function [state_x, state_tau, y, times, interp_x] = generate_data( flags, params )
%GENERATE_DATA Generate data for variable rate inference problems

if flags.app == 1
    [state_x, state_tau, y, times, interp_x] = generate_data_finance( params );
elseif flags.app == 1
    [state_x, state_tau, y, times, interp_x] = generate_data_tracking( params );
end

end

