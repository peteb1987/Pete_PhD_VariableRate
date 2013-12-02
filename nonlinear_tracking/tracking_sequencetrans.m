function [ cp_time, cp_param, cp_state ] = tracking_sequencetrans( model, start_time, stop_time, pre_cp_time, pre_cp_param, pre_cp_state )
%TRACKING_SEQUENCETRANS Sample a changepoint sequence from the prior.

% Make sure the preceding changepoint time is before the window.
assert(pre_cp_time<=start_time);

% Sample the first time using the inverse cdf method
min_time = pre_cp_time+model.period_trans_shift;
cp_time = min_time + inverse_cdf_gamma_rnd(model.period_trans_shape, model.period_trans_scale, start_time-min_time);
cp_param = tracking_paramtrans(model, []);
cp_state = tracking_evaluatestate(model, pre_cp_time, pre_cp_param, pre_cp_state, cp_time);

% Loop
num_new_cps = 0;
while cp_time(end) < stop_time
    
    % Sample next changepoint
    next_cp_time = cp_time(end) + tracking_periodtrans(model, cp_param(:,end));
    next_cp_param = tracking_paramtrans(model, cp_param(:,end));
    next_cp_state = tracking_evaluatestate(model, cp_time(end), cp_param(:,end), cp_state(:,end), next_cp_time);
    
    % Append it
    cp_time = [cp_time, next_cp_time];
    cp_param = [cp_param, next_cp_param];
    cp_state = [cp_state, next_cp_state];
    
    num_new_cps = num_new_cps + 1;
    
    %%% FUDGE TO STOP >1 CP in a window %%%
    if num_new_cps == 2
        break;
    end
    
end

% Remove last superfluous changepoint
cp_time(end) = [];
cp_param(:,end) = [];
cp_state(:,end) = [];

end

function tau = inverse_cdf_gamma_rnd(tau_shape, tau_scale, tau_min)

if tau_min > 0
    u_min = gamcdf(tau_min, tau_shape, tau_scale);
    u = unifrnd(u_min, 1);
    tau = gaminv(u, tau_shape, tau_scale);
else
    tau = gamrnd(tau_shape, tau_scale);
end

end
