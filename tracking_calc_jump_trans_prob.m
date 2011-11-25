function [ prob ] = tracking_calc_jump_trans_prob( params, t, next_t )
%TRACKING_CALC_JUMP_TRANS_PROB Calculate probability of one jump time given
%the previous one. If no second time is given, calculate the probability of
%no subsequent jump.

if (nargin == 3)&&(~isempty(next_t))
    prob = log(gampdf(next_t-t, params.rate_shape, params.rate_scale));
else
    prob = log(1-gamcdf(t, params.rate_shape, params.rate_scale));
end


end

