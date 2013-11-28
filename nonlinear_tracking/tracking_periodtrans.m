function [ period, prob ] = tracking_periodtrans( model, param, period )
% Sample and evaluate transition density for changepoint period.

% Set hyperparameters
tau_shape = model.period_trans_shape;
tau_scale = model.period_trans_scale;
tau_shift = model.period_trans_shift;

if (nargin < 3) || isempty(period)
    period = tau_shift + gamrnd(tau_shape, tau_scale);
end

if nargout > 1
    prob = log(gampdf(period-tau_shift, tau_shape, tau_scale));
end

end

