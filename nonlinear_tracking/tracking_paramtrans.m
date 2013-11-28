function [ param, prob ] = tracking_paramtrans( model, last_param, param )
% Sample and evaluate transition density for changepoint parameters.

% Set hyperparameters depending on whether this is the first frame or not
param_mn = model.param_trans_mn;
param_vr = model.param_trans_vr;

if (nargin < 3) || isempty(param)
    param = mvnrnd(param_mn, param_vr)';
end

if nargout > 1
    prob = loggausspdf(param, param_mn, param_vr);
end

end

