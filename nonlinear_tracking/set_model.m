% Set model parameters

% Basics
model.K = 1800;
model.fs = 30;
model.dw = 30;
model.dp = 1;
model.np = 2;
model.num_sens = 4;

% Observation model
model.y_obs_vr = 0.1^2;                         % Observation variance

% Template model
mw = zeros(model.dw,1);% model.w_prior_mn(10) = 1; model.w_prior_mn(15) = -1;
Pw = 0.2^2 * diag(0.5*(1-cos(2*pi*(1:model.dw)/(model.dw+1))).^2); %0.2^2*eye(model.dw);

if model.np == 1
    model.w_prior_vr = blkdiag(Pw, Pw, Pw, Pw);
    model.w_prior_mn = [mw; mw; mw; mw];
elseif model.np == 2
    model.w_prior_vr = blkdiag(Pw, Pw, Pw, Pw, Pw, Pw, Pw, Pw);
    model.w_prior_mn = [mw; mw; mw; mw; mw; mw; mw; mw];
end

% Beat offset model
model.p_prior_shape = 40;
model.p_prior_scale = 0.02;
model.p_trans_scale = 1E-3;

% Beat period model
model.tau_trans_shape = 1.5;
model.tau_trans_scale = 0.1/1.5;
