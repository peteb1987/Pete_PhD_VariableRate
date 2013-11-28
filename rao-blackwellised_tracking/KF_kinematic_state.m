function [ mu, P, lhood, A, Q ] = KF_kinematic_state( flags, params, k, cp, mu, P, times, observs  )
%KF_KINEMATIC_STATE Kalman filter the kinematic state given the nonlinear
%parameters.

Ns = cp.Ns;

% Find latest jump
ji = find(max(cp.tau(cp.tau<times(k)))==cp.tau(cp.tau<times(k)));

tau = cp.tau(ji);
m = cp.m(ji);
u = cp.u(ji);

t = times(k);
if k > 1
    last_t = times(k-1);
else
    last_t = 0;
end

% Transition matrices
if tau <= last_t
    [A, Q, ~] = construct_transmats(t-last_t, m, u, params.proc_var);
else
    [A, Q, ~] = construct_transmats(tau-last_t, cp.m(Ns-1), cp.u(Ns-1), params.proc_var);
    [Ainc, Qinc, Ajump] = construct_transmats(t-tau, m, u, params.proc_var);
    A = Ainc*Ajump*A;
    Q = Ainc*Ajump*Q*Ajump'*Ainc' + Qinc;
end

% Kalman Filter
[mu, P] = kf_predict(mu, P, A, Q);
[mu, P, ~, ~, ~, lin_lhood] = kf_update(mu, P, observs(:,k), params.H, params.R);
lhood = log(lin_lhood);

end