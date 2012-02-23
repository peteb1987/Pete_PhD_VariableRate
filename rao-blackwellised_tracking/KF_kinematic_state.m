function [ mu, P, lhood ] = KF_kinematic_state( flags, params, k, cp, mu, P, times, observs  )
%KF_KINEMATIC_STATE Kalman filter the kinematic state given the nonlinear
%parameters.

Ns = cp.Ns;
tau = cp.tau(Ns);
m = cp.m(Ns);
u = cp.u(Ns);

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
H = [1 0 0 0 0 0; 0 1 0 0 0 0];
[mu, P, ~, ~, ~, lin_lhood] = kf_update(mu, P, observs(:,k), H, params.R);
lhood = log(lin_lhood);

end