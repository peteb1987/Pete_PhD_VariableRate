function [ A, Q ] = diffusion_dynamics( F, C, dt )
%DIFFUSION_COVAR Calculate the state transition matrix and covariance for a
% brownian motion diffusion

n = size(F,1);

E = expm([F C; zeros(n) -F']*dt)*[zeros(n); eye(n)];

M1 = E(1:n, :);
M2 = E(n+1:end, :);

Q = M1/M2;

A = expm(F);

end