function ld = loggausspdf(x, m, S)
%LOGGAUSSPDF Calculates the log of a Gaussian pdf
%
%   LD = LOGGAUSSPDF(X,MU,SIGMA) returns the log-Gaussian density of a
%   point given its mean and covariance. X and MU are a Dx1 vectors. SIGMA 
%   is a DxD covariance matrix.
%

% Constant
log2pi = 1.83787706640935;

% Get size of data
[d,n] = size(x);

x = x - m;
ld = -0.5 * ( sum(x.*(S\x)) + d*log2pi +  2*sum(log(diag(chol(S)))) );  %sum(log(eig(S))) %log(det(S))

end

