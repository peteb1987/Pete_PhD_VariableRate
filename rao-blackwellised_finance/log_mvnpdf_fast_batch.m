function y = log_mvnpdf_fast_batch(X, Mu, Sigma)
%MVNPDF Multivariate normal probability density function (pdf) for fast
%batch calculation.
%
%   Y = MVNPDF(X,MU,SIGMA) returns a vector of multivariate normal density
%   values. X is a DxN array of points. Mu is a DxN array of means. Sigma
%   is a DxDxN array of covariance matrices.
%
%   Based on mvnpdf from the stats toolbox with all the error checking and
%   special cases stripped out. Do not use this while debugging.

% Get size of data
[d,n] = size(X);

% Center data
X0 = bsxfun(@minus, X, Mu);

% Create array of standardized data, and vector of log(sqrt(det(Sigma)))
xRinv = zeros(d,n,superiorfloat(X0,Sigma));
logSqrtDetSigma = zeros(1,n,class(Sigma));
for i = 1:n
    % Make sure Sigma is a valid covariance matrix
    R = chol(Sigma(:,:,i));

    xRinv(:,i) = R' \ X0(:,i);
    logSqrtDetSigma(i) = sum(log(diag(R)));
end

% The quadratic form is the inner products of the standardized data
quadform = sum(xRinv.^2, 1);

y = (-0.5*quadform - logSqrtDetSigma - d*log(2*pi)/2)';
