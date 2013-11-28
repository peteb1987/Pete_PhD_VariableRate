function [ y_mn ] = tracking_h( model, x )
%nlng_h Deterministic observation function for the nonlinear non-Gaussian
%benchmark model.

y_mn = zeros(model.do,1);

% Bearing
y_mn(1) = atan2(x(2), x(1));

% Range
y_mn(2) = sqrt( x(1)^2 + x(2)^2 );

end
