function [ H ] = tracking_obsjacobian( model, x )
%tracking_obsjacobian Calculate observation function jacobian for the
%tracking model.

rng_sq = x(1)^2 + x(2)^2;
rng = sqrt(rng_sq);

dh1_dx = [-x(2), x(1), 0, 0]/rng_sq;
dh2_dx = [x(1), x(2), 0, 0]/rng;

H = [dh1_dx; dh2_dx];
     
end

