function [ J ] = tracking_statejacobian( model, t0, u, x0, t, x )
%TRACKING_STATEJACOBIAN Calculate Jacobians for state with respect
%topreceeding changepoint parameters.

% Changepoint is at (t0, x0) with parameters u (with which we want to
% differentiate. Current point is (t,x).

% Unpack
aT = u(1);
aN = u(2);
r0 = x0(1:2);
s0 = sqrt(x0(3)^2+x0(4)^2);
psi0 = atan2(x0(4),x0(3));
r = x(1:2);
s = sqrt(x(3)^2+x(4)^2);
psi = atan2(x(4),x(3));

% Intermediate quantities
M0 = [cos(psi0), sin(psi0); sin(psi0), -cos(psi0)];
M = [cos(psi), sin(psi); sin(psi), -cos(psi)];
alpha = [2*aT; aN]/(4*aT^2+aN^2);

% Differentiate them
dM_dpsi = [-sin(psi), cos(psi); cos(psi), sin(psi)];
dalpha_da = [2*(aN^2-4*aT^2), -4*aT*aN   ;
             -8*aT*aN     , 4*aT^2-aN^2]/(( 4*aT^2+aN^2 )^2);

% Differentiate speed
ds_daT = t - t0;
ds_daN = 0;
ds_da = [ds_daT, ds_daN];

% Differentiate bearing
dpsi_daT = -(psi-psi0)/aT + (aN/aT)*ds_daT/s;
dpsi_daN =  (psi-psi0)/aN;
dpsi_da = [dpsi_daT, dpsi_daN];

% Differentiate position
dr_da = ( s^2*M - s0^2*M0 )*dalpha_da ...
       + (s^2)*dM_dpsi*alpha*dpsi_da ...
       + 2*s*M*alpha*ds_da;

% Concatenate it all
dx_da = [dr_da; ds_da; dpsi_da];

% Correct it to cartesian
J = [1, 0, 0,         0         ;
     0, 1, 0,         0         ;
     0, 0, cos(psi), -s*sin(psi);
     0, 0, sin(psi),  s*cos(psi)] * dx_da;

end

