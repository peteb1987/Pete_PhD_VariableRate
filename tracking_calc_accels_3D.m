function w_new = tracking_calc_accels_3D(flags, tau, next_tau, x, next_x)
%TRACKING_CALC_ACCELS_3D Calculate the acceleration vector given two bounding
%states and times, for the 3D model.

% Get old state
r = x(1:3);
v = x(4:6);
sdot = sqrt(sum(v.^2));

% Get new state
new_r = next_x(1:3);
new_v = next_x(4:6);
new_sdot = sqrt(sum(new_v.^2));

% Find unit vectors
et = v/sdot;
eb = cross(v, new_v);
eb = eb/sqrt(sum(eb.^2));
en = cross(eb, et);

R = [et, en, eb];

% Translate into planar motion
new_u = R\(new_r-r);
udot = R\v;
new_udot = R\new_v;

assert(abs(udot(2))<1E-10);
assert(abs(udot(3))<1E-10);
assert(abs(new_udot(3))<1E-10);

% Time difference
dt = next_tau - tau;

% Find change in dpsi
[dpsi, ~] = cart2pol(new_udot(1), new_udot(2));

% Find accelerations
aT = (new_sdot - sdot)/dt;
aN = aT*dpsi/log(new_sdot/sdot);

% Resolve aN into vertical and horizonal plane components
vert_norm = cross([0;0;1], et);
vert_norm = vert_norm/sqrt(sum(vert_norm.^2));
hoz_norm = cross(et, vert_norm);
hoz_norm = hoz_norm/sqrt(sum(hoz_norm.^2));
aNh = aN*(en'*vert_norm);
aNv = aN*(en'*hoz_norm);


w_new = [aT; aNv; aNh; 0; 0; 0];

x_nonlinonly = tracking_calc_next_state_3D(flags, x, next_tau-tau, w_new);

% assert(x_should(4:6)==next_x(4:6), 'Degenerate model, fool.');

dX = next_x(1:3) - x_nonlinonly(1:3);
aX = dX/dt;

w_new = [aT; aNv; aNh; aX];

end

