function [w_new] = determine_acceleration(flags, params, old_x, new_x, dt)
%DETERMINE_ACCELERATION Calculate the acceleration vector given two bounding
%states and times

sd = flags.space_dim;

% Get old state
old_r = old_x(1:sd);
old_v = old_x(sd+1:2*sd);
old_sdot = magn(old_v);

% Get new state
new_r = new_x(1:sd);
new_v = new_x(sd+1:2*sd);
new_sdot = magn(new_v);

% Translate to standard form (start at origin)
if flags.space_dim == 2
    % Calculate 2D bearing
    old_psi = atan2(old_v(2), old_v(1));
    new_psi = atan2(new_v(2), new_v(1));
    
elseif flags.space_dim == 3
    % Find unit vectors
    et = unit(old_v,1);
    eb = unit(cross(old_v, new_v),1);
    en = cross(eb, et);
    
    R = [et, en, eb];
    
    % Translate into planar motion
    new_udot = R\new_v;
    
    old_psi = 0;
    new_psi = atan2(new_udot(2), new_udot(1));
    
end

% Differences
dpsi = new_psi - old_psi;
while abs(dpsi)>pi
    dpsi=dpsi-sign(dpsi)*2*pi;
end
dsdot = new_sdot - old_sdot;

% Find accelerations
aT = dsdot/dt;
if abs(aT)>1E-10
    aNc = aT*dpsi/log(new_sdot/old_sdot);
else
    aNc = dpsi*old_sdot/dt;
end

if flags.space_dim == 2
    aN = aNc;
    
elseif flags.space_dim == 3
    % Resolve aN into vertical and horizonal plane components
    vert_norm = unit(cross(et, [0;0;1]),1);
    hoz_norm = unit(cross(vert_norm, et),1);
    aNh = aNc*(en'*vert_norm);
    aNv = aNc*(en'*hoz_norm);
    aN = [aNv; aNh];
end

% Make an acceleration vector with no linear accelerations
w_new = [aT; aN; zeros(sd,1)];

% Work out what the new state would be with this
x_nonlinonly = next_state(flags, params, old_x, w_new, dt);

assert(all(abs(x_nonlinonly(sd+1:end)-new_x(sd+1:end))<1E-10))

% Find the difference
dx = new_x(1:sd) - x_nonlinonly(1:sd);
aX = dx/dt;

% Put it all together
w_new = [aT; aN; aX];

end

