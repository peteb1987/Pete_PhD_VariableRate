function [ x ] = tracking_calc_next_state_batch_time_3D( flags, last_x, dt, w )
%TRACKING_CALC_NEXT_STATE_BATCH_TIME_3D Batch version - handles a vector
%for dt, with a different time difference in each element

K = length(dt);
ds = length(last_x);

dt(dt==0)=eps;

% All hell breaks lose if s hit 0, so lets limit it!
min_speed = 0.5;

% Get accelerations
aT = w(1);
aNv = w(2);
aNh = w(3);
aN = sqrt(aNv.^2+aNh.^2);
if flags.dyn_mod==6
    aX = w(4:6);
else
    aX = zeros(3,1);
end

% Get old state
last_r = last_x(1:3);
last_v = last_x(4:6);
sdot = sqrt(sum(last_v.^2));
psi = 0;

% Find angle of normal acceleration relative to vertical plane
phi = atan(aNv/aNh);

% Find unit vectors and 3D rotation matrix from plane to 3D space
et = last_v/sdot;
eb = zeros(3,1); u = sqrt(et(1)^2+et(2)^2);
eb(1) =  (et(2)*sin(phi)-et(1)*et(3)*cos(phi))/u;
eb(2) = (-et(1)*sin(phi)-et(2)*et(3)*cos(phi))/u;
eb(3) = cos(phi)*u;
en = cross(eb,et);
R = [et, en, eb];

%%% In-plane modelling - start with zero bearing and sdot speed in the _i_ direction

% Calculate new speed, and handle lower limiting
new_sdot = sdot + aT*dt;
if any(new_sdot<min_speed)
    new_sdot(new_sdot<min_speed) = min_speed;
    aT = (min_speed-sdot)/dt(end);
end

% Calculate new bearing
if aT~=0
    new_psi = (aN./aT).*log(new_sdot./sdot);
else
    new_psi = (aN.*dt)./sdot;
end

% Calculate new u (within plane position)
u = zeros(3,K);
SF1 = 4*aT.^2 + aN.^2;
if (aT~=0)&&(aN~=0)
    u(1,:) = ((new_sdot.^2)./SF1).*( aN.*sin(new_psi)+2*aT.*cos(new_psi)) - ((sdot^2)./SF1)*( aN.*sin(psi)+2*aT.*cos(psi));
    u(2,:) = ((new_sdot.^2)./SF1).*(-aN.*cos(new_psi)+2*aT.*sin(new_psi)) - ((sdot^2)./SF1)*(-aN.*cos(psi)+2*aT.*sin(psi));
elseif (aT==0)&&(aN~=0)
    u(1,:) = ((new_sdot.^2)./aN).*( sin(new_psi) - sin(psi) );
    u(2,:) = ((new_sdot.^2)./aN).*(-cos(new_psi) + cos(psi) );
elseif (aT~=0)&&(aN==0)
    u(1,:) = 0.5*dt.*cos(psi).*new_sdot;
    u(2,:) = 0.5*dt.*sin(psi).*new_sdot;
else
    u(1,:) = ( sdot*dt.*cos(psi) );
    u(2,:) = ( sdot*dt.*sin(psi) );
end

% Calculate cartesian in-plane velocity
udot = zeros(3,K);
[udot(1,:), udot(2,:)] = pol2cart(new_psi, new_sdot);

% Translate position and velocity into 3D space
r = bsxfun(@plus, R*u+aX*dt', last_r);
v = R*udot;

% Stack them up
x = [r; v];

end

