function [ old_x ] = previous_state( flags, params, new_x, w, dt )
%NEXT_STATE Deterministically calculate the previous state given the next
%state, the random variable, and the time difference using the various
%dynamic models

% If no time has passed (i.e. we're looking at an observation exactly
% after a jump) then return straight away.
if (dt==0)
    old_x=new_x;
    return
end

% Set minimum speed
min_speed = params.min_speed;

% Get accelerations
aT = w(1,:);
aN = w(2:flags.space_dim,:);
if flags.dyn_mod == 2
    aX = w(flags.space_dim+1:end,:);
else
    aX = zeros(flags.space_dim,Ns);
end

% Get new state
new_r = new_x(1:flags.space_dim,:) - aX*dt;
new_v = new_x(flags.space_dim+1:end,:);

% Transform to planar intrinisics
new_sdot = norm(new_v);

if flags.space_dim == 2
    % Calculate 2D bearing
    new_psi = atan2(new_v(2), new_v(1));
    aNc = aN(1,:);
    
elseif flags.space_dim == 3
    % Calculate declination and magnitude of normal acceleration (phi is
    % angle between the vector and a vertical plane.)
    [phi, aNc] = cart2pol(aN(1,:),aN(2,:));
    new_psi = 0;
    
end

%%% Solve planar differential equation %%%

% speed
aT = min(aT, (new_sdot-min_speed)/dt(end));
old_sdot = new_sdot - aT*dt;

% bearing
if abs(aT)>1E-10
    old_psi = new_psi - (aNc./aT).*log(new_sdot./old_sdot);
else
    old_psi = new_psi - (aNc.*dt)./old_sdot;
end

% displacement
SF = 4*aT.^2 + aNc.^2;

old_u = zeros(flags.space_dim,1);
if (aT~=0)&&(aNc~=0)
    old_u(1,:) = -((new_sdot.^2)./SF).*( aNc.*sin(new_psi)+2*aT.*cos(new_psi)) + ((old_sdot^2)./SF)*( aNc.*sin(old_psi)+2*aT.*cos(old_psi));
    old_u(2,:) = -((new_sdot.^2)./SF).*(-aNc.*cos(new_psi)+2*aT.*sin(new_psi)) + ((old_sdot^2)./SF)*(-aNc.*cos(old_psi)+2*aT.*sin(old_psi));
elseif (aT==0)&&(aNc~=0)
    old_u(1,:) = -((new_sdot.^2)./aNc).*( sin(new_psi) - sin(old_psi) );
    old_u(2,:) = -((new_sdot.^2)./aNc).*(-cos(new_psi) + cos(old_psi) );
elseif (aT~=0)&&(aNc==0)
    old_u(1,:) = -0.5*dt.*cos(old_psi).*new_sdot;
    old_u(2,:) = -0.5*dt.*sin(old_psi).*new_sdot;
else
    old_u(1,:) = -( old_sdot*dt.*cos(old_psi) );
    old_u(2,:) = -( old_sdot*dt.*sin(old_psi) );
end

if flags.space_dim == 3
    new_psi = -old_psi;
    old_psi = 0;
end

% Calculate cartesian in-plane velocity
old_udot = zeros(flags.space_dim,1);
[old_udot(1,:), old_udot(2,:)] = pol2cart(old_psi, old_sdot);

% Translate back to original coordinate system
if flags.space_dim == 2
    old_r = bsxfun(@plus, old_u, new_r);
    old_v = old_udot;
elseif flags.space_dim == 3
    
    % Calculate unit vectors for 3D intrinsics at start of sojourn
    ett = unit(new_v,1);
    dpsi = new_psi - old_psi;
    
    PC1 = sin(dpsi)*sin(phi);
    PC2 = (ett(2)^2 + ett(1)^2 - PC1^2);
    PC3 = sqrt(PC2);
    PC4 = cos(phi)*sin(dpsi);
    PC5 = (PC1^2 - 1);
    PC6 = (ett(1)^2 + ett(2)^2);

    et = zeros(3,1);
    et(1) = -(cos(dpsi)*(ett(1)*PC2 - ett(2)*PC1*PC3) + ett(3)*PC4*(ett(1)*PC3 - ett(2)*PC1))/(PC5*PC6);
    et(2) = -(cos(dpsi)*(ett(2)*PC2 + ett(1)*PC1*PC3) + ett(3)*PC4*(ett(2)*PC3 + ett(1)*PC1))/(PC5*PC6);
    et(3) = -( ett(3)*cos(dpsi) - PC4*PC3 )/PC5;
    
    en = zeros(3,1); u = sqrt(et(1)^2+et(2)^2);
    en(1,:) =  (et(2)*sin(phi)-et(1)*et(3)*cos(phi))/u;
    en(2,:) = (-et(1)*sin(phi)-et(2)*et(3)*cos(phi))/u;
    en(3,:) = cos(phi)*u;
    eb = cross(et,en);
    
    R = [et, en, eb];
    
    u_disp = -[[cos(dpsi) -sin(dpsi); sin(dpsi) cos(dpsi)]*(-old_u(1:2)); 0];
    
    old_r = R*u_disp + new_r;
    old_v = R*old_udot;

end

% Stack them up
old_x = [old_r; old_v];

% Check
assert(all(isreal(old_x)));

end

