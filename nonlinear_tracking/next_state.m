function [ new_x ] = next_state( flags, params, old_x, w, dt )
%NEXT_STATE Deterministically calculate the next state given the previous
%state, the random variable, and the time difference using the various
%dynamic models

% Two types of batch operation can be handled: time batch, for which dt is
% a row vector, and acceleration batch, for which w is a matrix of
% accelerations in columns.

K = numel(dt);
Ns = size(w,2);
assert((Ns==1)||(K==1), 'Cannot handle batches in both time and acceleration');
no_cols = max(Ns, K);

% If no time has passed (i.e. we're looking at an observation exactly
% after a jump) then return straight away. (This happens at t=0)
if (K==1)&&(dt<=0)
    new_x=old_x;
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

% Get old state
old_r = old_x(1:flags.space_dim,:);
old_v = old_x(flags.space_dim+1:end,:);

% Transform to planar intrinisics
old_sdot = norm(old_v);

if flags.space_dim == 2    
    % Calculate 2D bearing
    old_psi = atan(old_v(2)/old_v(1));
    aNc = aN(1,:);
    
elseif flags.space_dim == 3
    % Calculate declination and magnitude of normal acceleration (phi is
    % angle between the vector and a vertical plane.)
    [phi, aNc] = cart2pol(aN(2,:),-aN(1,:));
    
    % Calculate unit vectors for 3D intrinsics at start of sojourn
    et = unit(old_v);
    eb = zeros(3,Ns); u = sqrt(et(1,:)^2+et(2,:)^2);
    eb(1,:) =  (et(2,:)*sin(phi)-et(1,:)*et(3,:)*cos(phi))/u;
    eb(2,:) = (-et(1,:)*sin(phi)-et(2,:)*et(3,:)*cos(phi))/u;
    eb(3,:) = cos(phi)*u;
    en = cross(eb,et);
    R = [et, en, eb];
    
    old_psi = 0;
    
end

%%% Solve planar differential equation %%%

% speed
aT = max(aT, min_speed-old_sdot/dt(end));
new_sdot = old_sdot + aT*dt;

% bearing
if aT~=0
    new_psi = old_psi + (aNc./aT).*log(new_sdot./old_sdot);
else
    new_psi = old_psi + (aNc.*dt)./old_sdot;
end

% displacement
SF = 4*aT.^2 + aNc.^2;
new_u = zeros(flags.space_dim,no_cols);
if (aT~=0)&&(aNc~=0)
    new_u(1,:) = ((new_sdot.^2)./SF).*( aNc.*sin(new_psi)+2*aT.*cos(new_psi)) - ((old_sdot^2)./SF)*( aNc.*sin(old_psi)+2*aT.*cos(old_psi));
    new_u(2,:) = ((new_sdot.^2)./SF).*(-aNc.*cos(new_psi)+2*aT.*sin(new_psi)) - ((old_sdot^2)./SF)*(-aNc.*cos(old_psi)+2*aT.*sin(old_psi));
elseif (aT==0)&&(aNc~=0)
    new_u(1,:) = ((new_sdot.^2)./aNc).*( sin(new_psi) - sin(old_psi) );
    new_u(2,:) = ((new_sdot.^2)./aNc).*(-cos(new_psi) + cos(old_psi) );
elseif (aT~=0)&&(aNc==0)
    new_u(1,:) = 0.5*dt.*cos(old_psi).*new_sdot;
    new_u(2,:) = 0.5*dt.*sin(old_psi).*new_sdot;
else
    new_u(1,:) = ( old_sdot*dt.*cos(old_psi) );
    new_u(2,:) = ( old_sdot*dt.*sin(old_psi) );
end

% Calculate cartesian in-plane velocity
new_udot = zeros(3,no_cols);
[new_udot(1,:), new_udot(2,:)] = pol2cart(new_psi, new_sdot);

% Translate back to original coordinate system
if flags.space_dim == 2
    new_r = bsxfun(@plus, new_u, old_r);
    new_v = new_udot;
elseif flags.space_dim == 3
    new_r = bsxfun(@plus, R*new_u + aX*dt, old_r);
    new_v = R*new_udot;
end

% Stack them up
new_x = [new_r; new_v];

% Check
assert(all(isreal(new_x)));

end

