function [ x ] = tracking_calc_next_state_batch_accel_3D( flags, last_x, dt, w )
%TRACKING_CALC_NEXT_STATE_BATCH_ACCEL_3D Batch version - handles a matrix
%for w, with a different set of accelerations in each column

Ns = size(w,2);

% All hell breaks lose if s hit 0, so lets limit it!
min_speed = 0.5;

% Get accelerations
aT = w(1,:);
aNv = w(2,:);
aNh = w(3,:);
aN = sqrt(aNv.^2+aNh.^2);
if flags.dyn_mod==6
    aX = w(4:6,:);
else
    aX = zeros(3,Ns);
end

% Get old state
last_r = last_x(1:3,:);
last_v = last_x(4:6,:);
sdot = sqrt(sum(last_v.^2));
psi = 0;

% Find angle of normal acceleration relative to vertical plane
phi = atan(aNv./aNh);
phi(isnan(phi))=0;

% Find unit vectors and 3D rotation matrix from plane to 3D space
et = repmat( (last_v/sdot) ,1,Ns);
eb = zeros(3,Ns); u = sqrt(et(1)^2+et(2)^2);
eb(1,:) =  (et(2)*sin(phi)-et(1)*et(3)*cos(phi))/u;
eb(2,:) = (-et(1)*sin(phi)-et(2)*et(3)*cos(phi))/u;
eb(3,:) = cos(phi)*u;
en = cross(eb,et);
R = squeeze(num2cell(permute(cat(3, et, en, eb),[1,3,2]),[1,2]));
% R = cat(3, et', en', eb');

%%% In-plane modelling - start with zero bearing and sdot speed in the _i_ direction

% Calculate new speed, and handle lower limiting
new_sdot = sdot + aT*dt;
new_sdot(new_sdot<min_speed) = min_speed;
aT = (new_sdot-sdot)/dt;

% Flags for where aT and aN are zero
waTz = (aT==0);             % flag for where aT is 0
waNz = (aN==0);             % flag for where aN is 0
bnz = ((~waTz)&(~waNz));   % both not zero
oaTz = ((waTz)&(~waNz));   % only aT zero
oaNz = ((~waTz)&(waNz));   % only aN zero
bz = ((waTz)&(waNz));      % both zero

% Calculate new bearing
new_psi(~waTz) = psi + (aN(~waTz)./aT(~waTz)).*log(new_sdot(~waTz)/sdot);
new_psi(waTz) = psi + (aN(waTz).*dt)./sdot;

% Calculate new u (within plane position)
u = zeros(3,Ns);
SF1 = 4*aT.^2 + aN.^2;

% Calculate new x1 and x2 when...

% Neither aT nor aN are zero
u(1,bnz) = ((new_sdot(bnz).^2)./SF1(bnz)).*( aN(bnz).*sin(new_psi(bnz))+2*aT(bnz).*cos(new_psi(bnz))) - ((sdot.^2)./SF1(bnz)).*( aN(bnz).*sin(psi)+2*aT(bnz).*cos(psi));
u(2,bnz) = ((new_sdot(bnz).^2)./SF1(bnz)).*(-aN(bnz).*cos(new_psi(bnz))+2*aT(bnz).*sin(new_psi(bnz))) - ((sdot.^2)./SF1(bnz)).*(-aN(bnz).*cos(psi)+2*aT(bnz).*sin(psi));

% aT is zero but aN isn't
u(1,oaTz) = ((new_sdot(oaTz).^2)./aN(oaTz)).*( sin(new_psi(oaTz)) - sin(psi) );
u(2,oaTz) = ((new_sdot(oaTz).^2)./aN(oaTz)).*(-cos(new_psi(oaTz)) + cos(psi) );

% aN is zero but aT isn't
u(1,oaNz) = 0.5*dt*cos(psi).*new_sdot(oaNz);
u(2,oaNz) = 0.5*dt*sin(psi).*new_sdot(oaNz);

% Both aT and aN are zero
u(1,bz) = ( sdot*dt.*cos(psi) );
u(2,bz) = ( sdot*dt.*sin(psi) );


% Calculate cartesian in-plane velocity
udot = zeros(3,Ns);
[udot(1,:), udot(2,:)] = pol2cart(new_psi, new_sdot);

% Translate position and velocity into 3D space
u_cell = num2cell(u, 1)';
udot_cell = num2cell(udot, 1)';
aX_cell = num2cell(aX*dt, 1)';
r = cellfun(@(rot, u_ind, aX_ind) {rot*u_ind+last_r+aX_ind}, R, u_cell, aX_cell) ;
v = cellfun(@(rot, udot_ind) {rot*udot_ind}, R, udot_cell);

% Stack them up
x = [cell2mat(r'); cell2mat(v')];




end

