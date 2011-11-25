function [ x ] = tracking_calc_next_state_batch_accel( flags, last_x, dt, w )
%TRACKING_CALC_NEXT_STATE Calculate the next state given the dtrevious one,
%the resdtective times, and the random variables

if flags.dyn_mod >= 5
    [ x ] = tracking_calc_next_state_batch_accel_3D( flags, last_x, dt, w );
    return
end

% This version takes a set of w vectors (as columns of a matrix) and does a
% calculation for each. Useful for sigma points.

Ns = size(w,2);

% If no times has passed (i.e. we're looking at an observation exactly
% after a jump) then return straight away. (This happens at t=0)
if dt==0
    x=last_x;
    return
end

% All hell breaks lose if s hit 0, so lets limit it!
min_speed = 0.5;

% Get accelerations
aT = w(1,:);
aP = w(2,:);
aX1 = zeros(1,Ns);
aX2 = zeros(1,Ns);
aS = zeros(1,Ns);
aB = zeros(1,Ns);
if flags.dyn_mod == 2
    aX1 = w(3,:);
    aX2 = w(4,:);
elseif (flags.dyn_mod == 3)||(flags.dyn_mod == 4)
    aB = w(3,:);
    aS = w(4,:);
end

% Get old state
psi = last_x(3);
sdot = last_x(4);
x1 = last_x(1);
x2 = last_x(2);

% Create an array for the new state
x = zeros(size(w));
new_sdot = zeros(1,Ns);
new_psi = zeros(1,Ns);
new_x1 = zeros(1,Ns);
new_x2 = zeros(1,Ns);

% Calculate new speed, and handle lower limiting
new_sdot = sdot + aT*dt;
new_sdot(new_sdot<min_speed) = min_speed;
aT = (new_sdot-sdot)/dt;

% Precalculate scale factors
SF1 = 4*aT.^2 + aP.^2;

% Flags for where aT and aP are zero
waTz = (aT==0);             % flag for where aT is 0
waPz = (aP==0);             % flag for where aP is 0
bnz = ((~waTz)&(~waPz));   % both not zero
oaTz = ((waTz)&(~waPz));   % only aT zero
oaPz = ((~waTz)&(waPz));   % only aP zero
bz = ((waTz)&(waPz));      % both zero

% Calculate new bearing
new_psi(~waTz) = psi + (aP(~waTz)./aT(~waTz)).*log(new_sdot(~waTz)/sdot);
new_psi(waTz) = psi + (aP(waTz).*dt)./sdot;

if flags.dyn_mod == 4
    psi = psi + aB;
    sdot = sdot + aS;
    new_sdot = new_sdot + aS;
    new_psi = new_psi + aB;
else
    psi = repmat(psi, 1, Ns);
    sdot = repmat(sdot, 1, Ns);
end

% Calculate new x1 and x2 when...

% Neither aT nor aP are zero
new_x1(bnz) = x1 + ((new_sdot(bnz).^2)./SF1(bnz)).*( aP(bnz).*sin(new_psi(bnz))+2*aT(bnz).*cos(new_psi(bnz))) - ((sdot(bnz).^2)./SF1(bnz)).*( aP(bnz).*sin(psi(bnz))+2*aT(bnz).*cos(psi(bnz))) + aX1(bnz)*dt;
new_x2(bnz) = x2 + ((new_sdot(bnz).^2)./SF1(bnz)).*(-aP(bnz).*cos(new_psi(bnz))+2*aT(bnz).*sin(new_psi(bnz))) - ((sdot(bnz).^2)./SF1(bnz)).*(-aP(bnz).*cos(psi(bnz))+2*aT(bnz).*sin(psi(bnz))) + aX2(bnz)*dt;

% aT is zero but aP isn't
new_x1(oaTz) = x1 + ((new_sdot(oaTz).^2)./aP(oaTz)).*( sin(new_psi(oaTz)) - sin(psi(oaTz)) ) + aX1(oaTz)*dt;
new_x2(oaTz) = x2 + ((new_sdot(oaTz).^2)./aP(oaTz)).*(-cos(new_psi(oaTz)) + cos(psi(oaTz)) ) + aX2(oaTz)*dt;

% ap is zero but aT isn't
new_x1(oaPz) = x1 + 0.5*dt*cos(psi(oaPz)).*new_sdot(oaPz) + aX1(oaPz)*dt;
new_x2(oaPz) = x2 + 0.5*dt*sin(psi(oaPz)).*new_sdot(oaPz) + aX2(oaPz)*dt;

% Both aT and aP are zero
new_x1(bz) = x1 + ( sdot(bz)*dt.*cos(psi(bz)) ) + aX1(bz)*dt;
new_x2(bz) = x2 + ( sdot(bz)*dt.*sin(psi(bz)) ) + aX2(bz)*dt;

% Add extra noise to sdot and psi
if flags.dyn_mod == 3
    new_sdot = new_sdot + aS*dt;
    new_sdot(new_sdot<min_speed) = min_speed;
    new_psi = new_psi + aB*dt;
end

x(1,:) = new_x1;
x(2,:) = new_x2;
x(3,:) = new_psi;
x(4,:) = new_sdot;

assert(all(isreal(x(:))));

end

