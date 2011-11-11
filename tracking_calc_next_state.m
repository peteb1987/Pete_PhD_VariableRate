function [ x ] = tracking_calc_next_state( flags, last_x, dt, w )
%TRACKING_CALC_NEXT_STATE Calculate the next state given the dtrevious one,
%the resdtective times, and the random variables

% If no times has passed (i.e. we're looking at an observation exactly
% after a jump) then return straight away. (This happens at t=0)
if dt==0
    x=last_x;
    return
end

% All hell breaks lose if s hit 0, so lets limit it!
min_speed = 0.5;

% Get accelerations
aT = w(1);
aP = w(2);
aX1 = 0;
aX2 = 0;
aS = 0;
aB = 0;
if flags.dyn_mod == 2
    aX1 = w(3);
    aX2 = w(4);
elseif flags.dyn_mod == 3
    aB = w(3);
    aS = w(4);
end

% Get old state
psi = last_x(3);
sdot = last_x(4);
x1 = last_x(1);
x2 = last_x(2);

% Create an array for the new state
x = zeros(size(last_x));

% Calculate new speed, and handle lower limiting
new_sdot = sdot + aT*dt;
new_sdot(new_sdot<min_speed) = min_speed;
aT = (new_sdot-sdot)/dt;

SF1 = 4*aT.^2 + aP.^2;

% Calculate new bearing
if aT~=0
    new_psi = psi + (aP./aT).*log(new_sdot./sdot);
else
    new_psi = psi + (aP.*dt)./sdot;
end

% Calculate new x1 and x2
if (aT~=0)&&(aP~=0)
    new_x1 = x1 + ((new_sdot.^2)./SF1).*( aP.*sin(new_psi)+2*aT.*cos(new_psi)) - ((sdot^2)./SF1)*( aP.*sin(psi)+2*aT.*cos(psi)) + aX1*dt;
    new_x2 = x2 + ((new_sdot.^2)./SF1).*(-aP.*cos(new_psi)+2*aT.*sin(new_psi)) - ((sdot^2)./SF1)*(-aP.*cos(psi)+2*aT.*sin(psi)) + aX2*dt;
elseif (aT==0)&&(aP~=0)
    new_x1 = x1 + ((new_sdot.^2)./aP).*( sin(new_psi) - sin(psi) ) + aX1*dt;
    new_x2 = x2 + ((new_sdot.^2)./aP).*(-cos(new_psi) + cos(psi) ) + aX2*dt;
elseif (aT~=0)&&(aP==0)
    new_x1 = x1 + 0.5*dt*cos(psi).*new_sdot + aX1*dt;
    new_x2 = x2 + 0.5*dt*sin(psi).*new_sdot + aX2*dt;
else
    new_x1 = x1 + ( sdot*dt.*cos(psi) ) + aX1*dt;
    new_x2 = x2 + ( sdot*dt.*sin(psi) ) + aX2*dt;
end

new_sdot = new_sdot + aS*dt;
new_psi = new_psi + aB*dt;

x(1) = new_x1;
x(2) = new_x2;
x(3) = new_psi;
x(4) = new_sdot;

assert(all(isreal(x)));

end

