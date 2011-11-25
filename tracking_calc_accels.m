function w_new = tracking_calc_accels(flags, tau, next_tau, x, next_x)
%TRACKING_CALC_ACCELS Calculate the acceleration vector given two bounding
%states and times

if flags.dyn_mod >= 5
    w_new = tracking_calc_accels_3D(flags, tau, next_tau, x, next_x);
end

psi = x(3);
sdot = x(4);
x1 = x(1);
x2 = x(2);

new_psi = next_x(3);
new_sdot = next_x(4);
new_x1 = next_x(1);
new_x2 = next_x(2);

dt = next_tau - tau;

dx1 = new_x1 - x1;
dx2 = new_x2 - x2;
dpsi = new_psi - psi;

if flags.dyn_mod == 1
    
    aT = (new_sdot - sdot)/dt;
    aP = aT*dpsi/log(new_sdot/sdot);    
    
    w_new = [aT; aP];
    
    x_should = tracking_calc_next_state(flags, x, next_tau-tau, w_new);
    
    assert(x_should==next_x, 'Degenerate model, fool.');
    
elseif flags.dyn_mod == 2
    
    aT = (new_sdot - sdot)/dt;
    aP = aT*dpsi/log(new_sdot/sdot);
    SF1 = 4*aT.^2 + aP.^2;
    aX1 = ( new_x1 - x1 - (((new_sdot.^2)./SF1).*( aP.*sin(new_psi)+2*aT.*cos(new_psi)) - ((sdot^2)./SF1)*( aP.*sin(psi)+2*aT.*cos(psi))) )/dt;
    aX2 = ( new_x2 - x2 - (((new_sdot.^2)./SF1).*(-aP.*cos(new_psi)+2*aT.*sin(new_psi)) - ((sdot^2)./SF1)*(-aP.*cos(psi)+2*aT.*sin(psi))) )/dt;
    
    w_new = [aT; aP; aX1; aX2];
    
elseif flags.dyn_mod == 3
    
    error('No analytic solution');
    
elseif flags.dyn_mod == 4
    
    R = dx1/dx2;
    S = dx1^2+dx2^2;
    T = tan(psi);
    C = sqrt(sdot^4+new_sdot^4-2*sdot^2*new_sdot^2*cos(dpsi));
    A = (new_sdot^2*sin(dpsi)) / (new_sdot^2*cos(dpsi)-sdot^2);
    Q = (A*R+R*T+A*T-1) / (A*R*T+R+A+T);
    
    aT = 0.5*C/sqrt(S*(1+Q^2));
    aP = 2*Q*aT;
    aB = dpsi - (aP/aT)*log(new_sdot/sdot);
    aS = new_sdot - sdot - aT*dt;
    
    w_new = [aT; aP; aB; aS];
    
end

end

