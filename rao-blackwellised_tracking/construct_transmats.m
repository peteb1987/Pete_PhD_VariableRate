function [A, Q] = construct_transmats(dt, model, w, sigma)
%CONSTRUCT_TRANSMATS Construct the transition matrices for different
%tracking models

if (model == 1)||(w==0)
    
    A = [1, 0, dt,  0, dt^2/2,      0;
         0, 1,  0, dt,      0, dt^2/2;
         0, 0,  1,  0,     dt,      0;
         0, 0,  0,  1,      0,     dt;
         0, 0,  0,  0,      1,      0;
         0, 0,  0,  0,      0,      1];
    
    Q = [ dt^5/20,       0, dt^4/8,      0, dt^3/6,      0;
                0, dt^5/20,      0, dt^4/8,      0, dt^3/6;
           dt^4/8,       0, dt^3/3,      0, dt^2/2,      0;
                0,  dt^4/8,      0, dt^3/3,      0, dt^2/2;
           dt^3/6,       0, dt^2/2,      0,     dt,      0;
                0,  dt^3/6,      0, dt^2/2,      0,     dt] * sigma;

    
    
elseif model == 2
    
    A = [1, 0, dt,  0,     (1-cos(w*dt))/w^2, -(w*dt-sin(w*dt))/w^2;
         0, 1,  0, dt,  (w*dt-sin(w*dt))/w^2,     (1-cos(w*dt))/w^2;
         0, 0,  1,  0,           sin(w*dt)/w,      -(1-cos(w*dt))/w;
         0, 0,  0,  1,       (1-cos(w*dt))/w,           sin(w*dt)/w;
         0, 0,  0,  0,             cos(w*dt),            -sin(w*dt);
         0, 0,  0,  0,             sin(w*dt),             cos(w*dt)];
    
    temp1 = ( (w*dt)^3/3 + 2*w*dt - 4*sin(w*dt) + 2*w*dt*cos(w*dt) )/w^5;
    Q1 = [temp1, 0; 0, temp1];
    
    temp1 = ( 1 + (w*dt)^2/2 - cos(w*dt) - w*dt*sin(w*dt) )/w^4;
    temp2 = ( 2*w*dt - 3*sin(w*dt) + w*dt*cos(w*dt) )/w^4;
    Q2 = [temp1, -temp2; temp2, temp1];
    
    temp1 = ( -w*dt + 2*sin(w*dt) - w*dt*cos(w*dt) )/w^3;
    temp2 = ( 2 - 2*cos(w*dt) - w*dt*sin(w*dt) )/w^3;
    Q3 = [temp1, -temp2; temp2, temp1];
    
    temp1 = ( 2*w*dt - 2*sin(w*dt) )/w^3;
    Q4 = [temp1, 0; 0, temp1];
    
    temp1 = ( 1 - cos(w*dt) )/w^2;
    temp2 = ( w*dt - sin(w*dt) )/w^2;
    Q5 = [temp1, -temp2; temp2, temp1];

    Q6 = [dt 0; 0 dt];
    
    Q = [Q1  Q2' Q3';
         Q2  Q4  Q5';
         Q3  Q5  Q6 ] * sigma;

% Q = (Q+Q')/2;

%     A = [1, 0, dt,  0,    dt^2/2, -w*dt^3/6;
%          0, 1,  0, dt,  w*dt^3/6,    dt^2/2;
%          0, 0,  1,  0,        dt,  w*dt^2/2;
%          0, 0,  0,  1, -w*dt^2/2,        dt;
%          0, 0,  0,  0,         1,     -w*dt;
%          0, 0,  0,  0,      w*dt,         1];
%     
%     Q = [  dt^5/20,          0,     dt^4/8, w*dt^5/60,     dt^3/6, w*dt^4/12;
%                  0,    dt^5/20, -w*dt^5/60,    dt^4/8, -w*dt^4/12,    dt^3/6;
%             dt^4/8, -w*dt^5/60,     dt^3/3,         0,     dt^2/2,  w*dt^3/6;
%          w*dt^5/60,     dt^4/8,          0,    dt^3/3,  -w*dt^3/6,    dt^2/2;
%             dt^3/6, -w*dt^4/12,     dt^2/2, -w*dt^3/6,         dt,         0;
%          w*dt^4/12,     dt^3/6,   w*dt^3/6,    dt^2/2,          0,        dt] * sigma;

end


end

