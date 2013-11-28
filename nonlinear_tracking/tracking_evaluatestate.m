function [ state ] = tracking_evaluatestate(model, cp_time, cp_param, cp_state, time)
%TRACKING_EVALUATESTATE Evaluate the state given the preceeding changepoint

% Convert state to a intrinsic form
cp_pos = cp_state(1:2);
cp_head = atan2(cp_state(4),cp_state(3));
cp_speed = sqrt(cp_state(3)^2+cp_state(4)^2);
cp_head_mat = [cos(cp_head), sin(cp_head); sin(cp_head), -cos(cp_head)];

% Get accelerations
aT = cp_param(1);
aN = cp_param(2);
anorm = [2*aT; aN]/(4*aT^2 + aN^2);

% Speed update
speed = cp_speed + aT*(time-cp_time);
if speed < model.vmin
    speed = model.vmin;
    aT = (speed - cp_speed)/(time - cp_time);
end

% Remaining update equations
if aT > 0
    head = cp_head + (aN/aT)*log(speed/cp_speed);
else
    head = cp_head + aN*(time - cp_time)/speed;
end
head_mat = [cos(head), sin(head); sin(head), -cos(head)];
pos = cp_pos + (speed^2*head_mat - cp_speed^2*cp_head_mat)*anorm;

% Put it back together
state = [pos; speed*cos(head); speed*sin(head)];

end
