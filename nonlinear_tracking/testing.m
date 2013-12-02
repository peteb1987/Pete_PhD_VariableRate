%%

aT_range = -10:0.1:10;
aN_range = -10:0.1:10;
lhood=zeros(length(aT_range), length(aN_range)); 

figure, hold on

for aaT = 1:length(aT_range)
    for aaN = 1:length(aN_range)
        for kk = 1:length(time)
            x = tracking_evaluatestate(model, cp_time, [aT_range(aaT);aN_range(aaN)], cp_state, time(kk));
            lhood(aaT,aaN) = lhood(aaT,aaN) + tracking_likelihood(model, cp_time, [aT_range(aaT);aN_range(aaN)], cp_state, time(kk), observ(:,kk));
        end
    end
end

figure, surf(aT_range, aN_range, lhood); shading interp;

%%

a_range = -5:1:5;
lhood=zeros(length(a_range)); 

figure, hold on

for aa1 = 1:length(a_range)
    for aa2 = 1:length(a_range)
        traj = zeros(4,100);
        for kk = 1:100
            traj(:,kk) = tracking_evaluatestate(model, cp_time, [a_range(aa1);a_range(aa2)], cp_state, cp_time+0.02*kk);
        end
        plot(traj(1,:), traj(2,:));
%         lhood(aa1,aa2)=neg_log_lhood(model, time, observ, [], [], [], cp_time, [a_range(aa1);a_range(aa2)], cp_state);
    end
end

%%

close all

a_range = -2:0.01:2;
state_aT = zeros(4,length(a_range));
state_aN = zeros(4,length(a_range));

psi_test = zeros(1,length(a_range));

cp_speed = sqrt(cp_state(3)^2+cp_state(4)^2);
cp_head = atan2(cp_state(4),cp_state(3));

for aa = 1:length(a_range)
    state_aT(:,aa) = tracking_evaluatestate(model, cp_time, [a_range(aa);cp_param(2)], cp_state, time(end));
    psi_test(aa) = cp_head + cp_param(2)*(time(end)-cp_time)/cp_speed - 0.5*a_range(aa)*cp_param(2)*((time(end)-cp_time)/cp_speed)^2;
end
for aa = 1:length(a_range)
    state_aN(:,aa) = tracking_evaluatestate(model, cp_time, [cp_param(1);a_range(aa)], cp_state, time(end));
end

% speed
figure, plot(a_range, sqrt(state_aT(3,:).^2+state_aT(4,:).^2) );
figure, plot(a_range, sqrt(state_aN(3,:).^2+state_aN(4,:).^2) );

% heading
figure, plot(a_range, atan2(state_aT(4,:),state_aT(3,:)) );
figure, plot(a_range, atan2(state_aN(4,:),state_aN(3,:)) );

%%
delt = 1E-6;

% aT
s_delt = s0+(aT+delt)*(t-t0);
psi_delt = psi0 + (aN/(aT+delt))*log(s_delt/s0);
alpha_delt = [2*(aT+delt); aN]/(4*(aT+delt)^2+aN^2);
M_delt = [cos(psi_delt), sin(psi_delt); sin(psi_delt), -cos(psi_delt)];
r_delt = r0 + ((s_delt^2)*M_delt - (s0^2)*M0)*alpha_delt;

% (s_delt-s)/delt
% (psi_delt-psi)/delt
% (alpha_delt-alpha)/delt
(r_delt-r)/delt

% aN
s_delt = s0+aT*(t-t0);
psi_delt = psi0 + ((aN+delt)/aT)*log(s_delt/s0);
alpha_delt = [2*aT; (aN+delt)]/(4*aT^2+(aN+delt)^2);
M_delt = [cos(psi_delt), sin(psi_delt); sin(psi_delt), -cos(psi_delt)];
r_delt = r0 + ((s_delt^2)*M_delt - (s0^2)*M0)*alpha_delt;

% (s_delt-s)/delt
% (psi_delt-psi)/delt
% (alpha_delt-alpha)/delt
(r_delt-r)/delt

%%

% Plot data
[ox, oy] = pol2cart(observ(1,:), observ(2,:));
figure, hold on
plot(true_state(1,:), true_state(2,:), 'k--');
plot(ox, oy, 'r.')
plot(cp_state(1,:), cp_state(2,:), 'g*')
    
% Plot trajectories
for ii = 1:algo.N
    %     plot(pf(ii).cp_state(1,:), pf(ii).cp_state(2,:));
    traj = zeros(model.ds,kk);
    for cc = 1:kk+Lt
        cp_ind = most_recent_changepoint(pf(ii).cp_time, time(cc));
        pre_cp_time = pf(ii).cp_time(cp_ind);
        pre_cp_param = pf(ii).cp_param(:,cp_ind);
        pre_cp_state = pf(ii).cp_state(:,cp_ind);
        traj(:,cc) = tracking_evaluatestate(model, pre_cp_time, pre_cp_param, pre_cp_state, time(cc));
    end
    plot(traj(1,:), traj(2,:))
end