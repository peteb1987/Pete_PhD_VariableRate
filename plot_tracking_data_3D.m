screen_size = get(0, 'ScreenSize');
set(gcf, 'Position', [screen_size(3)/4 0 screen_size(3)/2 screen_size(4) ] );

hold on;

subplot(7,2,1:6), hold on

% Set plot area to contain the track (and the origin)
x1_max = max(max(interp_state(1,:))+20, 10);
x2_max = max(max(interp_state(2,:))+20, 10);
x3_max = max(max(interp_state(3,:))+20, 10);
x1_min = min(min(interp_state(1,:))-20, -10);
x2_min = min(min(interp_state(2,:))-20, -10);
x3_min = min(min(interp_state(3,:))-20, -10);

% Put an x on the origin
plot(0, 0, 'xk', 'markersize', 10);

% Plot track
plot3(interp_state(1,:), interp_state(2,:), interp_state(3,:), 'b');
xlim([x1_min, x1_max]), ylim([x2_min, x2_max]), zlim([x3_min, x3_max]);

% Plot state/jump points
plot3(state(1,:), state(2,:), state(3,:), 'g*');

% Plot observations
if flags.obs_mod == 1
    plot3(observ(1,:), observ(2,:), observ(2,:), 'r');
elseif flags.obs_mod == 2
    [x1, x2, x3] = sph2cart(observ(1,:), observ(2,:), observ(3,:));
    plot3(x1, x2, x3, 'r');
end

% Plot bearing and speed
subplot(7,2,7), hold on
plot(times, interp_state(4,:), 'b', 'linewidth', 3)
ylabel('x velocity')

subplot(7,2,9), hold on
plot(times, interp_state(5,:), 'b', 'linewidth', 3)
ylabel('y velocity')

subplot(7,2,11), hold on
plot(times, interp_state(6,:), 'b', 'linewidth', 3)
ylabel('z velocity')

subplot(7,2,8)
if flags.obs_vel
    plot(times, observ(4,:), 'r')
    title('Bearing Rate')
end

subplot(7,2,10)
if flags.obs_vel
    plot(times, observ(5,:), 'r')
    title('Elevation Rate')
end

subplot(7,2,12)
if flags.obs_vel
    plot(times, observ(6,:), 'r')
    title('Range Rate')
end


% Plot kernel density estimate
if exist('filt_jump_kdest','var')
    subplot(6,2,[11,12]), hold on
    plot(filt_kd_times, filt_jump_kdest, 'b');
    if flags.gen_data, for tt=1:length(tau), plot([tau(tt),tau(tt)], [0,1]','r'); end, end
end

% Draw
drawnow; shg;