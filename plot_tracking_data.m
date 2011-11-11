screen_size = get(0, 'ScreenSize');
set(gcf, 'Position', [screen_size(3)/4 0 screen_size(3)/2 screen_size(4) ] );

hold on;

subplot(6,2,1:6), hold on

% Set plot area to contain the track (and the origin)
x1_max = max(max(interp_state(1,:))+20, 10);
x2_max = max(max(interp_state(2,:))+20, 10);
x1_min = min(min(interp_state(1,:))-20, -10);
x2_min = min(min(interp_state(2,:))-20, -10);
xlim([x1_min, x1_max]), ylim([x2_min, x2_max])

% Put an x on the origin
plot(0, 0, 'xk', 'markersize', 10);

% Plot track
plot(interp_state(1,:), interp_state(2,:), 'b');

% Plot state/jump points
plot(state(1,:), state(2,:), 'g*');

% Plot observations
if flags.obs_mod == 1
    plot(observ(1,:), observ(2,:), 'r');
elseif flags.obs_mod == 2
    [x1, x2] = pol2cart(observ(1,:), observ(2,:));
    plot(x1, x2, 'r');
end

% Plot bearing and speed
subplot(6,2,7), hold on
plot(times, interp_state(3,:), 'b', 'linewidth', 3)
ylabel('Heading')
if params.obs_dim == 4
    subplot(6,2,8), hold on
    plot(times, observ(3,:), 'r');
    ylabel('Bearing Rate')
end

subplot(6,2,9), hold on
ylabel('Speed')
plot(times, interp_state(4,:), 'b', 'linewidth', 3)
if params.obs_dim == 4
    subplot(6,2,10), hold on
    plot(times, observ(4,:), 'r');
    ylabel('Range Rate')
end


% Plot kernel density estimate
if exist('filt_jump_kdest','var')
    subplot(6,2,[11,12]), hold on
    plot(filt_kd_times, filt_jump_kdest, 'b');
    if flags.gen_data, for tt=1:length(tau), plot([tau(tt),tau(tt)], [0,1]','r'); end, end
end

% Draw
drawnow; shg;