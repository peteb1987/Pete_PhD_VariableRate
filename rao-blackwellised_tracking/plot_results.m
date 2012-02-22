function plot_results( flags, params, fig, cp_x, cp_tau, times, x, observs, pts, kd_times, kd_est )
%PLOT_TRACKING_RESULTS Plot tracks, observations and optionally overlay
%particles

% Select figure
figure(fig); clf;

% Make it big
screen_size = get(0, 'ScreenSize');
set(gcf, 'Position', [screen_size(3)/4 0 screen_size(3)/2 screen_size(4) ] );

% Select first axis for x,y(,z) plot
subplot(7,2,1:6), hold on

% Set plot area to contain the track (and the origin)
x1_max = max(max(x(1,:))+20, 10);
x2_max = max(max(x(2,:))+20, 10);
x1_min = min(min(x(1,:))-20, -10);
x2_min = min(min(x(2,:))-20, -10);
xlim([x1_min, x1_max]), ylim([x2_min, x2_max])

% Put an x on the origin
plot(0, 0, 'xk', 'markersize', 10);

% Plot track
plot(x(1,:), x(2,:), 'b', 'linewidth', 3);

% Plot state/jump points
plot(cp_x(1,:), cp_x(2,:), 'g*');

% Plot observations
if flags.obs_mod == 1
    plot(observs(1,:), observs(2,:), 'r');
elseif flags.obs_mod == 2
    [x1, x2] = pol2cart(observs(1,:), observs(2,:));
    plot(x1, x2, 'r');
else
    error('unhandled option');
end

% Plot velocities on separate axes
subplot(7,2,7), hold on
plot(times, x(3,:), 'b', 'linewidth', 3)
ylabel('x velocity')

subplot(7,2,9), hold on
plot(times, x(4,:), 'b', 'linewidth', 3)
ylabel('y velocity')

% if flags.obs_vel
%     subplot(7,2,8)
%     plot(times, observs(flags.space_dim+1,:), 'r')
%     title('Bearing Rate')
%     
%     subplot(7,2,10)
%     plot(times, observs(flags.space_dim+2,:), 'r')
%     if flags.space_dim == 3
%         title('Elevation Rate')
%     elseif flags.space_dim == 2
%         title('Range Rate')
%     end
% 
% end

if ~isempty(pts)
    % Get interpolated states from particles
    pts_mu = cat(3,pts.mu);
    K = size(pts_mu, 2);
    x1 = squeeze(pts_mu(1,:,:));
    x2 = squeeze(pts_mu(2,:,:));
%     x1dot = squeeze(pts_mu(flags.space_dim+1,:,:));
%     x2dot = squeeze(pts_mu(flags.space_dim+2,:,:));
    
    % Overlay particles
    subplot(7,2,1:6), hold on
    plot(x1, x2);
    
%     subplot(7,2,7), hold on
%     plot(times(1:K), x1dot);
%     
%     subplot(7,2,9), hold on
%     plot(times(1:K), x2dot);
    
end

% Plot kernel density estimate
if ~isempty(kd_est)
    subplot(7,2,[13,14]), hold on
    plot(kd_times, kd_est, 'b');
    if flags.gen_data, for tt=1:length(cp_tau), plot([cp_tau(tt),cp_tau(tt)], [0,1]','r'); end, end
end

% Draw
drawnow; shg;


end

