function filter_results_movie2( flags, params, fig, x, tau, times, intx, observs, filt_pts )
%PLOT_TRACKING_RESULTS Plot tracks, observations and optionally overlay
%particles

% Select figure
figure(fig); clf; hold on;

% Make it big
screen_size = get(0, 'ScreenSize');
set(gcf, 'Position', [screen_size(3)/4 0 screen_size(3)/2 screen_size(4)/2 ] );

% Set plot area to contain the track (and the origin)
x1_max = max(max(intx(1,:))+20, 10);
x2_max = max(max(intx(2,:))+20, 10);
x1_min = min(min(intx(1,:))-20, -10);
x2_min = min(min(intx(2,:))-20, -10);
xlim([x1_min, x1_max]), ylim([x2_min, x2_max])

% Put an x on the origin
plot(0, 0, 'xk', 'markersize', 10);

% Plot track
plot(intx(1,:), intx(2,:), '-.k', 'linewidth', 2);

if ~isempty(x)
    % Plot state/jump points
    plot(x(1,:), x(2,:), 'g*');
end

% Plot observations
if flags.obs_mod == 1
    plot3(observs(1,:), observs(2,:), 'xr');
elseif flags.obs_mod == 2
    [x1, x2] = pol2cart(observs(1,:), observs(2,:));
    plot(x1, x2, 'xr');
else
    error('unhandled option');
end

xy_pts = [];
x1dot_pts = [];
x2dot_pts = [];
x3dot_pts = [];

vidObj = VideoWriter('vrpf_movie');
open(vidObj);

for k = 1:params.K
    
    % Get interpolated states from particles
    pts_intx = cat(3,filt_pts{k}.intx);
    x1 = squeeze(pts_intx(1,:,:));
    x2 = squeeze(pts_intx(2,:,:));
    x1dot = squeeze(pts_intx(flags.space_dim+1,:,:));
    x2dot = squeeze(pts_intx(flags.space_dim+2,:,:));
    
    % Overlay particles
    delete(xy_pts)
    xy_pts = plot(x1, x2, 'b');
 
    % Draw
    drawnow; shg;
    pause(0.001);
    
    currFrame = getframe;
    writeVideo(vidObj,currFrame);
    
    export_pdf(fig, ['animation/vrpf_' num2str(k) '.pdf'], 10, 10);
    
end

close(vidObj);

end

