function filter_results_movie( flags, params, fig, x, tau, times, intx, observs, filt_pts )
%PLOT_TRACKING_RESULTS Plot tracks, observations and optionally overlay
%particles

plot_results( flags, params, fig, x, tau, times, intx, observs, [], [], [] )

xy_pts = [];
x1dot_pts = [];
x2dot_pts = [];
x3dot_pts = [];

for k = 1:params.K
    
    % Get interpolated states from particles
    pts_intx = cat(3,filt_pts{k}.intx);
    x1 = squeeze(pts_intx(1,:,:));
    x2 = squeeze(pts_intx(2,:,:));
    x1dot = squeeze(pts_intx(flags.space_dim+1,:,:));
    x2dot = squeeze(pts_intx(flags.space_dim+2,:,:));
    if flags.space_dim == 3
        x3 = squeeze(pts_intx(3,:,:));
        x3dot = squeeze(pts_intx(flags.space_dim+3,:,:));
    end
    
    % Overlay particles
    subplot(7,2,1:6), hold on
    delete(xy_pts)
    if flags.space_dim == 2
        xy_pts = plot(x1, x2);
    else
        xy_pts = plot3(x1, x2, x3);
    end
    
    subplot(7,2,7), hold on
    delete(x1dot_pts)
    x1dot_pts = plot(times(1:k), x1dot);
    
    subplot(7,2,9), hold on
    delete(x2dot_pts)
    x2dot_pts = plot(times(1:k), x2dot);
    
    if flags.space_dim == 3
        subplot(7,2,11), hold on
        delete(x3dot_pts)
        x3dot_pts = plot(times(1:k), x3dot);
    end
    
    % Draw
    drawnow; shg;
    pause(0.01);
    
end

end

