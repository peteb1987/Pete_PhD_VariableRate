function plot_results( flags, params, fig, x, tau, times, intx, observs, pts, kd_times, kd_est )
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
x1_max = max(max(intx(1,:))+20, 10);
x2_max = max(max(intx(2,:))+20, 10);
x1_min = min(min(intx(1,:))-20, -10);
x2_min = min(min(intx(2,:))-20, -10);
xlim([x1_min, x1_max]), ylim([x2_min, x2_max])

% Put an x on the origin
plot(0, 0, 'xk', 'markersize', 10);

if flags.space_dim == 2
    % Plot track
    plot(intx(1,:), intx(2,:), 'b', 'linewidth', 3);
    
    if ~isempty(x)
        % Plot state/jump points
        plot(x(1,:), x(2,:), 'g*');
    end

    % Plot observations
    if flags.obs_mod == 1
        plot3(observs(1,:), observs(2,:), 'r');
    elseif flags.obs_mod == 2
        [x1, x2] = pol2cart(observs(1,:), observs(2,:));
        plot(x1, x2, 'r');
    else
        error('unhandled option');
    end
    
elseif flags.space_dim == 3
    % Plot track
    plot3(intx(1,:), intx(2,:), intx(3,:), 'b', 'linewidth', 3);
    
    if ~isempty(x)
        % Plot state/jump points
        plot3(x(1,:), x(2,:), x(3,:), 'g*');
    end
    
    % Plot observations
    if flags.obs_mod == 1
        plot3(observs(1,:), observs(2,:), observs(3,:), 'r');
    elseif flags.obs_mod == 2
        [x1, x2, x3] = sph2cart(observs(1,:), observs(2,:), observs(3,:));
        plot3(x1, x2, x3, 'r');
    else
        error('unhandled option');
    end
else
    error('unhandled option');
end

% Plot velocities on separate axes
subplot(7,2,7), hold on
plot(times, intx(flags.space_dim+1,:), 'b', 'linewidth', 3)
ylabel('x velocity')

subplot(7,2,9), hold on
plot(times, intx(flags.space_dim+2,:), 'b', 'linewidth', 3)
ylabel('y velocity')

if flags.space_dim == 3
    subplot(7,2,11), hold on
    plot(times, intx(flags.space_dim+3,:), 'b', 'linewidth', 3)
    ylabel('z velocity')
end

if flags.obs_vel
    subplot(7,2,8)
    plot(times, observs(flags.space_dim+1,:), 'r')
    title('Bearing Rate')
    
    subplot(7,2,10)
    plot(times, observs(flags.space_dim+2,:), 'r')
    if flags.space_dim == 3
        title('Elevation Rate')
    elseif flags.space_dim == 2
        title('Range Rate')
    end
    
    if flags.space_dim == 3
        subplot(7,2,12)
        plot(times, observs(flags.space_dim+3,:), 'r')
        title('Range Rate')
    end
end

if ~isempty(pts)
    Np = length(pts);
    
    % Get interpolated states from particles
    pts_intx = cat(3,pts.intx);
    K = size(pts_intx, 2);
    x1 = squeeze(pts_intx(1,:,:));
    x2 = squeeze(pts_intx(2,:,:));
    x1dot = squeeze(pts_intx(flags.space_dim+1,:,:));
    x2dot = squeeze(pts_intx(flags.space_dim+2,:,:));
    if flags.space_dim == 3
        x3 = squeeze(pts_intx(3,:,:));
        x3dot = squeeze(pts_intx(flags.space_dim+3,:,:));
    end
    
    if flags.obs_vel
        % Accumulate predicted mean observations
        pred_obs = zeros(Np, params.obs_dim, K);
        for ii = 1:Np
            cpi = 1;
            for k = 1:K
                if (cpi<pts(ii).Ns)&&(times(k)>pts(ii).tau(cpi+1))
                    cpi = cpi + 1;
                end
                pred_obs(ii,:,k) = observation_mean(flags, params, pts(ii).intx(:,k), pts(ii).w(:,cpi));
            end
        end
        pred_range_rate = squeeze(pred_obs(:,flags.space_dim+2,:));
    end
    
    % Overlay particles
    subplot(7,2,1:6), hold on
    if flags.space_dim == 2
        plot(x1, x2);
    else
        plot3(x1, x2, x3);
    end
    
    subplot(7,2,7), hold on
    plot(times(1:K), x1dot);
    
    subplot(7,2,9), hold on
    plot(times(1:K), x2dot);
    
    if flags.obs_vel
        subplot(7,2,10), hold on
        plot(times(1:K), pred_range_rate);
    end
    
    if flags.space_dim == 3
        subplot(7,2,11), hold on
        plot(times(1:K), x3dot);
    end
    
end

% Plot kernel density estimate
if ~isempty(kd_est)
    subplot(7,2,[13,14]), hold on
    plot(kd_times, kd_est, 'b');
    if flags.gen_data, for tt=1:length(tau), plot([tau(tt),tau(tt)], [0,1]','r'); end, end
end

% Draw
drawnow; shg;


end

