parts_xyz = [];
parts_xd = [];
parts_yd = [];
parts_zd = [];

for kk=params.K
    
    pts_intx = cat(3,filt_part_sets{kk}.intx);
    x1 = squeeze(pts_intx(1,:,:));
    x2 = squeeze(pts_intx(2,:,:));
    x3 = squeeze(pts_intx(3,:,:));
    x1dot = squeeze(pts_intx(4,:,:));
    x2dot = squeeze(pts_intx(5,:,:));
    x3dot = squeeze(pts_intx(6,:,:));
    
    subplot(7,2,1:6), hold on
    delete(parts_xyz)
    parts_xyz = plot3(x1, x2, x3);
    
    subplot(7,2,7), hold on
    delete(parts_xd)
    parts_xd = plot(times(1:kk), x1dot);

    subplot(7,2,9), hold on
    delete(parts_yd)
    parts_yd = plot(times(1:kk), x2dot);
    
    subplot(7,2,11), hold on
    delete(parts_zd)
    parts_zd = plot(times(1:kk), x3dot);
    
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
    
    pause(0.05)
end