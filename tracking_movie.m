figure(1)
clf

plot_tracking_data;

parts_xy = [];
parts_b = [];
parts_s = [];

for kk=1:params.K
    
    pts_intx = cat(3,filt_part_sets{kk}.intx);
    x1 = squeeze(pts_intx(1,:,:));
    x2 = squeeze(pts_intx(2,:,:));
    x3 = squeeze(pts_intx(3,:,:));
    x4 = squeeze(pts_intx(4,:,:));
    
    subplot(6,2,1:6)
    delete(parts_xy)
    parts_xy = plot(x1, x2);
    
    subplot(6,2,7)
    delete(parts_b)
    parts_b = plot(times(1:kk), x3);

    subplot(6,2,9)
    delete(parts_s)
    parts_s = plot(times(1:kk), x4);
    
    pause(0.05)
end