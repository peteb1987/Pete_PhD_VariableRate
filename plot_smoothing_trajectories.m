clf

plot_tracking_data;
    
pts_intx = cat(3,smooth_pts.intx);
x1 = squeeze(pts_intx(1,:,:));
x2 = squeeze(pts_intx(2,:,:));
x3 = squeeze(pts_intx(3,:,:));
x4 = squeeze(pts_intx(4,:,:));

subplot(6,2,1:6)
plot(x1, x2);

subplot(6,2,7)
plot(times, x3);

subplot(6,2,9)
plot(times, x4);