figure(1)
clf

plot_tracking_data;

parts_xy = [];
parts_b = [];
parts_s = [];

for kk=1:params.K
    subplot(6,2,1:6)
    delete(parts_xy)
    parts_xy = plot(filt_part_sets{kk}.pts_intx(:,1:kk,1)', filt_part_sets{kk}.pts_intx(:,1:kk,2)');
    
    subplot(6,2,7)
    delete(parts_b)
    parts_b = plot(times(1:kk), filt_part_sets{kk}.pts_intx(:,1:kk,3)');

    subplot(6,2,9)
    delete(parts_s)
    parts_s = plot(times(1:kk), filt_part_sets{kk}.pts_intx(:,1:kk,4)');
    
    pause(0.05)
end