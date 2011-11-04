figure(1)
% scrsz = get(0, 'ScreenSize');
% h=figure('Position', [1 scrsz(4) scrsz(3) scrsz(4)]);
% drawnow; pause(10);

parts = [];

for kk=1:params.K
	clf
    hold on, xlim([-200, 200]), ylim([-200, 200])
    plot(interp_state(1,:), interp_state(2,:), 'b');
    plot(state(1,:), state(2,:), 'g*');
    plot(observ(1,:), observ(2,:), 'r');
    parts = plot(filt_part_sets{kk}.pts_intx(:,1:kk,1)', filt_part_sets{kk}.pts_intx(:,1:kk,2)');
    pause(0.1)
end