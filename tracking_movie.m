figure(1)
clf
% scrsz = get(0, 'ScreenSize');
% h=figure('Position', [1 scrsz(4) scrsz(3) scrsz(4)]);
% drawnow; pause(10);

parts = [];

hold on, xlim([-200, 200]), ylim([-200, 200])
plot(interp_state(1,:), interp_state(2,:), 'b');
plot(state(1,:), state(2,:), 'g*');
if flags.obs_mod == 1
    plot(observ(1,:), observ(2,:), 'r');
elseif flags.obs_mod == 2
    [x1, x2] = pol2cart(observ(1,:), observ(2,:));
    plot(x1, x2, 'r');
end

for kk=1:params.K
    delete(parts)
    parts = plot(filt_part_sets{kk}.pts_intx(:,1:kk,1)', filt_part_sets{kk}.pts_intx(:,1:kk,2)');
    pause(0.1)
end