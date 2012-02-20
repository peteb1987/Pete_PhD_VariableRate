close all

%% Data
fig = figure(1);
plot_results(flags, params, times, observ, tau, type, interp_state, [], [], fig);
subplot(3,1,1); xlabel('time'), %legend('observed value', 'value', 'Location', 'SouthEast')
subplot(3,1,2); xlabel('time'), %legend('trend', 'Location', 'SouthEast')
subplot(3,1,3); xlabel('time'), ylabel(''); %legend('trend jumps', 'value jumps');

wid = 6; hei = 5;
remove_whitespace_for_saveas;
print(fig, '-dpdf', 'example_data.pdf');

%% filter kd estimate
fig = figure(2); hold on;
plot(filt_kd_est.times, filt_kd_est.jump_1_kd, 'b');
plot(filt_kd_est.times, filt_kd_est.jump_2_kd, 'g');
for tt=1:length(tau)
    if type(tt) == 1
        plot([tau(tt),tau(tt)], [0,1.4]','k');
    elseif type(tt) == 2
        plot([tau(tt),tau(tt)], [0,1.4]','r');
    end
end
set(gca, 'FontSize', 8);
%legend('value jump kernel density', 'trend jump kernel density', 'trend jumps', 'value jumps', 'Location', 'SouthWest');

wid = 6; hei = 2;
remove_whitespace_for_saveas;
print(fig, '-dpdf', 'example_filter_kdest.pdf');

%% smoother kd estimate
fig = figure(3); hold on;
plot(smooth_kd_est.times, smooth_kd_est.jump_1_kd, 'b');
plot(smooth_kd_est.times, smooth_kd_est.jump_2_kd, 'g');
for tt=1:length(tau)
    if type(tt) == 1
        plot([tau(tt),tau(tt)], [0,1.4]','k');
    elseif type(tt) == 2
        plot([tau(tt),tau(tt)], [0,1.4]','r');
    end
end
set(gca, 'FontSize', 8);
%legend('value jump kernel density', 'trend jump kernel density', 'trend jumps', 'value jumps', 'Location', 'SouthWest');

wid = 6; hei = 2;
remove_whitespace_for_saveas;
print(fig, '-dpdf', 'example_smoother_kdest.pdf');


%% filter state estimate
fig = figure(4); hold on;
pts = filt_part_sets{end};
pts_intmu = cat(3, pts.intmu);
pts_intP = cat(4, pts.intP);

subplot(2,1,1), hold on
plot(times, squeeze(pts_intmu(1,:,:)));
plot(times, squeeze(pts_intmu(1,:,:)) + 2* sqrt(squeeze(pts_intP(1,1,:,:))), '-.');
plot(times, squeeze(pts_intmu(1,:,:)) - 2* sqrt(squeeze(pts_intP(1,1,:,:))), '-.');
plot(times, interp_state(1,:), 'b', 'LineWidth', 2);
xlabel('time'); ylabel('value');
ylim([-0.01 0.04]);

subplot(2,1,2), hold on
plot(times, squeeze(pts_intmu(2,:,:)));
plot(times, squeeze(pts_intmu(2,:,:)) + 2* sqrt(squeeze(pts_intP(2,2,:,:))), '-.');
plot(times, squeeze(pts_intmu(2,:,:)) - 2* sqrt(squeeze(pts_intP(2,2,:,:))), '-.');
plot(times, interp_state(2,:), 'b', 'LineWidth', 2);
xlabel('time'); ylabel('trend');
ylim([-0.2 0.3]);

wid = 6; hei = 4;
remove_whitespace_for_saveas;
print(fig, '-dpdf', 'example_filter_state.pdf');

%% filter-smoother state estimate
fig = figure(5); hold on;
pts = kiti_pts;
pts_intmu = cat(3, pts.intmu);
pts_intP = cat(4, pts.intP);

subplot(2,1,1), hold on
plot(times, squeeze(pts_intmu(1,:,:)));
plot(times, squeeze(pts_intmu(1,:,:)) + 2* sqrt(squeeze(pts_intP(1,1,:,:))), '-.');
plot(times, squeeze(pts_intmu(1,:,:)) - 2* sqrt(squeeze(pts_intP(1,1,:,:))), '-.');
plot(times, interp_state(1,:), 'b', 'LineWidth', 2);
xlabel('time'); ylabel('value');
ylim([-0.01 0.04]);

subplot(2,1,2), hold on
plot(times, squeeze(pts_intmu(2,:,:)));
plot(times, squeeze(pts_intmu(2,:,:)) + 2* sqrt(squeeze(pts_intP(2,2,:,:))), '-.');
plot(times, squeeze(pts_intmu(2,:,:)) - 2* sqrt(squeeze(pts_intP(2,2,:,:))), '-.');
plot(times, interp_state(2,:), 'b', 'LineWidth', 2);
xlabel('time'); ylabel('trend');
ylim([-0.2 0.3]);

wid = 6; hei = 4;
remove_whitespace_for_saveas;
print(fig, '-dpdf', 'example_filtersmoother_state.pdf');

%% smoother state estimate
fig = figure(6); hold on;
pts = smooth_pts;
pts_intmu = cat(3, pts.intmu);
pts_intP = cat(4, pts.intP);

subplot(2,1,1), hold on
plot(times, squeeze(pts_intmu(1,:,:)));
plot(times, squeeze(pts_intmu(1,:,:)) + 2* sqrt(squeeze(pts_intP(1,1,:,:))), '-.');
plot(times, squeeze(pts_intmu(1,:,:)) - 2* sqrt(squeeze(pts_intP(1,1,:,:))), '-.');
plot(times, interp_state(1,:), 'b', 'LineWidth', 2);
xlabel('time'); ylabel('value');
ylim([-0.01 0.04]);

subplot(2,1,2), hold on
plot(times, squeeze(pts_intmu(2,:,:)));
plot(times, squeeze(pts_intmu(2,:,:)) + 2* sqrt(squeeze(pts_intP(2,2,:,:))), '-.');
plot(times, squeeze(pts_intmu(2,:,:)) - 2* sqrt(squeeze(pts_intP(2,2,:,:))), '-.');
plot(times, interp_state(2,:), 'b', 'LineWidth', 2);
xlabel('time'); ylabel('trend');
ylim([-0.2 0.3]);

wid = 6; hei = 4;
remove_whitespace_for_saveas;
print(fig, '-dpdf', 'example_smoother_state.pdf');