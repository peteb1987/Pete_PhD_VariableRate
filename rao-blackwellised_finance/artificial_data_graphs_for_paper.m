close all

%% Data
fig = figure(1); clf;

subplot(2,1,1); hold on;
plot(times, interp_state(1,:), 'b', 'LineWidth', 2);
plot(times, observ, 'r', 'LineWidth', 1);
xlabel('time'), ylabel('x');
ylim([-0.01, 0.05]);

subplot(2,1,2); hold on;
plot(times, interp_state(2,:), 'b', 'LineWidth', 2);
xlabel('time'), ylabel('x-dot');
ylim([-0.2, 0.3]);

for tt=1:length(tau)
    if type(tt) == 1
        subplot(2,1,1);
        plot([tau(tt),tau(tt)], [-0.01,0.05]',':k');
    elseif type(tt) == 2
        subplot(2,1,2);
        plot([tau(tt),tau(tt)], [-0.2,0.3]',':k');
    end
end

wid = 4.5; hei = 4;
remove_whitespace_for_saveas;
print(fig, '-dpdf', 'example_data.pdf');

%% filter kd estimate
fig = figure(2); clf;

subplot(2,1,1); hold on;
plot(filt_kd_est.times, filt_kd_est.jump_1_kd, 'b');
xlabel('time'); ylabel('x jump k.d. estimate');

subplot(2,1,2); hold on;
plot(filt_kd_est.times, filt_kd_est.jump_2_kd, 'b');
xlabel('time'); ylabel('x-dot jump k.d. estimate');

for tt=1:length(tau)
    if type(tt) == 1
        subplot(2,1,1);
        plot([tau(tt),tau(tt)], [0,2]',':k');
    elseif type(tt) == 2
        subplot(2,1,2);
        plot([tau(tt),tau(tt)], [0,2]',':k');
    end
end

wid = 4.5; hei = 4;
remove_whitespace_for_saveas;
print(fig, '-dpdf', 'example_filter_kdest.pdf');

%% smoother kd estimate
fig = figure(3); clf;

subplot(2,1,1); hold on;
plot(smooth_kd_est.times, smooth_kd_est.jump_1_kd, 'b');
xlabel('time'); ylabel('x jump k.d. estimate');

subplot(2,1,2); hold on;
plot(smooth_kd_est.times, smooth_kd_est.jump_2_kd, 'b');
xlabel('time'); ylabel('x-dot jump k.d. estimate');

for tt=1:length(tau)
    if type(tt) == 1
        subplot(2,1,1);
        plot([tau(tt),tau(tt)], [0,2]',':k');
    elseif type(tt) == 2
        subplot(2,1,2);
        plot([tau(tt),tau(tt)], [0,2]',':k');
    end
end

wid = 4.5; hei = 4;
remove_whitespace_for_saveas;
print(fig, '-dpdf', 'example_smoother_kdest.pdf');


%% filter state estimate
fig = figure(4); clf;

pts = filt_part_sets{end};
pts_intmu = cat(3, pts.intmu);
pts_intP = cat(4, pts.intP);

% subplot(2,1,1), hold on
% plot(times, squeeze(pts_intmu(1,:,:)));
% plot(times, squeeze(pts_intmu(1,:,:)) + 2* sqrt(squeeze(pts_intP(1,1,:,:))), '-.');
% plot(times, squeeze(pts_intmu(1,:,:)) - 2* sqrt(squeeze(pts_intP(1,1,:,:))), '-.');
% plot(times, interp_state(1,:), 'b', 'LineWidth', 2);
% xlabel('time'); ylabel('value');
% ylim([-0.01 0.05]);
% 
% subplot(2,1,2), hold on
hold on;
plot(times, squeeze(pts_intmu(2,:,:)));
plot(times, squeeze(pts_intmu(2,:,:)) + 2* sqrt(squeeze(pts_intP(2,2,:,:))), '-.');
plot(times, squeeze(pts_intmu(2,:,:)) - 2* sqrt(squeeze(pts_intP(2,2,:,:))), '-.');
plot(times, interp_state(2,:), 'b', 'LineWidth', 2);
xlabel('time'); ylabel('trend');
ylim([-0.2 0.3]);

wid = 4.5; hei = 2;
remove_whitespace_for_saveas;
print(fig, '-dpdf', 'example_filter_state.pdf');

%% filter-smoother state estimate
fig = figure(5); clf;

pts = kiti_pts;
pts_intmu = cat(3, pts.intmu);
pts_intP = cat(4, pts.intP);

% subplot(2,1,1), hold on
% plot(times, squeeze(pts_intmu(1,:,:)));
% plot(times, squeeze(pts_intmu(1,:,:)) + 2* sqrt(squeeze(pts_intP(1,1,:,:))), '-.');
% plot(times, squeeze(pts_intmu(1,:,:)) - 2* sqrt(squeeze(pts_intP(1,1,:,:))), '-.');
% plot(times, interp_state(1,:), 'b', 'LineWidth', 2);
% xlabel('time'); ylabel('value');
% ylim([-0.01 0.05]);
% 
% subplot(2,1,2), hold on
hold on
plot(times, squeeze(pts_intmu(2,:,:)));
plot(times, squeeze(pts_intmu(2,:,:)) + 2* sqrt(squeeze(pts_intP(2,2,:,:))), '-.');
plot(times, squeeze(pts_intmu(2,:,:)) - 2* sqrt(squeeze(pts_intP(2,2,:,:))), '-.');
plot(times, interp_state(2,:), 'b', 'LineWidth', 2);
xlabel('time'); ylabel('trend');
ylim([-0.2 0.3]);

wid = 4.5; hei = 2;
remove_whitespace_for_saveas;
print(fig, '-dpdf', 'example_filtersmoother_state.pdf');

%% smoother state estimate
fig = figure(6); clf;

pts = smooth_pts;
pts_intmu = cat(3, pts.intmu);
pts_intP = cat(4, pts.intP);

% subplot(2,1,1), hold on
% plot(times, squeeze(pts_intmu(1,:,:)));
% plot(times, squeeze(pts_intmu(1,:,:)) + 2* sqrt(squeeze(pts_intP(1,1,:,:))), '-.');
% plot(times, squeeze(pts_intmu(1,:,:)) - 2* sqrt(squeeze(pts_intP(1,1,:,:))), '-.');
% plot(times, interp_state(1,:), 'b', 'LineWidth', 2);
% xlabel('time'); ylabel('value');
% ylim([-0.01 0.05]);
% 
% subplot(2,1,2), hold on
hold on
plot(times, squeeze(pts_intmu(2,:,:)));
plot(times, squeeze(pts_intmu(2,:,:)) + 2* sqrt(squeeze(pts_intP(2,2,:,:))), '-.');
plot(times, squeeze(pts_intmu(2,:,:)) - 2* sqrt(squeeze(pts_intP(2,2,:,:))), '-.');
plot(times, interp_state(2,:), 'b', 'LineWidth', 2);
xlabel('time'); ylabel('trend');
ylim([-0.2 0.3]);

wid = 4.5; hei = 2;
remove_whitespace_for_saveas;
print(fig, '-dpdf', 'example_smoother_state.pdf');

%%

figure, plot(times, mean(cat(1,results.VRPS.unique_over_time), 1))
close all
figure, hold on,
plot(times, mean(cat(1,results.VRPS.unique_over_time), 1), 'g')
plot(times, mean(cat(1,results.kita.unique_over_time), 1), '--b')

wid = 4; hei = 3;
remove_whitespace_for_saveas;
print(gcf, '-dpdf', ['finance_unique_particles.pdf']);