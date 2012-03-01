close all

%% Data
fig = figure(1); clf;

hold on;
plot(times, observ, 'r', 'LineWidth', 1);
xlabel('time'), ylabel('x');
%ylim([-0.01, 0.05]);

wid = 4.5; hei = 2;
remove_whitespace_for_saveas;
print(fig, '-dpdf', 'fx_data.pdf');

%% filter kd estimate
fig = figure(2); clf;

subplot(2,1,1); hold on;
plot(filt_kd_est.times, filt_kd_est.jump_1_kd, 'b');
xlabel('time'); ylabel('x jump k.d. estimate');
ylim([0,2]);

subplot(2,1,2); hold on;
plot(filt_kd_est.times, filt_kd_est.jump_2_kd, 'b');
xlabel('time'); ylabel('x-dot jump k.d. estimate');
ylim([0,2]);

wid = 4.5; hei = 4;
remove_whitespace_for_saveas;
print(fig, '-dpdf', 'fx_filter_kdest.pdf');

%% smoother kd estimate
fig = figure(3); clf;

subplot(2,1,1); hold on;
plot(smooth_kd_est.times, smooth_kd_est.jump_1_kd, 'b');
xlabel('time'); ylabel('x jump k.d. estimate');
ylim([0,2]);

subplot(2,1,2); hold on;
plot(smooth_kd_est.times, smooth_kd_est.jump_2_kd, 'b');
xlabel('time'); ylabel('x-dot jump k.d. estimate');
ylim([0,2]);

wid = 4.5; hei = 4;
remove_whitespace_for_saveas;
print(fig, '-dpdf', 'fx_smoother_kdest.pdf');