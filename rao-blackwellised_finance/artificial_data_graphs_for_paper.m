close all

%% Data
fig = figure(1); clf;

subplot(2,1,1); hold on;
plot(times, interp_state(1,:), 'b', 'LineWidth', 2);
plot(times, observ, 'r', 'LineWidth', 1);
xlabel('time'), ylabel('$z$', 'interpreter', 'latex');
ylim([-0.01, 0.05]);

subplot(2,1,2); hold on;
plot(times, interp_state(2,:), 'b', 'LineWidth', 2);
xlabel('time'), ylabel('$\dot{z}$', 'interpreter', 'latex');
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
xlabel('time'); ylabel('$\dot{z}$', 'interpreter', 'latex');
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
xlabel('time'); ylabel('$\dot{z}$', 'interpreter', 'latex');
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
xlabel('time'); ylabel('$\dot{z}$', 'interpreter', 'latex');
ylim([-0.2 0.3]);

wid = 4.5; hei = 2;
remove_whitespace_for_saveas;
print(fig, '-dpdf', 'example_smoother_state.pdf');

%% Multimodality illustration

wid = 3; hei = 2;

pts = kiti_pts;
pts_intmu = cat(3, pts.intmu);
pts_intP = cat(4, pts.intP);

figure, hold on
plot(times, squeeze(pts_intmu(2,:,:)));
plot(times, interp_state(2,:), 'b', 'LineWidth', 2);
xlabel('time'); ylabel('$\dot{z}$', 'interpreter', 'latex');
xlim([0.55, 0.8]);
ylim([-0.05 0.1]);

pos = get(gca, 'Position');
set(gca, 'Position', pos + [0.07 0.07 0 0]);
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), wid, hei]);

export_pdf(gcf, 'finance_filter_multimodality.pdf');

pts = smooth_pts;
pts_intmu = cat(3, pts.intmu);
pts_intP = cat(4, pts.intP);

figure, hold on
plot(times, squeeze(pts_intmu(2,:,:)));
plot(times, interp_state(2,:), 'b', 'LineWidth', 2);
xlabel('time'); ylabel('$\dot{z}$', 'interpreter', 'latex');
xlim([0.55, 0.8]);
ylim([-0.05 0.1]);

pos = get(gca, 'Position');
set(gca, 'Position', pos + [0.07 0.07 0 0]);
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1), pos(2), wid, hei]);

export_pdf(gcf, 'finance_smoother_multimodality.pdf');

%%

figure, hold on,
VRPS_test = mean(cat(1,results.VRPS.unique_over_time), 1);
kita_test = mean(cat(1,results.kita.unique_over_time), 1);
end_test = 1000;
plot(times(:,1:end_test), VRPS_test(:,1:end_test), 'g', 'linewidth', 2)
plot(times(:,1:end_test), kita_test(:,1:end_test), '--b', 'linewidth', 2)
xlabel('time'); ylabel('Number of unique particles');
% orient landscape

% wid = 4; hei = 3;
% remove_whitespace_for_saveas;
% print(gcf, '-dpdf', ['finance_unique_particles.pdf']);
export_pdf(gcf, 'finance_unique_particles.pdf', 4, 4, 'inches');

%% MENEES

figure, hold on

load('results_100rpts_NF1000_filteronly.mat');
plot(times, mean(cat(1,results.filt.ENEES)), 'r:', 'linewidth', 2)
plot(times, mean(cat(1,results.kita.ENEES)), 'b--', 'linewidth', 2)

load('results_100rpts.mat')
plot(times, mean(cat(1,results.VRPS.ENEES)), 'g', 'linewidth', 2)

xlabel('time'); ylabel('Mean TNEES');
ylim([0.5 1.05]);
legend('Filter', 'Filter-smoother', 'Smoother', 'Location', 'SouthEast');

plot( [times(1) times(end)], [0.54 0.54], 'k:' );
ylim([0.5 1.1]);

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 6 3]);

export_pdf(gcf, 'finance_mtnees.pdf', 4, 2, 'inches');

%% RMSE

figure, hold on

load('results_100rpts_NF1000_filteronly.mat');
plot(times, mean(cell2mat(cat(2,arrayfun(@(x) {x.mean_rmse.value_over_time}, results.filt))')), 'r:', 'linewidth', 2)
plot(times, mean(cell2mat(cat(2,arrayfun(@(x) {x.mean_rmse.value_over_time}, results.kita))')), 'b--', 'linewidth', 2)

mean(mean(cell2mat(cat(2,arrayfun(@(x) {x.mean_rmse.value_over_time}, results.filt))')))
mean(mean(cell2mat(cat(2,arrayfun(@(x) {x.mean_rmse.value_over_time}, results.kita))')))

load('results_100rpts.mat')
plot(times, mean(cell2mat(cat(2,arrayfun(@(x) {x.mean_rmse.value_over_time}, results.VRPS))')), 'g', 'linewidth', 2)

mean(mean(cell2mat(cat(2,arrayfun(@(x) {x.mean_rmse.value_over_time}, results.VRPS))')))

xlabel('time'); ylabel('Mean value error');
legend('Filter', 'Filter-smoother', 'Smoother', 'Location', 'SouthEast');
ylim([0 0.6E-3]);

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 6 3]);

export_pdf(gcf, 'finance_value_rmse.pdf', 4, 2, 'inches');




figure, hold on

load('results_100rpts_NF1000_filteronly.mat');
plot(times, mean(cell2mat(cat(2,arrayfun(@(x) {x.mean_rmse.trend_over_time}, results.filt))')), 'r:', 'linewidth', 2)
plot(times, mean(cell2mat(cat(2,arrayfun(@(x) {x.mean_rmse.trend_over_time}, results.kita))')), 'b--', 'linewidth', 2)

mean(mean(cell2mat(cat(2,arrayfun(@(x) {x.mean_rmse.trend_over_time}, results.filt))')))
mean(mean(cell2mat(cat(2,arrayfun(@(x) {x.mean_rmse.trend_over_time}, results.kita))')))

load('results_100rpts.mat')
plot(times, mean(cell2mat(cat(2,arrayfun(@(x) {x.mean_rmse.trend_over_time}, results.VRPS))')), 'g', 'linewidth', 2)

mean(mean(cell2mat(cat(2,arrayfun(@(x) {x.mean_rmse.trend_over_time}, results.VRPS))')))

xlabel('time'); ylabel('Mean trend error');
legend('Filter', 'Filter-smoother', 'Smoother', 'Location', 'SouthEast');

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 6 3]);

export_pdf(gcf, 'finance_trend_rmse.pdf', 4, 2, 'inches');