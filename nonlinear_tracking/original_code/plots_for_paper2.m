clup

%% Problem

load('2DSmoother.mat')
plot_2D_simulated_for_paper( flags, params, times, true_intx, observs, true_tau, true_x, [] )
xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');

wid = 4; hei = 3;
format_graph_for_pdf;
print(gcf, '-dpdf', ['simulated_problem.pdf']);

%% Filter

load('2DFilter.mat')
plot_2D_simulated_for_paper( flags, params, times, true_intx, [], [], [], filt_pts )
xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');

wid = 4; hei = 3;
format_graph_for_pdf;
print(gcf, '-dpdf', ['2Dfilter.pdf']);

%% Smoother

load('2DSmoother.mat')
plot_2D_simulated_for_paper( flags, params, times, true_intx, [], [], [], smooth_pts )
xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');

wid = 4; hei = 3;
format_graph_for_pdf;
print(gcf, '-dpdf', ['2Dsmoother.pdf']);

%% Unique Particles

load('filtsmooth_comparison_results100.mat')
figure, hold on;
plot(times, mean([results.kita_Nup], 1), 'b--', 'linewidth', 2)
plot(times, mean([results.smooth_Nup], 1), 'g', 'linewidth', 2)
xlabel('time'); ylabel('Number of unique particles');

export_pdf(gcf, 'tracking_unique_particles.pdf', 4, 4, 'inches');

% wid = 3; hei = 2;
% format_graph_for_pdf;
% print(gcf, '-dpdf', ['tracking_unique_particles.pdf']);

%% MENEES

figure, hold on

load('filtsmooth_comparison_results_Nf500_Ns0_100rpts.mat')
results.filt_ENEES(:,1) = NaN;

plot(times, mean(results.filt_ENEES), 'r:', 'linewidth', 2)
plot(times, mean(results.kita_ENEES), 'b--', 'linewidth', 2)

load('filtsmooth_comparison_results_Nf50_Ns50_100rpts.mat')

plot(times, mean(results.smooth_ENEES), 'g', 'linewidth', 2)
xlabel('time'); ylabel('Mean TNEES');
legend('Filter', 'Filter-smoother', 'Smoother', 'Location', 'SouthEast');

plot( [times(1) times(end)], [0.73 0.73], 'k:' );
ylim([0.5 1.1]);

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 6 3]);

export_pdf(gcf, 'tracking_MTNEES.pdf', 4, 2, 'inches');

%% RMSE

figure, hold on

load('filtsmooth_comparison_results_Nf500_Ns0_100rpts.mat')

plot(times, mean(cat(1,results.filt_rmse.pos_over_time)), 'r:', 'linewidth', 2)
plot(times, mean(cat(1,results.kita_rmse.pos_over_time)), 'b--', 'linewidth', 2)

mean(mean(cat(1,results.filt_rmse.pos_over_time)))
mean(mean(cat(1,results.kita_rmse.pos_over_time)))

load('filtsmooth_comparison_results_Nf50_Ns50_100rpts.mat')

plot(times, mean(cat(1,results.smooth_rmse.pos_over_time)), 'g', 'linewidth', 2)

mean(mean(cat(1,results.smooth_rmse.pos_over_time)))

xlabel('time'); ylabel('Mean position error');
legend('Filter', 'Filter-smoother', 'Smoother', 'Location', 'SouthEast');

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 6 3]);

export_pdf(gcf, 'tracking_pos_RMSE.pdf', 4, 2, 'inches');

figure, hold on

load('filtsmooth_comparison_results_Nf500_Ns0_100rpts.mat')

plot(times, mean(cat(1,results.filt_rmse.vel_over_time)), 'r:', 'linewidth', 2)
plot(times, mean(cat(1,results.kita_rmse.vel_over_time)), 'b--', 'linewidth', 2)

load('filtsmooth_comparison_results_Nf50_Ns50_100rpts.mat')

plot(times, mean(cat(1,results.smooth_rmse.vel_over_time)), 'g', 'linewidth', 2)
xlabel('time'); ylabel('Mean velocity error');
legend('Filter', 'Filter-smoother', 'Smoother', 'Location', 'SouthEast');

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 6 3]);

export_pdf(gcf, 'tracking_vel_RMSE.pdf', 4, 2, 'inches');

%% Varying number of smoothing particles

num_pts = [5 15 50 150];
pos_rmse = zeros(1,4);
vel_rmse = zeros(1,4);
menees = zeros(1,4);

% 5 particles
load('filtsmooth_comparison_results_Nf50_Ns5_100rpts');
pos_rmse(1) = mean(mean(cat(1,results.smooth_rmse.pos_over_time)));
vel_rmse(1) = mean(mean(cat(1,results.smooth_rmse.vel_over_time)));
menees(1) = mean(mean(cat(1,results.smooth_ENEES)));

% 15 particles
load('filtsmooth_comparison_results_Nf50_Ns15_100rpts');
pos_rmse(2) = mean(mean(cat(1,results.smooth_rmse.pos_over_time)));
vel_rmse(2) = mean(mean(cat(1,results.smooth_rmse.vel_over_time)));
menees(2) = mean(mean(cat(1,results.smooth_ENEES)));

% 50 particles
load('filtsmooth_comparison_results_Nf50_Ns50_100rpts');
pos_rmse(3) = mean(mean(cat(1,results.smooth_rmse.pos_over_time)));
vel_rmse(3) = mean(mean(cat(1,results.smooth_rmse.vel_over_time)));
menees(3) = mean(mean(cat(1,results.smooth_ENEES)));

% 150 particles
load('filtsmooth_comparison_results_Nf50_Ns150_100rpts');
pos_rmse(4) = mean(mean(cat(1,results.smooth_rmse.pos_over_time)));
vel_rmse(4) = mean(mean(cat(1,results.smooth_rmse.vel_over_time)));
menees(4) = mean(mean(cat(1,results.smooth_ENEES)));

figure, hold on
plot(num_pts, menees, 'linewidth', 2)
xlabel('Number of smoothing particles');
ylabel('MTNEES');

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 3 3]);

export_pdf(gcf, 'MTNEES_vs_numpts.pdf', 2, 2, 'inches');

figure, hold on
plot(num_pts, pos_rmse, 'linewidth', 2)
xlabel('Number of smoothing particles');
ylabel('Mean Absolute Position Error');

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 3 3]);

export_pdf(gcf, 'pos_RMSE_vs_numpts.pdf', 2, 2, 'inches');

figure, hold on
plot(num_pts, vel_rmse, 'linewidth', 2)
xlabel('Number of smoothing particles');
ylabel('Mean Absolute Velocity Error');

set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 3 3]);

export_pdf(gcf, 'vel_RMSE_vs_numpts.pdf', 2, 2, 'inches');