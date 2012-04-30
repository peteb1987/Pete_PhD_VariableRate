clup

%% Problem

load('2DSmoother.mat')
plot_2D_simulated_for_paper( flags, params, times, true_intx, observs, true_tau, true_x, [] )
xlabel('x'); ylabel('y');

wid = 4; hei = 3;
format_graph_for_pdf;
print(gcf, '-dpdf', ['simulated_problem.pdf']);

%% Filter

load('2DSmoother.mat')
plot_2D_simulated_for_paper( flags, params, times, true_intx, [], [], [], filt_pts )
xlabel('x'); ylabel('y');

wid = 4; hei = 3;
format_graph_for_pdf;
print(gcf, '-dpdf', ['2Dfilter.pdf']);

%% Smoother

load('2DSmoother.mat')
plot_2D_simulated_for_paper( flags, params, times, true_intx, [], [], [], smooth_pts )
xlabel('x'); ylabel('y');

wid = 4; hei = 3;
format_graph_for_pdf;
print(gcf, '-dpdf', ['2Dsmoother.pdf']);

%% Unique Particles

load('filtsmooth_comparison_results.mat')
figure, hold on;
plot(times, mean([results.kita_Nup], 1), 'b--')
plot(times, mean([results.smooth_Nup], 1), 'g')
xlabel('time'); ylabel('Number of unique particles');

wid = 3; hei = 2;
format_graph_for_pdf;
print(gcf, '-dpdf', ['tracking_unique_particles.pdf']);