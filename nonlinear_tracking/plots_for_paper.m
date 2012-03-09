clup

%% Problem

load('M1withoutVelresults.mat')
plot_2D_benchmark_for_paper( flags, params, times, true_intx, observs, [] )

wid = 4; hei = 3;
format_graph_for_pdf;
print(gcf, '-dpdf', ['benchmark_problem.pdf']);

%% M1 without

load('M1withoutVelresults.mat')
plot_2D_benchmark_for_paper( flags, params, times, true_intx, [], pts )

wid = 4; hei = 3;
format_graph_for_pdf;
print(gcf, '-dpdf', ['M1withoutVel.pdf']);

%% M1 with

load('M1withVelresults.mat')
plot_2D_benchmark_for_paper( flags, params, times, true_intx, [], pts )

wid = 4; hei = 3;
format_graph_for_pdf;
print(gcf, '-dpdf', ['M1withVel.pdf']);

%% M2 without

load('M2withoutVelresults.mat')
plot_2D_benchmark_for_paper( flags, params, times, true_intx, [], pts )

wid = 4; hei = 3;
format_graph_for_pdf;
print(gcf, '-dpdf', ['M2withoutVel.pdf']);

%% M2 with

load('M2withVelresults.mat')
plot_2D_benchmark_for_paper( flags, params, times, true_intx, [], pts )

wid = 4; hei = 3;
format_graph_for_pdf;
print(gcf, '-dpdf', ['M2withVel.pdf']);

%% 3D Problem

load('3DCartesian.mat')
plot_3D_benchmark_for_paper( flags, params, times, true_intx, observs, [] )

wid = 4; hei = 3;
format_graph_for_pdf;
print(gcf, '-dpdf', ['3Dbenchmark_problem.pdf']);

%% 3D Cartesian

load('3DCartesian.mat')
plot_3D_benchmark_for_paper( flags, params, times, true_intx, [], pts )

wid = 4; hei = 3;
format_graph_for_pdf;
print(gcf, '-dpdf', ['3DCartesian.pdf']);

%% 3D Intrinsic

load('3DIntrinsic.mat')
plot_3D_benchmark_for_paper( flags, params, times, true_intx, [], pts )

wid = 4; hei = 3;
format_graph_for_pdf;
print(gcf, '-dpdf', ['3DIntrinsic.pdf']);