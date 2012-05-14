clup

%% Problem

load('final_M1withoutVelresults.mat')
plot_2D_benchmark_for_paper( flags, params, times, true_intx, observs, [] )

export_pdf(gcf, 'benchmark_problem.pdf', 4, 3, 'inches');

%% M1 without

load('final_M1withoutVelresults.mat')
plot_2D_benchmark_for_paper( flags, params, times, true_intx, [], pts )

export_pdf(gcf, 'M1withoutVel.pdf', 4, 3, 'inches');

%% M1 with

load('final_M1withVelresults.mat')
plot_2D_benchmark_for_paper( flags, params, times, true_intx, [], pts )

export_pdf(gcf, 'M1withVel.pdf', 4, 3, 'inches');

%% M2 without

load('final_M2withoutVelresults.mat')
plot_2D_benchmark_for_paper( flags, params, times, true_intx, [], pts )

export_pdf(gcf, 'M2withoutVel.pdf', 4, 3, 'inches');

%% M2 with

load('final_M2withVelresults.mat')
plot_2D_benchmark_for_paper( flags, params, times, true_intx, [], pts )

export_pdf(gcf, 'M2withVel.pdf', 4, 3, 'inches');

%% 3D Problem

load('3DCartesian_moreRM.mat')
plot_3D_benchmark_for_paper( flags, params, times, true_intx, observs, [] )

export_pdf(gcf, '3Dbenchmark_problem.pdf', 4, 3, 'inches');

%% 3D Cartesian

load('3DCartesian_moreRM.mat')
plot_3D_benchmark_for_paper( flags, params, times, true_intx, [], pts )

export_pdf(gcf, '3DCartesian.pdf', 4, 3, 'inches');

%% 3D Intrinsic

load('3DIntrinsic_moreRM.mat')
plot_3D_benchmark_for_paper( flags, params, times, true_intx, [], pts )

export_pdf(gcf, '3DIntrinsic.pdf', 4, 3, 'inches');
