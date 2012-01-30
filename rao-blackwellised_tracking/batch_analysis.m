num_test = 10;

comb_filt.mospa = zeros(num_test,1);
comb_filt.unique_sequences = zeros(num_test,1);
comb_filt.unique_times = zeros(num_test,1);
comb_filt.MMSE_value_RMSE = zeros(num_test,1);
comb_filt.MMSE_trend_RMSE = zeros(num_test,1);
comb_filt.MAP_value_RMSE = zeros(num_test,1);
comb_filt.MAP_trend_RMSE = zeros(num_test,1);

comb_kiti.mospa = zeros(num_test,1);
comb_kiti.unique_sequences = zeros(num_test,1);
comb_kiti.unique_times = zeros(num_test,1);
comb_kiti.MMSE_value_RMSE = zeros(num_test,1);
comb_kiti.MMSE_trend_RMSE = zeros(num_test,1);
comb_kiti.MAP_value_RMSE = zeros(num_test,1);
comb_kiti.MAP_trend_RMSE = zeros(num_test,1);

comb_VRPS.mospa = zeros(num_test,1);
comb_VRPS.unique_sequences = zeros(num_test,1);
comb_VRPS.unique_times = zeros(num_test,1);
comb_VRPS.MMSE_value_RMSE = zeros(num_test,1);
comb_VRPS.MMSE_trend_RMSE = zeros(num_test,1);
comb_VRPS.MAP_value_RMSE = zeros(num_test,1);
comb_VRPS.MAP_trend_RMSE = zeros(num_test,1);

for test = 1:num_test
    
    filename  = ['RBVRPS_results_' num2str(test) '.mat'];
    load(filename);
    
    comb_filt.mospa(test) = filt.mospa;
    comb_filt.unique_sequences(test) = filt.unique_sequences;
    comb_filt.unique_times(test) = filt.unique_times;
    comb_filt.MMSE_value_RMSE(test) = filt.mean_rmse.value;
    comb_filt.MMSE_trend_RMSE(test) = filt.mean_rmse.trend;
    comb_filt.MAP_value_RMSE(test) = filt.MAP_rmse.value;
    comb_filt.MAP_trend_RMSE(test) = filt.MAP_rmse.trend;
    
    comb_kiti.mospa(test) = filt.mospa;
    comb_kiti.unique_sequences(test) = filt.unique_sequences;
    comb_kiti.unique_times(test) = filt.unique_times;
    comb_kiti.MMSE_value_RMSE(test) = filt.mean_rmse.value;
    comb_kiti.MMSE_trend_RMSE(test) = filt.mean_rmse.trend;
    comb_kiti.MAP_value_RMSE(test) = filt.MAP_rmse.value;
    comb_kiti.MAP_trend_RMSE(test) = filt.MAP_rmse.trend;
    
    comb_VRPS.mospa(test) = filt.mospa;
    comb_VRPS.unique_sequences(test) = filt.unique_sequences;
    comb_VRPS.unique_times(test) = filt.unique_times;
    comb_VRPS.MMSE_value_RMSE(test) = filt.mean_rmse.value;
    comb_VRPS.MMSE_trend_RMSE(test) = filt.mean_rmse.trend;
    comb_VRPS.MAP_value_RMSE(test) = filt.MAP_rmse.value;
    comb_VRPS.MAP_trend_RMSE(test) = filt.MAP_rmse.trend;
end

fprintf(1, 'MAP est. OSPA, No. unique sequences, No. unique jumps, MMSE est. value RMSE, MMSE est. trend RMSE, MAP est. value RMSE, MAP est. trend RMSE   \n');
fprintf(1, 'Filter:          %f  %f  %f  %f  %f  %f  %f.\n', mean(comb_filt.mospa), mean(comb_filt.unique_sequences), mean(comb_filt.unique_times), mean(comb_filt.MMSE_value_RMSE), mean(comb_filt.MMSE_trend_RMSE), mean(comb_filt.MAP_value_RMSE), mean(comb_filt.MAP_trend_RMSE));
fprintf(1, 'Filter-smoother: %f  %f  %f  %f  %f  %f  %f.\n', mean(comb_kiti.mospa), mean(comb_kiti.unique_sequences), mean(comb_kiti.unique_times), mean(comb_kiti.MMSE_value_RMSE), mean(comb_kiti.MMSE_trend_RMSE), mean(comb_kiti.MAP_value_RMSE), mean(comb_kiti.MAP_trend_RMSE));
fprintf(1, 'SMoother:        %f  %f  %f  %f  %f  %f  %f.\n', mean(comb_VRPS.mospa), mean(comb_VRPS.unique_sequences), mean(comb_VRPS.unique_times), mean(comb_VRPS.MMSE_value_RMSE), mean(comb_VRPS.MMSE_trend_RMSE), mean(comb_VRPS.MAP_value_RMSE), mean(comb_VRPS.MAP_trend_RMSE));
