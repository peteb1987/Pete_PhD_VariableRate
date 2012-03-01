function [ mNs, MOSPA_ospa, mean_rmse, MAP_rmse ] = performance_measures( params, pts, times, observs, true_tau, true_intx )
%PERFORMANCE_MEASURES Calculates various performance measures for RBVRPF/S
%output

MOSPA_ospa = [];
mean_rmse = [];
MAP_rmse = [];

% Count particles
Np = length(pts);

% Calculate mean number of jumps
mNs = mean(cat(1,pts.Ns));

% Calculate MMSE estimate RMSE
if ~isempty(true_intx)
    intmu = mean(cat(3, pts.intmu),3);
    value_error = abs(true_intx(1,:) - intmu(1,:));
    mean_rmse.value_over_time = value_error;
    mean_rmse.value = sqrt(mean(value_error(:).^2));
    trend_error = abs(true_intx(2,:) - intmu(2,:));
    mean_rmse.trend_over_time = trend_error;
    mean_rmse.trend = sqrt(mean(trend_error(:).^2));
end

MAP_pt = pick_max_particle(params, pts, times, observs);

% Find minimum OSPA estimate
MOSPA_pt = pick_MOSPA_particle(params, pts);

% Calculate MAP estimate RMSE
if ~isempty(true_intx)
    MAP_value_error = abs(bsxfun(@minus, true_intx(1,:), MAP_pt.intmu(1,:)));
    MAP_rmse.value_over_time = sqrt(mean(MAP_value_error.^2, 1));
    MAP_rmse.value = sqrt(mean(MAP_value_error(:).^2));
    MAP_trend_error = abs(bsxfun(@minus, true_intx(2,:), MAP_pt.intmu(2,:)));
    MAP_rmse.trend_over_time = sqrt(mean(MAP_trend_error.^2, 1));
    MAP_rmse.trend = sqrt(mean(MAP_trend_error(:).^2));
end

% Calculate MAP OSPA between particle jump sequences and real sequence
if ~isempty(true_tau)
    MAP_ospa = OSPA(true_tau, MAP_pt.tau, 1, 0.01)
    MOSPA_ospa = OSPA(true_tau, MOSPA_pt.tau, 1, 0.01)
end

end

