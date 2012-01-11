function [ mNs, mospa, rmse, rmse_over_time ] = performance_measures( pts, times, true_tau, true_intx )
%PERFORMANCE_MEASURES Calculates various performance measures for RBVRPF/S
%output

mospa = [];
rmse = [];
rmse_over_time = [];

% Count particles
Np = length(pts);

% Calculate mean number of jumps
mNs = mean(cat(1,pts.Ns));

% Calculate mean OSPA between particle jump sequences and real sequence
if ~isempty(true_tau)
    ospas = zeros(Np,1);
    for ii = 1:Np
        ospas(ii) = OSPA(true_tau, pts(ii).tau(1:pts(ii).Ns), 1, 0.01);
    end
    mospa = mean(ospas);
end

% Calculate price RMSE
if ~isempty(true_intx)
    intmu = cat(3, pts.intmu);
    error = abs(bsxfun(@minus, true_intx(1,:), squeeze(intmu(1,:,:))'));
    rmse_over_time = sqrt(mean(error.^2, 1));
    rmse = sqrt(mean(error(:).^2));
end

% figure, plot(times, rmse_over_time)

end

