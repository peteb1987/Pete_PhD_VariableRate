function [ mNs, mospa, rmse ] = performance_measures( pts, true_tau, true_intx )
%PERFORMANCE_MEASURES Calculates various performance measures for RBVRPF/S
%output

% Count particles
Np = length(pts.pts_Ns);

% Calculate mean OSPA between particle jump sequences and real sequence
ospas = zeros(Np,1);
for ii = 1:Np
    ospas(ii) = OSPA(true_tau, pts.pts_tau(ii, 1:pts.pts_Ns(ii)), 1, 0.01);
end
mospa = mean(ospas);

% Calculate mean number of jumps
mNs = mean (pts.pts_Ns);

% Calculate price RMSE
error = abs(bsxfun(@minus, true_intx(1,:), squeeze(pts.pts_intmu(:,1,:))));
rmse_over_time = sqrt(mean(error.^2, 1));
rmse = sqrt(mean(error(:).^2));

figure, plot(rmse_over_time)

end

