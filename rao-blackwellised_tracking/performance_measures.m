function [ mospa ] = performance_measures( pts, true_tau )
%PERFORMANCE_MEASURES Calculates various performance measures for RBVRPF/S
%output

Np = length(pts.pts_Ns);

ospas = zeros(Np,1);

for ii = 1:Np
    ospas(ii) = OSPA(true_tau, pts.pts_tau(ii, 1:pts.pts_Ns(ii)), 1, 0.01);
end

mospa = mean(ospas);

end

