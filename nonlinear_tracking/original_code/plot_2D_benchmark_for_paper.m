function plot_2D_benchmark_for_paper( flags, params, times, intx, observs, pts )
%PLOT_BENCHMARK_FOR_PAPER

figure, hold on,

% % Put an x on the origin
% plot(0, 0, 'xk', 'markersize', 10);

xlim([3E4, 7.5E4]); ylim([-1E4 3.5E4]);

% Plot track
plot(intx(1,:), intx(2,:), 'k:', 'linewidth', 1);
if ~isempty(observs)
    [x1, x2] = pol2cart(observs(1,:), observs(2,:));
    plot(x1, x2, 'rx');
end

% Plot particles
if ~isempty(pts)
    Np = length(pts);
    
    % Get interpolated states from particles
    pts_intx = cat(3,pts.intx);
    K = size(pts_intx, 2);
    x1 = squeeze(pts_intx(1,:,:));
    x2 = squeeze(pts_intx(2,:,:));
    
    plot(x1, x2, 'b');
    
    % Put inferred changpoints in
    used_tau = [];
    for ii = 1:Np
        for cpi = 1:pts(ii).Ns
            if ~any(pts(ii).tau(cpi)==used_tau)
                plot(pts(ii).x(1,cpi), pts(ii).x(2,cpi), 'g*')
                used_tau = [used_tau pts(ii).tau(cpi)];
            end
        end
    end
    
end

end

