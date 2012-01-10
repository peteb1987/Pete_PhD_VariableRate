function plot_results( flags, params, times, observ, true_tau, true_intx, pts, kd_est, fig )
%PLOT_RESULTS Plot RBVRPF results

if nargin > 4
    figure(fig); clf;
else
    fig = figure;
end

pts_intmu = cat(3, pts.intmu);

subplot(3,1,1), hold on
plot(times, squeeze(pts_intmu(1,:,:)));
ylabel('x');
plot(times, observ, 'r', 'LineWidth', 2); 
if ~isempty(true_intx)
    plot(times, true_intx(1,:), 'b', 'LineWidth', 2);
end

subplot(3,1,2), hold on
plot(times, squeeze(pts_intmu(2,:,:)));
ylabel('x-dot');
if ~isempty(true_intx)
    plot(times, true_intx(2,:), 'b', 'LineWidth', 2);
end

subplot(3,1,3), hold on
plot(kd_est.times, kd_est.jump_1_kd, 'b');
plot(kd_est.times, kd_est.jump_2_kd, 'g');
ylabel('jump kd estimate');
legend('x', 'x-dot');
if flags.gen_data
    for tt=1:length(true_tau)
        plot([true_tau(tt),true_tau(tt)], [0,1]','r');
    end
end

drawnow;

end

