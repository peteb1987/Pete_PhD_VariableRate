function plot_results( flags, params, part_sets, fig )
%PLOT_RESULTS Plot RBVRPF results

if nargin > 3
    figure(fig); clf;
else
    fig = figure;
end

subplot(3,1,1), hold on
plot(times, part_sets{params.K}.pts_intmu(:,:,1)');
ylabel('x');
plot(times, interp_state(1,:), 'b', times, observ, 'r', 'LineWidth', 2);

subplot(3,1,2), hold on
plot(times, part_sets{params.K}.pts_intmu(:,:,2)');
ylabel('x-dot');
if flags.gen_data
    plot(times, interp_state(2,:), 'b', 'LineWidth', 2);
end

subplot(3,1,3), hold on
plot(filt_kd_times, filt_x_jump_kdest, 'b');
plot(filt_kd_times, filt_xdot_jump_kdest, 'g');
ylabel('jump kd estimate');
legend('x', 'x-dot');
if flags.gen_data
    for tt=1:length(tau)
        plot([tau(tt),tau(tt)], [0,1]','r');
    end
end

drawnow;

end

