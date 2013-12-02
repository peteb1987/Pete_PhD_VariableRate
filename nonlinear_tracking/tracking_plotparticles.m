function [ fig ] = tracking_plotparticles( algo, model, fig, K, time, pf )
%TRACKING_PLOTPARTICLES

if isempty(fig)
    fig = figure;
else
    figure(fig);
end
hold on;

% Plot trajectories
if ~isempty(pf)
    for ii = 1:length(pf)
        traj = zeros(model.ds,K);
        for kk = 1:K
            cp_ind = most_recent_changepoint(pf(ii).cp_time, time(kk));
            pre_cp_time = pf(ii).cp_time(cp_ind);
            pre_cp_param = pf(ii).cp_param(:,cp_ind);
            pre_cp_state = pf(ii).cp_state(:,cp_ind);
            traj(:,kk) = tracking_evaluatestate(model, pre_cp_time, pre_cp_param, pre_cp_state, time(kk));
        end
        plot(traj(1,:), traj(2,:))
        plot(pf(ii).cp_state(1,:), pf(ii).cp_state(2,:), 'g.');
    end
end

end

