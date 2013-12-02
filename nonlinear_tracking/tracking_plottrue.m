function [ fig ] = tracking_plottrue( model, time, state, observ, cp_time, cp_param, cp_state )
%TRACKING_PLOTTRUE

fig = figure;
hold on

% Transform the observations
[ox, oy] = pol2cart(observ(1,:), observ(2,:));

% Plot
plot(state(1,:), state(2,:), 'k--');
plot(ox, oy, 'r.')
plot(cp_state(1,:), cp_state(2,:), 'g*')

end

