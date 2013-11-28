% Base script for variable rate tracking

%% Preliminaries

if ~exist('test', 'var') || ~isfield(test,'flag_batch') || (~test.flag_batch)
    
    clup
    dbstop if error
    %     dbstop if warning
    
    % Set flag to non-batch
    test.flag_batch = false;
    
    %%% SETTINGS %%%
    
    % DEFINE RANDOM SEED
    rand_seed = 0;
    
    % Set display options
    display.text = true;
    display.plot_during = false;
    display.plot_after = true;
    
    % Tests to run
    
end

%% Setup
fprintf('   Random seed: %u.\n', rand_seed);

% Set random seed
rng(rand_seed);

% Set model parameters
[model] = tracking_setmodel(test);

% Set algorithm parameters
[algo] = tracking_setalgo(test);

% Generate data
[time, cp_time, cp_param, cp_state, state, observ] = tracking_generatedata(model);

%% Plot data

[ox, oy] = pol2cart(observ(1,:), observ(2,:));
figure, hold on
plot(state(1,:), state(2,:), 'k--');
plot(ox, oy, 'r.')
plot(cp_state(1,:), cp_state(2,:), 'g*')

%% Filtering

% Run the particle filter
[pf, pf_diagnostics] = vr_particle_filter(display, algo, model, time, observ);


% %% Analysis
% pos_rmse = cell(num_alg,1);
% vel_rmse = cell(num_alg,1);
% tnees = cell(num_alg,1);
% nus = cell(num_alg,1);
% for aa = 1:num_alg
%     [pos_rmse{aa}, vel_rmse{aa}, tnees{aa}, nus{aa}] = particle_smoother_analysis(model, state, ps{aa});
% end
% 
% %% Display
% if display.plot_after
%     close all
%     
%     colours = 'kbgrm';
%     
%     if model.ds == 4
%         
%         [ox, oy] = pol2cart(observ(1,:), observ(2,:));
%         figure, hold on
%         plot(state(1,:), state(2,:), 'k--');
%         plot(ox, oy, 'r-*')
%         
%         for aa = 1:num_alg
%             figure, hold on
%             plot(state(1,:), state(2,:), 'k--');
%             for ii = 1:length(ps{aa}), plot(ps{aa}(ii).traj_state(1,:), ps{aa}(ii).traj_state(2,:), colours(aa)); end
%         end
%         
%     elseif model.ds == 6
%         
%         [ox, oy, oz] = sph2cart(observ(1,:), observ(2,:), observ(3,:));
%         figure, hold on
%         plot3(state(1,:), state(2,:), state(3,:));
%         plot3(ox, oy, oz, 'r-*')
%         
%         for aa = 1:num_alg
%             figure, hold on
%             plot3(state(1,:), state(2,:), state(3,:));
%             for ii = 1:length(ps{aa}), plot3(ps{aa}(ii).traj_state(1,:), ps{aa}(ii).traj_state(2,:), ps{aa}(ii).traj_state(3,:), colours(aa)); end
%         end
%         
%     end
%     
%     figure, hold on
%     for aa = 1:num_alg
%         plot(pos_rmse{aa}, colours(aa));
%     end
%     
%     figure, hold on
%     for aa = 1:num_alg
%         plot(vel_rmse{aa}, colours(aa));
%     end
%     
%     figure, hold on
%     for aa = 1:num_alg
%         plot(tnees{aa}, colours(aa));
%     end
%     
%     figure, hold on
%     for aa = 1:num_alg
%         plot(nus{aa}, colours(aa));
%     end
%     
% end
