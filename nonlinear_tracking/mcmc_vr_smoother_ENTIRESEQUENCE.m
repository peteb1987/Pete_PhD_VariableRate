function [ pts ] = mcmc_vr_smoother( flags, params, filt_part_sets, times, observ )
%MCMC_VR_SMOOTHER Smoother for variable rate models using MCMC method.

% Set some local variables
K = params.K; S = params.S; Np = params.Np;
ds = params.state_dim; dr = params.rnd_dim; do = params.obs_dim;

% Create a cell array for the smoothed particles
% pts = cell(S,1);
pts = initialise_particles(flags, params, S, observ);

% Loop through smoothing particles
for ii = 1:S
    
    % Initialise particle with a random final frame filtering particle
    old_pt = filt_part_sets{K}(unidrnd(Np));
    
    % Append a state just after the final time
    old_pt.Ns = old_pt.Ns+1;
    old_pt.tau = [old_pt.tau, times(end)+eps(times(end))];
    old_pt.x = [old_pt.x, old_pt.intx(:,end)];
    old_pt.w = [old_pt.w, zeros(dr,1)];
    old_pt.tau_prob = [old_pt.tau_prob; 0];
    old_pt.w_prob = [old_pt.w_prob; 0];
    
    % Loop backwards through time
    for k = K:-1:2
        
        t = times(k);
        
        % Propose a new past
        jj = unidrnd(Np);
        copy_pt = filt_part_sets{k}(jj);
        ppsl = log(unidpdf(jj, Np));
        rev_ppsl = ppsl;
        
        % Find the indexes for the start and end of the bridge region
        copy_start_tau_idx = find_nearest(copy_pt.tau, t, false);
        copy_stop_tau_idx = find_nearest(copy_pt.tau, t, true);
        copy_start_t_idx = find_nearest(times, copy_pt.tau(copy_start_tau_idx), true);
        copy_stop_t_idx = find_nearest(times, copy_pt.tau(copy_stop_tau_idx), false);
        old_start_tau_idx = find_nearest(old_pt.tau, t, false);
        old_stop_tau_idx = find_nearest(old_pt.tau, t, true);
        old_start_t_idx = find_nearest(times, old_pt.tau(old_start_tau_idx), true);
        old_stop_t_idx = find_nearest(times, old_pt.tau(old_stop_tau_idx), false);
        
        % Construct the new particle
        new_pt.x  = [copy_pt.x(:,1:copy_start_tau_idx) old_pt.x(:,old_stop_tau_idx:end)];
        new_pt.w  = [copy_pt.w(:,1:copy_start_tau_idx) old_pt.w(:,old_stop_tau_idx:end)];
        new_pt.tau  = [copy_pt.tau(:,1:copy_start_tau_idx) old_pt.tau(:,old_stop_tau_idx:end)];
        new_pt.Ns = size(new_pt.x,2);
        new_pt.intx = [copy_pt.intx(:,1:k) old_pt.intx(:,k+1:end)];
        new_pt.lhood = [copy_pt.lhood(1:k); old_pt.lhood(k+1:end)];
        new_pt.tau_prob  = [copy_pt.tau_prob(1:copy_start_tau_idx); old_pt.tau_prob(old_stop_tau_idx:end)];
        new_pt.w_prob  = [copy_pt.w_prob(1:copy_start_tau_idx); old_pt.w_prob(old_stop_tau_idx:end)];
        
        % Find the indexes for the new particle
        new_start_tau_idx = find_nearest(new_pt.tau, t, false);
        new_stop_tau_idx = find_nearest(new_pt.tau, t, true);
        new_start_t_idx = find_nearest(times, new_pt.tau(new_start_tau_idx), true);
        new_stop_t_idx = find_nearest(times, new_pt.tau(new_stop_tau_idx), false);
        
        % Calculate jump transition probabilities
        [~, new_jump_trans_prob] = sample_jump_time(flags, params, new_pt.tau(new_start_tau_idx), t, new_pt.tau(new_stop_tau_idx));
        [~, old_jump_trans_prob] = sample_jump_time(flags, params, old_pt.tau(old_start_tau_idx), t, old_pt.tau(old_stop_tau_idx));
        
        % Calculate accelerations and associated probabilities
        w_new = determine_acceleration(flags, params, new_pt.x(:,new_start_tau_idx), new_pt.x(:,new_stop_tau_idx), new_pt.tau(new_stop_tau_idx)-new_pt.tau(new_start_tau_idx));
        w_old = determine_acceleration(flags, params, old_pt.x(:,old_start_tau_idx), old_pt.x(:,old_stop_tau_idx), old_pt.tau(old_stop_tau_idx)-old_pt.tau(old_start_tau_idx));
        new_pt.w(:,new_start_tau_idx) = w_new;
        new_accel_prob = log(mvnpdf(w_new', zeros(1,dr), params.Q));
        old_accel_prob = log(mvnpdf(w_old', zeros(1,dr), params.Q));
        new_pt.tau_prob(new_start_tau_idx) = new_jump_trans_prob;
        new_pt.w_prob(new_start_tau_idx) = new_accel_prob;
        
        % Calculate bridging likelihoods
        [new_pt.intx(:,new_start_t_idx:new_stop_t_idx), new_pt.lhood(new_start_t_idx:new_stop_t_idx)] = ...
            interpolate_state(flags, params, new_pt.tau(new_start_tau_idx), new_pt.x(:,new_start_tau_idx), w_new, times(new_start_t_idx:new_stop_t_idx), observ(:,new_start_t_idx:new_stop_t_idx));
        [old_pt.intx(:,old_start_t_idx:old_stop_t_idx), old_pt.lhood(old_start_t_idx:old_stop_t_idx)] = ...
            interpolate_state(flags, params, old_pt.tau(old_start_tau_idx), old_pt.x(:,old_start_tau_idx), w_old, times(old_start_t_idx:old_stop_t_idx), observ(:,old_start_t_idx:old_stop_t_idx));
        new_bridge_lhood = sum(new_pt.lhood(new_start_t_idx:new_stop_t_idx));
        old_bridge_lhood = sum(old_pt.lhood(old_start_t_idx:old_stop_t_idx));
        
        % Evaluate past probabilities (stored)
        new_past_lhood = sum(new_pt.lhood(1:new_start_t_idx-1));
        old_past_lhood = sum(old_pt.lhood(1:old_start_t_idx-1));
        new_past_trans_prob = sum(new_pt.tau_prob(1:new_start_tau_idx)) + sum(new_pt.w_prob(1:new_start_tau_idx-1));
        old_past_trans_prob = sum(old_pt.tau_prob(1:old_start_tau_idx)) + sum(old_pt.w_prob(1:old_start_tau_idx-1));
        
        % Calculate MH acceptance probability
        acc_prob = (new_past_lhood+new_past_trans_prob+new_bridge_lhood+new_accel_prob+new_jump_trans_prob) ...
                  -(old_past_lhood+old_past_trans_prob+old_bridge_lhood+old_accel_prob+old_jump_trans_prob) ...
                  +(rev_ppsl-ppsl);
        
        % Test for acceptance
        if log(rand) < acc_prob
            
            % Keep proposed particle
            old_pt = new_pt;
            
        end
        
    end
    
    % Cut off the end jump
    old_pt.Ns = old_pt.Ns - 1;
    old_pt.tau(end) = [];
    old_pt.x(:,end) = [];
    old_pt.w(:,end) = [];
    old_pt.w_prob(end) = [];
    old_pt.tau_prob(end) = [];
    
    % Store the particle
    pts(ii) = old_pt;
    
    % Output
    fprintf('*** Completed MCMC iteration %d.\n', ii);
    
end

end

