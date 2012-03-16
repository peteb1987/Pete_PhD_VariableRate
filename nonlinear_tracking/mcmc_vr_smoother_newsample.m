function [ pts ] = mcmc_vr_smoother_newsample( flags, params, filt_part_sets, filt_weight_sets, times, observ )
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
    ind = randsample(Np, 1, true, exp(filt_weight_sets{end}));
%     ind = unidrnd(Np);
    old_pt = filt_part_sets{K}(ind);
    
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
        
        for m = 1:params.M
            
            % Propose a new past
            %             jj = unidrnd(Np);
            jj = randsample(Np, 1, true, exp(filt_weight_sets{k}));
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
            
            if (new_start_tau_idx>1)&&(old_start_tau_idx>1)
                
                % Propose a change in the previous changepoint time
                tau_replace = inf;
                while (tau_replace>times(k))||(tau_replace<new_pt.tau(new_start_tau_idx-1))
                    tau_replace = normrnd(new_pt.tau(new_start_tau_idx), params.ppsl_move_time_sd);
                end
                %             ppsl_tau_prob = log(normpdf(tau_replace, last_tau, params.ppsl_move_time_sd));
                new_pt.tau(new_start_tau_idx) = tau_replace;
                
                % Calculate previous changepoint probabilities
                [~, new_past_cptime_prob] = sample_jump_time(flags, params, new_pt.tau(new_start_tau_idx-1), [], new_pt.tau(new_start_tau_idx));
                [~, old_past_cptime_prob] = sample_jump_time(flags, params, old_pt.tau(old_start_tau_idx-1), [], old_pt.tau(old_start_tau_idx));
                new_pt.tau_prob(new_start_tau_idx) = new_past_cptime_prob;
                
                % Propose a new ante-previous parameter
                [w_replace, ~, ~] = acceleration_proposal(flags, params, new_pt.tau(new_start_tau_idx-1), new_pt.x(:,new_start_tau_idx-1), new_pt.w(:,new_start_tau_idx-1), times(1:new_start_t_idx-1), observ(:,1:new_start_t_idx-1));
                new_pt.w(:,new_start_tau_idx-1) = w_replace;
                new_pt.x(:,new_start_tau_idx) = next_state(flags, params,  new_pt.x(:,new_start_tau_idx-1), new_pt.w(:,new_start_tau_idx-1), new_pt.tau(new_start_tau_idx)-new_pt.tau(new_start_tau_idx-1));
                
                % Calculate ante-previous parameter probabilities
                new_past_param_prob = log(mvnpdf(new_pt.w(:,new_start_tau_idx)', zeros(1,dr), params.Q));
                old_past_param_prob = log(mvnpdf(old_pt.w(:,old_start_tau_idx)', zeros(1,dr), params.Q));
                new_pt.w_prob(new_start_tau_idx) = new_past_param_prob;
                
                % Update time indexes
                new_start_t_idx = find_nearest(times, new_pt.tau(new_start_tau_idx), true);
                new_past_start_t_idx = find_nearest(times, new_pt.tau(new_start_tau_idx-1), true);
                old_past_start_t_idx = find_nearest(times, old_pt.tau(old_start_tau_idx-1), true);
                
                % Likelihood bits
                new_past_replaced_lhood = sum(new_pt.lhood(new_past_start_t_idx:new_start_t_idx));
                old_past_replaced_lhood = sum(old_pt.lhood(old_past_start_t_idx:old_start_t_idx));
                
                [new_pt.intx(:,new_past_start_t_idx:new_start_t_idx-1), new_pt.lhood(new_past_start_t_idx:new_start_t_idx-1)] = ...
                    interpolate_state(flags, params, new_pt.tau(new_start_tau_idx-1), new_pt.x(:,new_start_tau_idx-1), w_replace, times(new_past_start_t_idx:new_start_t_idx-1), observ(:,new_past_start_t_idx:new_start_t_idx-1));
                w_old = old_pt.w(:,old_start_tau_idx-1);
                [old_pt.intx(:,old_past_start_t_idx:old_start_t_idx-1), old_pt.lhood(old_past_start_t_idx:old_start_t_idx-1)] = ...
                    interpolate_state(flags, params, old_pt.tau(old_start_tau_idx-1), old_pt.x(:,old_start_tau_idx-1), w_old, times(old_past_start_t_idx:old_start_t_idx-1), observ(:,old_past_start_t_idx:old_start_t_idx-1));
                new_past_bridge_lhood = sum(new_pt.lhood(new_past_start_t_idx:new_start_t_idx-1));
                old_past_bridge_lhood = sum(old_pt.lhood(old_past_start_t_idx:old_start_t_idx-1));
                
            else
                
                new_past_cptime_prob = 0;
                old_past_cptime_prob = 0;
                new_past_replaced_lhood = 0;
                old_past_replaced_lhood = 0;
                new_past_bridge_lhood = 0;
                old_past_bridge_lhood = 0;
            
            end
            
            % Evaluate replaced filtering probabilities (stored)
            new_replaced_lhood = sum(new_pt.lhood(new_start_t_idx:k));
            old_replaced_lhood = sum(old_pt.lhood(old_start_t_idx:k));
            new_replaced_param_prob = new_pt.w_prob(new_start_tau_idx);
            old_replaced_param_prob = old_pt.w_prob(old_start_tau_idx);
            
            % Calculate changepoint transition probabilities
            [~, new_cptime_prob] = sample_jump_time(flags, params, new_pt.tau(new_start_tau_idx), t, new_pt.tau(new_stop_tau_idx));
            [~, old_cptime_prob] = sample_jump_time(flags, params, old_pt.tau(old_start_tau_idx), t, old_pt.tau(old_stop_tau_idx));
            
            % Calculate accelerations and associated probabilities
            w_new = determine_acceleration(flags, params, new_pt.x(:,new_start_tau_idx), new_pt.x(:,new_stop_tau_idx), new_pt.tau(new_stop_tau_idx)-new_pt.tau(new_start_tau_idx));
            w_old = determine_acceleration(flags, params, old_pt.x(:,old_start_tau_idx), old_pt.x(:,old_stop_tau_idx), old_pt.tau(old_stop_tau_idx)-old_pt.tau(old_start_tau_idx));
            new_param_prob = log(mvnpdf(w_new', zeros(1,dr), params.Q));
            old_param_prob = log(mvnpdf(w_old', zeros(1,dr), params.Q));
            
            new_pt.w(:,new_start_tau_idx) = w_new;
            new_pt.tau_prob(new_start_tau_idx) = new_cptime_prob;
            new_pt.w_prob(new_start_tau_idx) = new_param_prob;
            
            % Calculate bridging likelihoods
            [new_pt.intx(:,new_start_t_idx:new_stop_t_idx), new_pt.lhood(new_start_t_idx:new_stop_t_idx)] = ...
                interpolate_state(flags, params, new_pt.tau(new_start_tau_idx), new_pt.x(:,new_start_tau_idx), w_new, times(new_start_t_idx:new_stop_t_idx), observ(:,new_start_t_idx:new_stop_t_idx));
            [old_pt.intx(:,old_start_t_idx:old_stop_t_idx), old_pt.lhood(old_start_t_idx:old_stop_t_idx)] = ...
                interpolate_state(flags, params, old_pt.tau(old_start_tau_idx), old_pt.x(:,old_start_tau_idx), w_old, times(old_start_t_idx:old_stop_t_idx), observ(:,old_start_t_idx:old_stop_t_idx));
            new_bridge_lhood = sum(new_pt.lhood(new_start_t_idx:new_stop_t_idx));
            old_bridge_lhood = sum(old_pt.lhood(old_start_t_idx:old_stop_t_idx));
            
            % Calculate MH acceptance probability
            acc_prob = (new_bridge_lhood+new_param_prob+new_cptime_prob+new_past_cptime_prob+new_past_param_prob+new_past_bridge_lhood) ...
                      -(old_bridge_lhood+old_param_prob+old_cptime_prob+old_past_cptime_prob+old_past_param_prob+old_past_bridge_lhood) ...
                      +(old_replaced_lhood+old_replaced_param_prob+old_past_replaced_lhood) ...
                      -(new_replaced_lhood+new_replaced_param_prob+new_past_replaced_lhood);
            
            % Test for acceptance
            if log(rand) < acc_prob
                
                % Keep proposed particle
                old_pt = new_pt;
                
            end
            
        end
        
%         % Do an MCMC move on the future! Like a sort of smoothing RM step
%         if times(k) < old_pt.tau(end-1)
%             
%             % Now do an MH move on the next changepoint
%             old_start_tau_idx = find_nearest(old_pt.tau, t, false);
%             old_stop_tau_idx = find_nearest(old_pt.tau, t, true);
%             
%             last_tau = old_pt.tau(old_start_tau_idx);
%             next_tau = old_pt.tau(old_stop_tau_idx);
%             future_tau = old_pt.tau(old_stop_tau_idx+1);
%             
%             last_x = old_pt.x(:,old_start_tau_idx);
%             last_w = old_pt.w(:,old_start_tau_idx);
%             
%             % tau proposal
%             tau_replace = inf;
%             while (tau_replace>future_tau-params.dt)||(tau_replace<times(k))
%                 tau_replace = normrnd(next_tau, params.ppsl_move_time_sd);
%             end
%             ppsl_tau_prob = log(normpdf(tau_replace, next_tau, params.ppsl_move_time_sd));
%             
%             % Propose a move in the last accelerations
%             [w_replace, ppsl_w_prob, ~] = acceleration_proposal(flags, params, last_tau, last_x, last_w, times(times<tau_replace), observ(:,times<tau_replace));
%             
%             % Calculate state at this jump time
%             x_replace = next_state(flags, params, last_x, w_replace, tau_replace-last_tau);
%             
%             % Calculate reverse proposal probabilities
%             rev_ppsl_tau_prob = log(normpdf(next_tau, tau_replace, params.ppsl_move_time_sd));
%             [~, ~, rev_ppsl_w_prob] = acceleration_proposal(flags, params, last_tau, last_x, last_w, times(times<next_tau), observ(:,times<next_tau));
%             
%             % Proposal probabilities
%             ppsl_prob = ppsl_tau_prob+ppsl_w_prob;
%             rev_ppsl_prob = rev_ppsl_tau_prob + rev_ppsl_w_prob;
%             
%             % Build new pt
%             new_pt = old_pt;
%             new_pt.tau(old_stop_tau_idx) = tau_replace;
%             new_pt.x(:,old_stop_tau_idx) = x_replace;
%             new_pt.w(:,old_start_tau_idx) = w_replace;
%             new_pt.w(:,old_stop_tau_idx) = determine_acceleration(flags, params, new_pt.x(:,old_stop_tau_idx), new_pt.x(:,old_stop_tau_idx+1), new_pt.tau(old_stop_tau_idx+1)-new_pt.tau(old_stop_tau_idx));
%             
%             % Find the indexes for the new particle
%             start_tau_idx = old_start_tau_idx;
%             stop_tau_idx = old_stop_tau_idx;
%             
%             % Calculate jump transition probabilities
%             [~, new_jump1] = sample_jump_time(flags, params, new_pt.tau(start_tau_idx), t, new_pt.tau(stop_tau_idx));
%             [~, old_jump1] = sample_jump_time(flags, params, old_pt.tau(start_tau_idx), t, old_pt.tau(stop_tau_idx));
%             [~, new_jump2] = sample_jump_time(flags, params, new_pt.tau(start_tau_idx), t, new_pt.tau(stop_tau_idx));
%             [~, old_jump2] = sample_jump_time(flags, params, old_pt.tau(start_tau_idx), t, old_pt.tau(stop_tau_idx));
%             old_jump_trans_prob = old_jump1+old_jump2;
%             new_jump_trans_prob = new_jump1+new_jump2;
%             
%             % Calculate accelerations and associated probabilities
%             new_accel_prob = log(mvnpdf(new_pt.w(:,start_tau_idx)', zeros(1,dr), params.Q)) + log(mvnpdf(new_pt.w(:,stop_tau_idx)', zeros(1,dr), params.Q));
%             old_accel_prob = log(mvnpdf(old_pt.w(:,start_tau_idx)', zeros(1,dr), params.Q)) + log(mvnpdf(old_pt.w(:,stop_tau_idx)', zeros(1,dr), params.Q));
%             
%             % Calculate bridging likelihoods
%             start_t_idx = find_nearest(times, new_pt.tau(start_tau_idx), true);
%             stop_t_idx = find_nearest(times, new_pt.tau(stop_tau_idx), false);
%             [new_pt.intx(:,start_t_idx:stop_t_idx), new_pt.lhood(start_t_idx:stop_t_idx)] = ...
%                 interpolate_state(flags, params, new_pt.tau(start_tau_idx), new_pt.x(:,start_tau_idx), new_pt.w(:,start_tau_idx), times(start_t_idx:stop_t_idx), observ(:,start_t_idx:stop_t_idx));
%             [old_pt.intx(:,start_t_idx:stop_t_idx), old_pt.lhood(start_t_idx:stop_t_idx)] = ...
%                 interpolate_state(flags, params, old_pt.tau(start_tau_idx), old_pt.x(:,start_tau_idx), old_pt.w(:,start_tau_idx), times(start_t_idx:stop_t_idx), observ(:,start_t_idx:stop_t_idx));
%             new_bridge_lhood = sum(new_pt.lhood(start_t_idx:stop_t_idx));
%             old_bridge_lhood = sum(old_pt.lhood(start_t_idx:stop_t_idx));
%             
%             start_t_idx = find_nearest(times, new_pt.tau(start_tau_idx+1), true);
%             stop_t_idx = find_nearest(times, new_pt.tau(stop_tau_idx+1), false);
%             [new_pt.intx(:,start_t_idx:stop_t_idx), new_pt.lhood(start_t_idx:stop_t_idx)] = ...
%                 interpolate_state(flags, params, new_pt.tau(start_tau_idx+1), new_pt.x(:,start_tau_idx+1), new_pt.w(:,start_tau_idx+1), times(start_t_idx:stop_t_idx), observ(:,start_t_idx:stop_t_idx));
%             [old_pt.intx(:,start_t_idx:stop_t_idx), old_pt.lhood(start_t_idx:stop_t_idx)] = ...
%                 interpolate_state(flags, params, old_pt.tau(start_tau_idx+1), old_pt.x(:,start_tau_idx+1), old_pt.w(:,start_tau_idx+1), times(start_t_idx:stop_t_idx), observ(:,start_t_idx:stop_t_idx));
%             new_bridge_lhood = new_bridge_lhood + sum(new_pt.lhood(start_t_idx:stop_t_idx));
%             old_bridge_lhood = old_bridge_lhood + sum(old_pt.lhood(start_t_idx:stop_t_idx));
%             
%             % Calculate MH acceptance probability
%             acc_prob = (new_bridge_lhood+new_accel_prob+new_jump_trans_prob) ...
%                 -(old_bridge_lhood+old_accel_prob+old_jump_trans_prob) ...
%                 +(rev_ppsl_prob+ppsl_prob);
%             
%             % Test for acceptance
%             if log(rand) < acc_prob
%                 
%                 % Keep proposed particle
%                 old_pt = new_pt;
%                 
%             end
%             
%         end
        
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

