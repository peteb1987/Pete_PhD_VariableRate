function [ pt ] = resample_move( flags, params, k, pt, times, observs, type )
%RESAM_MOVE_TRACKING Conduct resample move on a particle, altering current
%accelerations and last jump time

% type can be:  1 = change to previous accelerations
%               2 = change to previous accelerations and jump time
%               3 = add/remove a jump using RJMCMC

% Get current state of particle
Ns = pt.Ns;
last_tau = pt.tau(Ns);
last_x = pt.x(:,Ns);
last_w = pt.w(:,Ns);

switch type
    
    case 1
        %%% Change only the last accelerations
        
        % Find frame in which to start likelihood calculations/interpolation
        start_idx = find(min(times(times>last_tau))==times);
        
        % Propose a change in the acclerations
        [w_replace, ppsl_prob, rev_ppsl_prob] = acceleration_proposal(flags, params, last_tau, last_x, last_w, times(1:k), observs(:,1:k));
        
        % Calculate likelihoods
        old_lhood = pt.lhood(start_idx:end);
        [new_intx, new_lhood] = interpolate_state(flags, params, last_tau, last_x, w_replace, times(start_idx:k), observs(:,start_idx:k));
        
        % Calculate transition (i.e. acceleration) probabilities
        old_accel_prob = log(mvnpdf(last_w', zeros(1,params.rnd_dim), params.Q));
        new_accel_prob = log(mvnpdf(w_replace', zeros(1,params.rnd_dim), params.Q));
        
        % MH acceptance
        acc_prob = +(sum(new_lhood)+new_accel_prob) ...
                   -(sum(old_lhood)+old_accel_prob) ...
                   +(rev_ppsl_prob-ppsl_prob) ;
               
        if log(rand)<acc_prob
            pt.w(:,Ns) = w_replace;
            pt.intx(:,start_idx:k) = new_intx;
            pt.lhood(start_idx:k) = new_lhood;
            pt.w_prob(Ns) = new_accel_prob;
        end
        
        
    case 2
        %%% Change both the last accelerations and the last state time
        
        penult_tau = pt.tau(Ns-1);
        penult_x = pt.x(:,Ns-1);
        penult_w = pt.w(:,Ns-1);
        
        % Propose a move in the last jump time
        tau_replace = normrnd(last_tau, params.ppsl_move_time_sd);
        ppsl_tau_prob = log(normpdf(tau_replace, last_tau, params.ppsl_move_time_sd));
        
        if (tau_replace>times(k))||(tau_replace<penult_tau)
            % Invalid jump time proposed. Return immediately.
            return
        end
        
        % Calculate state at this jump time
        x_replace = next_state(flags, params, penult_x, penult_w, tau_replace-penult_tau);
        
        % Propose a move in the last accelerations
        [w_replace, ppsl_w_prob, ~] = acceleration_proposal(flags, params, tau_replace, x_replace, last_w, times(1:k), observs(:,1:k));
        
        % Calculate reverse proposal probabilities
        rev_ppsl_tau_prob = log(normpdf(last_tau, tau_replace, params.ppsl_move_time_sd));
        [~, ~, rev_ppsl_w_prob] = acceleration_proposal(flags, params, last_tau, last_x, last_w, times(1:k), observs(:,1:k));
        
        % Proposal probabilities
        ppsl_prob = ppsl_tau_prob+ppsl_w_prob;
        rev_ppsl_prob = rev_ppsl_tau_prob + rev_ppsl_w_prob;
        
        % Find frame in which to start likelihood calculations/interpolation
        last_start_idx = find(min(times(times>last_tau))==times);
        replace_start_idx = find(min(times(times>tau_replace))==times);
        start_idx = min(last_start_idx,replace_start_idx);
        
        % Calculate likelihoods
        old_lhood = pt.lhood(start_idx:k);
        new_intx = zeros(size(pt.intx));
        new_lhood = zeros(size(pt.lhood));
        [new_intx(:,replace_start_idx:k), new_lhood(replace_start_idx:k)] = interpolate_state(flags, params, tau_replace, x_replace, w_replace, times(replace_start_idx:k), observs(:,replace_start_idx:k));
        if replace_start_idx>last_start_idx
            [new_intx(:,last_start_idx:replace_start_idx-1), new_lhood(last_start_idx:replace_start_idx-1)] = interpolate_state(flags, params, penult_tau, penult_x, penult_w, times(last_start_idx:replace_start_idx-1), observs(:,last_start_idx:replace_start_idx-1));
        end
        
        % Calculate transition (i.e. acceleration, and jump time) probabilities
        old_accel_prob = log(mvnpdf(last_w', zeros(1,params.rnd_dim), params.Q));
        new_accel_prob = log(mvnpdf(w_replace', zeros(1,params.rnd_dim), params.Q));
        [~,new_trans1] = sample_jump_time(flags, params, penult_tau, [], tau_replace);
        [~,new_trans2] = sample_jump_time(flags, params, tau_replace, [], []);
        new_trans_prob = new_trans1 + new_trans2;
        [~,old_trans_1] = sample_jump_time(flags, params, penult_tau, [], last_tau);
        [~,old_trans_2] = sample_jump_time(flags, params, last_tau, [], []);
        old_trans_prob = old_trans_1 + old_trans_2;
        
        % MH acceptance
        acc_prob = +(sum(new_lhood)+new_accel_prob+new_trans_prob) ...
            -(sum(old_lhood)+old_accel_prob+old_trans_prob) ...
            +(rev_ppsl_prob-ppsl_prob) ;
        
        if log(rand)<acc_prob
            pt.tau(Ns) = tau_replace;
            pt.w(:,Ns) = w_replace;
            pt.x(:,Ns) = x_replace;
            pt.intx(:,start_idx:k) = new_intx(:,start_idx:k);
            pt.lhood(start_idx:k) = new_lhood(start_idx:k);
            pt.w_prob(Ns) = new_accel_prob;
            [~, pt.tau_prob(Ns)] = sample_jump_time(flags, params, penult_tau, [], tau_replace);
        end
    
end

end

