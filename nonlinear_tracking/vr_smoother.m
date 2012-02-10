function [ pts ] = vr_smoother( flags, params, filt_part_sets, times, observ )
%MCMC_VR_SMOOTHER Smoother for variable rate models using direct backwards
%conditional sampling

% Set some local variables
K = params.K; S = params.S; Np = params.Np;
ds = params.state_dim; dr = params.rnd_dim; do = params.obs_dim;

% Create a cell array for the smoothed particles
% pts = cell(S,1);
pts = initialise_particles(flags, params, S, observ);

% Loop through smoothing particles
for ii = 1:S
    
    % Initialise particle with a random final frame filtering particle
    pt = filt_part_sets{K}(unidrnd(Np));
    
    % Append a state just after the final time
    pt.Ns = pt.Ns+1;
    pt.tau = [pt.tau, times(end)+eps(times(end))];
    pt.x = [pt.x, pt.intx(:,end)];
    pt.w = [pt.w, zeros(dr,1)];
    pt.tau_prob = [];
    pt.w_prob = [];
    
    % Loop backwards through time
    for k = K:-1:2

        t = times(k);
        
        % Create an array of bacward sampling weights
%         back_weights = zeros(Np,1);
        history_lhood_arr = zeros(Np,1);
        history_trans_prob_arr = zeros(Np,1);
        current_lhood_arr = zeros(Np,1);
        current_trans_prob_arr = zeros(Np,1);
        
        % Find indexes
        stop_tau_idx = find_nearest(pt.tau, t, true);
        stop_t_idx = find_nearest(times, pt.tau(stop_tau_idx), false);
        
        % Loop through filtering particles
        for jj = 1:Np
            
            filt_pt = filt_part_sets{k}(jj);
            
            % Find indexes
            start_tau_idx = find_nearest(filt_pt.tau, t, false);
            start_t_idx = find_nearest(times, filt_pt.tau(start_tau_idx), true);
            
            % Transition probability
            [~, jump_time_prob] = sample_jump_time(flags, params, filt_pt.tau(start_tau_idx), [], pt.tau(stop_tau_idx));
            w = determine_acceleration(flags, params, filt_pt.x(:,start_tau_idx), pt.x(:,stop_tau_idx), pt.tau(stop_tau_idx)-filt_pt.tau(start_tau_idx));
            innov_prob = log(mvnpdf(w', zeros(1,dr), params.Q));
            current_trans_prob_arr(jj) = jump_time_prob + innov_prob;
            
            % Likelihood
            [~, lhood] = interpolate_state(flags, params, filt_pt.tau(start_tau_idx), filt_pt.x(:,start_tau_idx), w, times(start_t_idx:stop_t_idx), observ(:, start_t_idx:stop_t_idx));
            current_lhood_arr(jj) = sum(lhood);
            
            % History probabilities
            history_lhood_arr(jj) = sum(filt_pt.lhood(1:start_t_idx-1));
            history_trans_prob_arr(jj) = sum(filt_pt.tau_prob(1:start_tau_idx)) + sum(filt_pt.w_prob(1:start_tau_idx-1));
            
        end
        
        % Calculate weights
        back_weights = history_lhood_arr + history_trans_prob_arr + current_lhood_arr + current_trans_prob_arr;
        
        % Sample weights
        back_weights = back_weights - max(back_weights);
        jj = randsample(Np, 1, true, exp(back_weights));
        filt_pt = filt_part_sets{k}(jj);
        
        % Construct new particle
        start_tau_idx = find_nearest(filt_pt.tau, t, false);
        start_t_idx = find_nearest(times, filt_pt.tau(start_tau_idx), true);
        pt.tau = [filt_pt.tau(1:start_tau_idx) pt.tau(stop_tau_idx:end)];
        pt.Ns = size(pt.tau,2);
        pt.x = [filt_pt.x(:,1:start_tau_idx) pt.x(:,stop_tau_idx:end)];
        pt.w = [filt_pt.w(:,1:start_tau_idx) pt.w(:,stop_tau_idx:end)];
        w = determine_acceleration(flags, params, filt_pt.x(:,start_tau_idx), pt.x(:,start_tau_idx+1), pt.tau(start_tau_idx+1)-filt_pt.tau(start_tau_idx));
        pt.w(:, start_tau_idx) = w;
        [pt.intx(:,start_t_idx:stop_t_idx), ~] = interpolate_state(flags, params, pt.tau(start_tau_idx), pt.x(:,start_tau_idx), pt.w(:,start_tau_idx), times(start_t_idx:stop_t_idx), observ(:, start_t_idx:stop_t_idx));
        
    end
    
    % Output
    fprintf('*** Completed trajectory %d.\n', ii);
    
end


end

