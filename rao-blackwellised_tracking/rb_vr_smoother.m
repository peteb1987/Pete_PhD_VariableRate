function [ smooth_pts ] = rb_vr_smoother( flags, params, times, observs, filt_part_sets, filt_weight_sets)
%RB_VR_SMOOTHER Takes the variable rate filtering output and runs a particle
%smoothing algorithm on it.

% Set local variables for parameters
K = params.K; Ns = params.Ns;
ds = params.state_dim;

% Resample S particles from final filtering distribution
Nchild = systematic_resample(exp(filt_weight_sets{K}), Ns);
jj = 0;

% Loop through filter particles
for ii = 1:length(filt_part_sets{K})
    
    % How many offspring?
    Ni = Nchild(ii);
    
    for ch = 1:Ni
        
        jj = jj + 1;
        
        % Copy the particle
        smooth_pts(jj,1) = filt_part_sets{K}(ii);
        
    end
    
end

% Loop through trajectories
for ii = 1:Ns
    
    pt = smooth_pts(ii);
    cp = pt.cp;
    
    % Append an extra jump just after the data finishes
    cp.Ns = cp.Ns + 1;
    cp.tau = [cp.tau, params.T+params.dt/10];
    cp.u = [cp.u, 0];
    cp.m = [cp.m, 0];
    
    % Initialise arrays for the backward filter
    back_mu = zeros(ds,K);
    back_P = zeros(ds,ds,K);
    back_mu(:,K) = squeeze(pt.mu(:,end));
    back_P(:,:,K) = squeeze(pt.P(:,:,end));
        
    % Loop through time backwards
    for k = K-1:-1:2
        
        t = times(k);
        next_t = times(k+1);
        Np = length(filt_part_sets{k});
        
        % Find latest jump before t+1 and earliest jump after t
        [prev_jump, prev_ji] = max(cp.tau(cp.tau<next_t));
        [next_jump] = min(cp.tau(cp.tau>t)); next_ji = find(cp.tau==next_jump);
        
        % Get k+1 backward estimated mean and covariance
        next_b_mu = back_mu(:,k+1);
        next_b_P = squeeze(back_P(:,:,k+1));
        
        % Transition matrices
        if next_jump >= next_t
            % No new jump
            [A, Q, ~] = construct_transmats(next_t-t, cp.m(next_ji-1), cp.u(next_ji-1), params.proc_var);
        else
            [A, Q, Ajump] = construct_transmats(next_t-next_jump, cp.m(next_ji), cp.u(next_ji), params.proc_var);
            [Ainc, Qinc, ~] = construct_transmats(next_jump-t, cp.m(next_ji-1), cp.u(next_ji-1), params.proc_var);
            A = A*Ajump*Ainc;
            Q = A*Ajump*Qinc*Ajump'*A' + Q;
        end
        
        [pred_b_mu, pred_b_P] = kf_predict(next_b_mu, next_b_P, inv(A), A\Q/A');
        [b_mu, b_P] = kf_update(pred_b_mu, pred_b_P, observs(:,k), params.H, params.R);
        
        % Store
        back_mu(:,k) = b_mu;
        back_P(:,:,k) = (b_P+b_P')/2;

        % Fetch forward KF results
        f_mu = cell2mat(arrayfun(@(x) {x.mu(:,k)'}, filt_part_sets{k}))';
        f_P = cell2mat(permute(arrayfun(@(x) {x.P(:,:,k)}, filt_part_sets{k}), [2 3 1]));
        
        % Linear bit
        sum_P = bsxfun(@plus, f_P, pred_b_P);
        linear_weights = log_mvnpdf_fast_batch(f_mu, pred_b_mu, sum_P);
        
        % Nonlinear bit
        Np = length(filt_part_sets{k});
        nonlinear_weights = zeros(Np, 1);
        for jj = 1:Np
            [prev_jump, prev_ji] = max(filt_part_sets{k}(jj).cp.tau(filt_part_sets{k}(jj).cp.tau<next_t));
            [ ~, ~, ~, cp_trans_prob ] = sample_jump_time( flags, params, prev_jump, t, next_jump, filt_part_sets{k}(jj).cp.m(prev_ji), filt_part_sets{k}(jj).cp.u(prev_ji) );
            nonlinear_weights(jj) = log(cp_trans_prob);
        end
        
        % Construct weights
        cond_weights = filt_weight_sets{k} + linear_weights;
        
        % Sample
        ind = randsample(length(cond_weights), 1, true, exp(cond_weights));
        
        % Find jump sequence cut points
        [~, stop_past_ji] = max(filt_part_sets{k}(ind).cp.tau(filt_part_sets{k}(ind).cp.tau<t));
        if ~isempty(stop_past_ji), stop_past_ji = filt_part_sets{k}(ind).cp.Ns; end
        
        % Update trajectory
        pt.cp.tau = [filt_part_sets{k}(ind).cp.tau(1:stop_past_ji), pt.cp.tau(next_ji:end)];
        pt.cp.u = [filt_part_sets{k}(ind).cp.u(1:stop_past_ji), pt.cp.u(next_ji:end)];
        pt.cp.m = [filt_part_sets{k}(ind).cp.m(1:stop_past_ji), pt.cp.m(next_ji:end)];
        pt.mu = [filt_part_sets{k}(ind).mu(:,1:k-1), pt.mu(:,k:end)];
        pt.P = cat(3, filt_part_sets{k}(ind).P(:,:,1:k-1), pt.P(:,:,k:end));
        pt.cp.Ns = length(pt.cp.tau);
        
    end
    
    % Kalman smooth
    A_arr = zeros(ds,ds,K);
    Q_arr = zeros(ds,ds,K);
    mu_arr = zeros(ds,K);
    P_arr = zeros(ds,ds,K);
    last_mu = params.start_state;
    last_P = params.start_var;
    
    ji = 2;
    for k = 1:K
        
        t = times(k);
        
        % Kalman Filter
        [ mu, P, lhood ] = KF_kinematic_state( flags, params, k, pt.cp, last_mu, last_P, times, observs  );
        
        % Store
        A_arr(:,:,k) = A;
        Q_arr(:,:,k) = Q;
        mu_arr(:,k) = mu;
        P_arr(:,:,k) = P;
        
        last_mu = mu;
        last_P = P;
        
    end
    
    % RTS smooth
    [mu, P] = rts_smooth(mu_arr, P_arr, A_arr, Q_arr);
    
    pt.mu = mu;
    pt.P = P;
    smooth_pts(ii) = pt;
    
    % Output
    fprintf('*** Completed trajectory %d.\n', ii);
    
end

end


