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
    
    % Append an extra jump just after the data finishes
    pt.cp.Ns = pt.cp.Ns + 1;
    pt.cp.tau = [pt.cp.tau, params.T+params.dt/10];
    pt.cp.u = [pt.cp.u, 0];
    pt.cp.m = [pt.cp.m, 0];
    
%     % Initialise arrays for the backward filter
%     back_mu = zeros(ds,K);
%     back_P = zeros(ds,ds,K);
%     back_mu(:,K) = squeeze(pt.mu(:,end));
%     back_P(:,:,K) = squeeze(pt.P(:,:,end));
    back_lam = zeros(ds,K);
    back_Ups = zeros(ds,ds,K);
    back_norm = zeros(1,K);
    back_lam(:,K) = (params.H'/params.R)*observs(:,K);
    back_Ups(:,:,K) = (params.H'/params.R)*params.H;
    back_norm(K) = det(2*pi*params.R)^(-1/2) * exp(-( (observs(:,K)'/params.R)*observs(:,K) )/2);
        
    % Loop through time backwards
    for k = K-1:-1:2
        
        cp = pt.cp;
        
        t = times(k);
        next_t = times(k+1);
        Np = length(filt_part_sets{k});
        
        % Find latest jump before t+1 and earliest jump after t
        [prev_jump, prev_ji] = max(cp.tau(cp.tau<next_t));
        [next_jump] = min(cp.tau(cp.tau>t)); next_ji = find(cp.tau==next_jump);
        
%         % Get k+1 backward estimated mean and covariance
%         next_b_mu = back_mu(:,k+1);
%         next_b_P = squeeze(back_P(:,:,k+1));
        next_b_lam = back_lam(:,k+1);
        next_b_Ups = squeeze(back_Ups(:,:,k+1));
        next_b_norm = back_norm(k+1);
        
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
        
%         [pred_b_mu, pred_b_P] = kf_predict(next_b_mu, next_b_P, inv(A), A\Q/A');
%         [b_mu, b_P] = kf_update(pred_b_mu, pred_b_P, observs(:,k), params.H, params.R);
%         
%         % Store
%         back_mu(:,k) = b_mu;
%         back_P(:,:,k) = (b_P+b_P')/2;
        
        invQ = inv(Q);
        pred_b_Ups = A'*( invQ-(invQ/(next_b_Ups+invQ))*invQ )*A;
        pred_b_lam = A'*( (invQ/(next_b_Ups+invQ))*next_b_lam );
        pred_b_norm = next_b_norm * (det( inv(next_b_Ups+invQ) )/det( Q ))^(1/2) * exp(-( next_b_lam'*( inv(next_b_Ups+invQ) )*next_b_lam )/2);
        
        b_lam = pred_b_lam + (params.H'/params.R)*observs(:,k);
        b_Ups = pred_b_Ups + (params.H'/params.R)*params.H;
        b_norm = pred_b_norm * det(2*pi*params.R)^(-0.5) *exp(-( (observs(:,k)'/params.R)*observs(:,k) )/2);
        
        assert(~isinf(b_norm));
        assert(~isnan(b_norm));
        
        back_lam(:,k) = b_lam;
        back_Ups(:,:,k) = b_Ups;
        back_norm(k) = b_norm;
        
        % Fetch forward KF results
        forw_means = cell2mat(arrayfun(@(x) {x.mu(:,k)'}, filt_part_sets{k}))';
        forw_covs = cell2mat(permute(arrayfun(@(x) {x.P(:,:,k)}, filt_part_sets{k}), [2 3 1]));
        
%         % Linear bit
%         sum_P = bsxfun(@plus, f_P, pred_b_P);
%         linear_weights = log_mvnpdf_fast_batch(f_mu, pred_b_mu, sum_P);

        linear_weights = zeros(Np,1);
        % Linear bit
        for jj = 1:Np
            
            f_m = forw_means(:,jj);
            f_P = squeeze(forw_covs(:,:,jj));
            
%             nu = f_m'*( inv(f_P) - f_P\inv(pred_b_Ups+inv(f_P))/f_P )*f_m - 2*( (pred_b_lam'/(pred_b_Ups+inv(f_P)))/f_P )*f_m - (pred_b_lam'/(pred_b_Ups+inv(f_P))*pred_b_lam);
%             
%             linear_weights(jj) = log( sqrt( det( inv(pred_b_Ups+inv(f_P)) )/det( f_P ) ) ) - nu/2;
            
            Gam = chol(f_P + eps*eye(params.state_dim))';
            temp1 = Gam'*(pred_b_lam - pred_b_Ups*f_m);
            temp2 = Gam'*pred_b_Ups*Gam + eye(params.state_dim);
            nu = f_m'*pred_b_Ups*f_m - 2*pred_b_lam'*f_m - temp1'*(temp2\temp1);
            
            linear_weights(jj) = -0.5*log(det(temp2)) - nu/2;
            
        end
        
        % Nonlinear bit
        Np = length(filt_part_sets{k});
        nonlinear_weights = zeros(Np, 1);
        for jj = 1:Np
            [prev_jump, prev_ji] = max(filt_part_sets{k}(jj).cp.tau(filt_part_sets{k}(jj).cp.tau<next_t));
            [ ~, ~, ~, cp_trans_prob ] = sample_jump_time( flags, params, prev_jump, t, next_jump, filt_part_sets{k}(jj).cp.m(prev_ji), filt_part_sets{k}(jj).cp.u(prev_ji) );
            nonlinear_weights(jj) = cp_trans_prob;
        end
        
        % Construct weights
        cond_weights = filt_weight_sets{k} + linear_weights;
        cond_weights = cond_weights - logsumexp(cond_weights);
        
        % Sample
        ind = randsample(length(cond_weights), 1, true, exp(cond_weights));
        
        % Find jump sequence cut points
        [~, stop_past_ji] = max(filt_part_sets{k}(ind).cp.tau(filt_part_sets{k}(ind).cp.tau<t));
        if ~isempty(stop_past_ji), stop_past_ji = filt_part_sets{k}(ind).cp.Ns; end
        
        % Update trajectory
        pt.cp.tau = [filt_part_sets{k}(ind).cp.tau(1:stop_past_ji), pt.cp.tau(next_ji:end)];
        pt.cp.u = [filt_part_sets{k}(ind).cp.u(1:stop_past_ji), pt.cp.u(next_ji:end)];
        pt.cp.m = [filt_part_sets{k}(ind).cp.m(1:stop_past_ji), pt.cp.m(next_ji:end)];
%         pt.mu = [filt_part_sets{k}(ind).mu(:,1:k-1), pt.mu(:,k:end)];
%         pt.P = cat(3, filt_part_sets{k}(ind).P(:,:,1:k-1), pt.P(:,:,k:end));
        pt.cp.Ns = length(pt.cp.tau);
        
        assert(all( ismember(pt.cp.tau,unique(pt.cp.tau)) ));
        
    end
    
    % Kalman smooth
    A_arr = zeros(ds,ds,K);
    Q_arr = zeros(ds,ds,K);
    mu_arr = zeros(ds,K);
    P_arr = zeros(ds,ds,K);
    last_mu = params.start_state;
    last_P = params.start_var;
    
    for k = 1:K
        
        % Kalman Filter
        [ mu, P, lhood, A, Q ] = KF_kinematic_state( flags, params, k, pt.cp, last_mu, last_P, times, observs  );
        
        % Store
        if k>1
            A_arr(:,:,k-1) = A;
            Q_arr(:,:,k-1) = Q;
        end
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


