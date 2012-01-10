function [ smooth_pts] = rb_vr_smoother( flags, params, times, observ, filt_part_sets, filt_weight_sets)
%VR_SMOOTHER Takes the variable rate filtering output and runs a particle
%smoothing algorithm on it.

% Set local variables for parameters
K = params.K; S = params.S;
ds = params.state_dim; do = params.obs_dim;
F = params.F; C = params.C; H = params.H; R = params.R;

% Resample S particles from final filtering distribution
Nchild = systematic_resample(exp(filt_weight_sets{K}), S);
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
for ii = 1:S
    
    pt = smooth_pts(ii);
    
    % Get trajectory from the array - append an extra jump just after the data finishes
    pt.Ns = pt.Ns + 1;
    pt.tau = [pt.tau, params.T];
    pt.type = [pt.type, 0];
    
    % Initialise arrays for the backward filter
    back_intmu = zeros(ds,K);
    back_intP = zeros(ds,ds,K);
    back_intmu(:,K) = squeeze(pt.intmu(:,end));
    back_intP(:,:,K) = squeeze(pt.intP(:,:,end));
        
    % Loop through time backwards
    for k = K-1:-1:2
        
        t = times(k);
        next_t = times(k+1);
        Np = length(filt_part_sets{k});
        
        % Find latest jump before t+1 and earliest jump after t
        [prev_jump, prev_ji] = max(pt.tau(pt.tau<next_t));
        [next_jump] = min(pt.tau(pt.tau>t)); next_ji = find(pt.tau==next_jump);
        
        % Get k+1 backward estimated mean and covariance
        next_b_mu = back_intmu(:,k+1);
        next_b_P = squeeze(back_intP(:,:,k+1));
        
        % Run backwards filter prediction
        [A, Q] = lti_disc(F,eye(2),C,next_t-t);
        if (prev_jump>t)&&(pt.type(prev_ji)==1)
            Q = Q + [params.x_start_sd^2 0; 0 0];
        elseif (prev_jump>t)&&(pt.type(prev_ji)==2)
            Q = Q + [0 0; 0 params.xdot_start_sd^2];
        end
        [pred_b_mu, pred_b_P] = kf_predict(next_b_mu, next_b_P, inv(A), A\Q/A');
        [b_mu, b_P] = kf_update(pred_b_mu, pred_b_P, observ(:,k)', H, R);
        
        % Store
        back_intmu(:,k) = b_mu';
        back_intP(:,:,k) = (b_P+b_P')/2;
        
        % Construct weights
        cond_weights = zeros(Np,1);
        for jj = 1:Np
            
            % Fetch forward KF results
            f_mu = filt_part_sets{k}(jj).intmu(:,k);
            f_P = filt_part_sets{k}(jj).intP(:,:,k);
            
            % Linear part
            cond_weights(jj) = cond_weights(jj) + log(mvnpdf(f_mu', b_mu', f_P+b_P));
            
            % Nonlinear part not required because we're using exponentially
            % distributed jump times, so there is no dependence on the
            % past.
            
        end
        
        % Sample
        ind = randsample(length(cond_weights), 1, true, exp(cond_weights));
        
        % Find jump sequence cut points
        [~, stop_past_ji] = max(filt_part_sets{k}(ind).tau(filt_part_sets{k}(ind).tau<t));
        if ~isempty(stop_past_ji), stop_past_ji = filt_part_sets{k}(ind).Ns; end
        
        % Update trajectory
        pt.tau = [filt_part_sets{k}(ind).tau(1:stop_past_ji), pt.tau(next_ji:end)];
        pt.type = [filt_part_sets{k}(ind).type(1:stop_past_ji), pt.type(next_ji:end)];
        pt.intmu = [filt_part_sets{k}(ind).intmu(:,1:k-1), pt.intmu(:,k:end)];
        pt.P = cat(3, filt_part_sets{k}(ind).intP(:,:,1:k-1), pt.intP(:,:,k:end));
        pt.Ns = length(pt.tau);
        
    end
    
    % Kalman smooth
    A_arr = zeros(ds,ds,K);
    Q_arr = zeros(ds,ds,K);
    mu_arr = zeros(ds,K);
    P_arr = zeros(ds,ds,K);
    last_t = 0;
    last_mu = [observ(:,1)'; 0];
    last_P = [params.x_start_sd^2, 0; 0, params.xdot_sd^2];
    mu_arr(:,1) = last_mu;
    P_arr(:,:,1) = last_P;
    
    ji = 2;
    for k = 2:K
        
        t = times(k);
        
        % Create transition matrices
        [A, Q] = lti_disc(F, eye(2), C, t-last_t);
        
        % See if a jump happened
        if (ji<length(pt.tau))&&(t>pt.tau(ji))
            if pt.type(ji)==1
                Q = Q + [params.x_jump_sd^2, 0; 0, 0];
            elseif pt.type(ji)==2
                Q = Q + [0, 0; 0, params.xdot_jump_sd^2];
            end
            ji = ji + 1;
        end
                
        % Store
        A_arr(:,:,k-1) = A;
        Q_arr(:,:,k-1) = Q;
        mu_arr(:,k) = last_mu;
        P_arr(:,:,k) = last_P;
        
        % Kalman Filter
        [pred_m, pred_P] = kf_predict(last_mu, last_P, A, Q);
        [last_mu, last_P] = kf_update(pred_m, pred_P, observ(:,k)', H, R);

        last_t = t;
        
    end
    
    [intmu, intP] = rts_smooth(mu_arr, P_arr, A_arr, Q_arr);
    
    smooth_pts(ii).Ns = pt.Ns;
    smooth_pts(ii).tau = pt.tau;
    smooth_pts(ii).type = pt.type;
    smooth_pts(ii).intmu = intmu;
    smooth_pts(ii).intP = intP;
    
    % Output
    fprintf('*** Completed trajectory %d.\n', ii);
    
end

end


