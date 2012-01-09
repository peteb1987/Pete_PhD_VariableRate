function [ smooth_part_sets] = rb_vr_smoother( flags, params, times, observ, filt_part_sets)
%VR_SMOOTHER Takes the variable rate filtering output and runs a particle
%smoothing algorithm on it.

% Set local variables for parameters
K = params.K; S = params.S;
ds = params.state_dim; do = params.obs_dim;
F = params.F; C = params.C; H = params.H; R = params.R;

% Initialise particle set array
smooth_part_sets = cell(K,1);

% Initialise Kth set
% Initialise array
max_Ns = max(filt_part_sets{K}.pts_Ns);
smooth_part_sets{K}.pts_Ns = zeros(S, 1);
smooth_part_sets{K}.pts_tau = zeros(S, max_Ns);
smooth_part_sets{K}.pts_type = zeros(S, max_Ns);
% smooth_part_sets{K}.pts_mu = zeros(S, max_Ns, ds);
% smooth_part_sets{K}.pts_P = zeros(S, max_Ns, ds, ds);
smooth_part_sets{K}.pts_intmu = zeros(S, K, ds);
smooth_part_sets{K}.pts_intP = zeros(S, K, ds, ds);

Nchild = systematic_resample(exp(filt_part_sets{K}.pts_weights), S);
Ncum = 0;

% Loop through filter particles
for ii = 1:filt_part_sets{K}.Np
    
    % How many offspring?
    Ni = Nchild(ii);
    
    % Make Ni copies of the particle
    Ns = filt_part_sets{K}.pts_Ns(ii);
    smooth_part_sets{K}.pts_Ns(Ncum+1:Ncum+Ni) = repmat(filt_part_sets{K}.pts_Ns(ii), [Ni, 1]);
    smooth_part_sets{K}.pts_tau(Ncum+1:Ncum+Ni,1:Ns) = repmat(filt_part_sets{K}.pts_tau(ii,1:Ns), [Ni, 1]);
    smooth_part_sets{K}.pts_type(Ncum+1:Ncum+Ni,1:Ns) = repmat(filt_part_sets{K}.pts_type(ii,1:Ns), [Ni, 1]);
%     smooth_part_sets{K}.pts_mu(Ncum+1:Ncum+Ni,1:Ns,:) = repmat(filt_part_sets{K}.pts_mu(ii,1:Ns,:), [Ni, 1, 1]);
%     smooth_part_sets{K}.pts_P(Ncum+1:Ncum+Ni,1:Ns,:,:) = repmat(filt_part_sets{K}.pts_P(ii,1:Ns,:,:), [Ni, 1, 1, 1]);
    smooth_part_sets{K}.pts_intmu(Ncum+1:Ncum+Ni,:,:) = repmat(filt_part_sets{K}.pts_intmu(ii,:,:), [Ni, 1, 1]);
    smooth_part_sets{K}.pts_intP(Ncum+1:Ncum+Ni,:,:,:) = repmat(filt_part_sets{K}.pts_intP(ii,:,:,:), [Ni, 1, 1, 1]);
    
    Ncum = Ncum + Ni;
    
end

% Loop through trajectories
for ii = 1:S

    % Get trajectory from the array - append an extra jump just after the data finishes
    Ns = smooth_part_sets{K}.pts_Ns(ii);
    tau = squeeze(smooth_part_sets{K}.pts_tau(ii,1:Ns));
    type = squeeze(smooth_part_sets{K}.pts_type(ii,1:Ns));
%     mu = permute(smooth_part_sets{K}.pts_mu(ii,1:Ns,:),[2,3,1]);
%     P = permute(smooth_part_sets{K}.pts_P(ii,1:Ns,:,:),[2,3,4,1]);
%     intmu = squeeze(smooth_part_sets{K}.pts_intmu(ii,:,:));
%     intP = squeeze(smooth_part_sets{K}.pts_intP(ii,:,:,:));
    
    % Initialise arrays for the backward filter
    back_intmu = zeros(K,ds);
    back_intP = zeros(K,ds,ds);
    back_intmu(K,:) = squeeze(smooth_part_sets{K}.pts_intmu(ii,K,:));
    back_intP(K,:,:) = squeeze(smooth_part_sets{K}.pts_intP(ii,K,:,:));
    
    % Initialise arrays for the forward filter values
    forw_intmu = zeros(K,ds);
    forw_intP = zeros(K,ds,ds);
    
    % Initialise arrays for the combined state estimate
    intmu = zeros(K,ds);
    intP = zeros(K,ds,ds);
    intmu(K,:) = squeeze(smooth_part_sets{K}.pts_intmu(ii,K,:));
    intP(K,:,:) = squeeze(smooth_part_sets{K}.pts_intP(ii,K,:,:));
    
    % Loop through time backwards
    for k = K-1:-1:1
        
        t = times(k);
        
        Np = filt_part_sets{k}.Np;
        
        % Get k+1 backward estimated mean and covariance
        b_mu = back_intmu(k+1,:)';
        b_P = squeeze(back_intP(k+1,:,:));
        
        % Run backwards filter
        tt = t;
        [prev_jump, prev_ji] = max(tau(tau<t));
        while (k>1)&&(tt>prev_jump)&&(prev_jump>times(k-1))
            % Diffuse
            [A, Q] = lti_disc(F,eye(2),C,tt-prev_jump);
            [b_mu, b_P] = kf_predict(b_mu, b_P, inv(A), A\Q/A');
            % Jump
            if type(prev_ji)==1
                % x jump
                b_P = b_P + [params.x_jump_sd^2 0; 0 0];
            elseif type(prev_ji)==2
                % xdot jump
                b_P = b_P + [0 0; 0 params.xdot_jump_sd^2];
            end
            tt = prev_jump;
            prev_ji = prev_ji - 1;
            prev_jump = tau(prev_ji);
        end
        if k>1
            % Diffuse
            [A, Q] = lti_disc(F,eye(2),C,tt-times(k-1));
            [b_mu, b_P] = kf_predict(b_mu, b_P, inv(A), A\Q/A');
        end
        % KF update
        pred_b_mu = b_mu; pred_b_P = b_P;
        [b_mu, b_P] = kf_update(b_mu, b_P, observ(:,k)', H, R);
        % Store
        back_intmu(k,:) = b_mu';
        back_intP(k,:,:) = (b_P+b_P')/2;
        
        % Find forward estimated mean and covariance
        f_mu = permute(filt_part_sets{k}.pts_intmu(:,k,:),[3,1,2]);
        f_P = permute(filt_part_sets{k}.pts_intP(:,k,:,:),[3,4,1,2]);

        % Calculate Gaussian weight term
        comb_P = f_P + repmat(pred_b_P, [1 1 Np]);
%         state_match_prob = log(mvnpdf(f_mu', pred_b_mu', comb_P));
        state_match_prob = log_mvnpdf_fast_batch(f_mu, repmat(pred_b_mu,1,Np), comb_P);
        
        % Calculate jump transition probabilities
        [trans_prob, prev_ji, next_ji] = rb_transition_probability(params, t, filt_part_sets{k}.pts_tau, tau, filt_part_sets{k}.pts_type, type);
        
        % Calculate conditional weights
        conditional_weights = filt_part_sets{k}.pts_weights + trans_prob + state_match_prob;
        conditional_weights = conditional_weights-max(conditional_weights);
        temp = exp(conditional_weights); temp = temp/sum(temp);
        conditional_weights = log(temp);
        
        % Sample
        ind = randsample(length(conditional_weights), 1, true, exp(conditional_weights));
        
        % Update trajectory
        tau = [filt_part_sets{k}.pts_tau(ind,1:prev_ji(ind)), tau(next_ji:end)];
        type = [filt_part_sets{k}.pts_type(ind,1:prev_ji(ind)), type(next_ji:end)];
%         mu = [permute(filt_part_sets{k}.pts_mu(ind,1:prev_ji(ind),:),[2,3,1]); mu(next_ji:end,:)];
%         P = [permute(filt_part_sets{k}.pts_P(ind,1:prev_ji(ind),:,:),[2,3,4,1]); P(next_ji:end,:,:)];
        Ns = length(tau);
        
% %         % Calculate interpolated values
%         f_P_chosen = squeeze(f_P(:,:,ind));
%         f_mu_chosen = squeeze(f_mu(:,ind));
%         forw_intmu(k,:) = f_mu_chosen;
%         forw_intP(k,:,:) = f_P_chosen;
%         av_P = inv(inv(pred_b_P)+inv(f_P_chosen));
%         intmu(k,:) = av_P*(pred_b_P\pred_b_mu+f_P_chosen\f_mu_chosen);
%         intP(k,:,:) = av_P;
        
        % Store trajectory
        smooth_part_sets{k}.pts_tau(ii,1:Ns) = tau;
        smooth_part_sets{k}.pts_type(ii,1:Ns) = type;
%         smooth_part_sets{k}.pts_mu(ii,1:Ns,:) = mu;
%         smooth_part_sets{k}.pts_P(ii,1:Ns,:,:) = P;
        smooth_part_sets{k}.pts_Ns(ii) = Ns;
%         smooth_part_sets{k}.pts_intmu(ii,:,:) = intmu;
%         smooth_part_sets{k}.pts_intP(ii,:,:,:) = intP;
        
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
        if (ji<length(tau))&&(t>tau(ji))
            if type(ji)==1
                Q = Q + [params.x_jump_sd^2, 0; 0, 0];
            elseif type(ji)==2
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
    
    smooth_part_sets{1}.pts_intmu(ii,:,:) = intmu;
    smooth_part_sets{1}.pts_intP(ii,:,:,:) = intP;
    
    % Output
    fprintf('*** Completed trajectory %d.\n', ii);
    
end

end


