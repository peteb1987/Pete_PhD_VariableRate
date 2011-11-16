function [ part_sets ] = vr_filter_tracking( flags, params, times, observ )
%VR_INFERENCE Run variable rate particle filter and smoother algorithms

% params is a structure with all the model and algorithm parameters
% times is a vector specifying the time grid
% observ is a vector/matrix of observations, first index is time

% The first index in any particle array is the particle index.
% The second index in any particle array is the observation time (for interpolated variables) or the jump time

% Offset time so first observation occurs at t=0
times = times - times(1);

% Set local variables for commonly-used parameters
ds = params.state_dim;      % State dimensionality
do = params.obs_dim;        % Observation dimensionality
dr = params.rnd_dim;        % Random variable dimensionality
Np = params.Np;             % Number of particles (this is overwritten in each iteration of the loop)

% Create an array to store particle sets for each time frame
part_sets = cell(params.K, 1);
MH_type_counts = zeros(4,1);

% Count the obervations
K = numel(times);
T = times(K);
assert(size(observ,2)==K);

% Initialise the particle set
MNJ = 2;                                    % MAXIMUM NUMBER OF JUMPS - this is a parameter for array creation
last_pts_weights = log(ones(Np, 1)/Np);             % Particle weights
last_pts_Ns = ones(Np, 1);                          % Number of states per particle
last_pts_tau = zeros(Np, MNJ);                      % Particle jump times
last_pts_intx = zeros(Np, 1, ds);                   % States interpolated at observation times
last_pts_lhood = zeros(Np, 1, 1);                  % Likelihoods of each observation given current states and accelerations

% Initialise first points
[last_pts_x, last_pts_w] = initialise_state_tracking(flags, params, Np, MNJ, observ);

% Initialise interpolated variables
last_pts_intx(:,1,:) = last_pts_x(:,1,:);
last_t = 0;
last_Np = Np;
last_MNJ = MNJ;

% Store
part_sets{1}.Np = last_Np;
part_sets{1}.pts_weights = last_pts_weights;
part_sets{1}.pts_Ns = last_pts_Ns;
part_sets{1}.pts_tau = last_pts_tau;
part_sets{1}.pts_x = last_pts_x;
part_sets{1}.pts_w = last_pts_w;
part_sets{1}.pts_intx = last_pts_intx;
part_sets{1}.pts_lhood = last_pts_lhood;

% Loop through observations
for k = 2:K
    
    t = times(k);
    
    % Create particle arrays
    MNJ = max(last_MNJ, 2 + max(last_pts_Ns));
    pts_weights = log(ones(Np, 1)/Np);             % Particle weights
    pts_Ns = ones(Np, 1);                          % Number of states per particle
    pts_tau = zeros(Np, MNJ);                      % Particle jump times
    pts_x = zeros(Np, MNJ, ds);                    % Particle states
    pts_w = zeros(Np, MNJ, dr);                    % Random variables
    pts_intx = zeros(Np, k, ds);                   % States interpolated at observation times
    pts_lhood = zeros(Np, k, 1);
    
    % New particle counter
    jj = 0;
    
    % Calculate number of children for each particle
    round_weights = max(1/last_Np, exp(last_pts_weights));
    %     round_weights(last_pts_weights-max(last_pts_weights)<-2000)=0;
    round_weights = round_weights/sum(round_weights);
    Nchild = systematic_resample(round_weights, params.Np);
    %     Nchild = max(1, floor(last_Np*exp(last_pts_weights)));
    
    % Loop through particles
    for ii = 1:last_Np
        
        % Calculate number of children of this particle
        Ni = Nchild(ii);
        % %         if (last_Np>params.Np)&&(Ni==1)&&(rand>(last_Np*exp(last_pts_weights(ii))))
        %         if (Ni==1)&&(rand>(last_Np*exp(last_pts_weights(ii))))
        %             Ni=0;
        %         end
        
        last_tau = last_pts_tau(ii, last_pts_Ns(ii));
        last_x = squeeze(last_pts_x(ii, last_pts_Ns(ii),:));
        last_w = squeeze(last_pts_w(ii, last_pts_Ns(ii),:));
        
        non_jumping_kids = 0;
        
        % Loop through children
        for ch = 1:Ni
            
            % Sample next jump time to see if one has happened since t-1
            [tau, w, x ] = sample_next_state_tracking( flags, params, last_t, last_tau, last_x, last_w );
            
            % If a jump has occured, add it to the state and update weight
            if tau < t
                
                jj = jj + 1;
                
                % Copy old particle bits accross
                pts_Ns(jj) = last_pts_Ns(ii);
                pts_tau(jj,1:last_MNJ) = last_pts_tau(ii,:);
                pts_x(jj,1:last_MNJ,:) = last_pts_x(ii,:,:);
                pts_w(jj,1:last_MNJ,:) = last_pts_w(ii,:,:);
                pts_intx(jj,1:k-1,:) = last_pts_intx(ii,:,:);
                pts_lhood(jj,1:k-1,:) = last_pts_lhood(ii,:,:);
                
                % Run a RM step
                if 1%%&&(rand<0.25)&&((t-last_tau)>(params.dt*30))
                    
                    kk = k-1;
                    
                    % Propose a resample-move step to adjust the last w
                    [w_replace, ppsl_prob, rev_ppsl_prob] = tracking_acceleration_proposal(flags, params, last_x, last_tau, last_w, times(1:kk), observ(:,1:kk));
                    
                    % Calculate likelihoods
                    start_idx = find(min(times(times>last_tau))==times);
                    old_lhood = sum(pts_lhood(jj,start_idx:end));
                    [new_lhood, ppsl_intx] = tracking_calc_likelihood(flags, params, last_x, last_tau, w_replace, times(1:kk), observ(:,1:kk));
                    old_accel_prob = log(mvnpdf(last_w', zeros(1,params.rnd_dim), params.Q));
                    new_accel_prob = log(mvnpdf(w_replace', zeros(1,params.rnd_dim), params.Q));
                    
                    % MH acceptance
                    acc_prob = (sum(new_lhood)+new_accel_prob)-(sum(old_lhood)+old_accel_prob)+(rev_ppsl_prob-ppsl_prob);
                    if log(rand)<acc_prob
                        pts_w(jj,pts_Ns(jj),:) = w_replace;
                        pts_intx(jj,start_idx:kk,:) = ppsl_intx(1,start_idx:kk,:);
                        pts_lhood(jj,start_idx:kk,:) = new_lhood;
                        x = tracking_calc_next_state(flags, last_x, tau-last_tau, w_replace);
                        MH_type_counts(1) = MH_type_counts(1)+1;
                    end
                    
                end
                
                % Add new bits
                pts_weights(jj) = last_pts_weights(ii)-log(Ni);
                pts_Ns(jj) = pts_Ns(jj) + 1;
                pts_tau(jj,pts_Ns(jj)) = tau;
                pts_x(jj,pts_Ns(jj),:) = x;
                pts_w(jj,pts_Ns(jj),:) = w;
                
                [pred_lhood, pts_intx(jj,k,:)] = interp_and_lhood_tracking(flags, params, tau, t, last_w, x, observ(:,k));
                pts_lhood(jj,k,:) = pred_lhood;
                
                %                 assert(all(eig(squeeze(pts_intP(jj,k,:,:)))>=-1E-6))
                
                % Update weight
                pts_weights(jj) = pts_weights(jj) + pred_lhood;
                
            else
                
                % Count it as a non-jumping particle
                non_jumping_kids = non_jumping_kids + 1;
                
            end
            
        end
        
        
        if non_jumping_kids>0
            
            for ch = 1:non_jumping_kids
                
                % Copy a non-jumping copy of the particle and update weight
                jj = jj + 1;
                %                 pts_weights(jj) = log(non_jumping_kids)+last_pts_weights(ii)-log(Ni);
                pts_weights(jj) = last_pts_weights(ii)-log(Ni);
                pts_Ns(jj) = last_pts_Ns(ii);
                pts_tau(jj,1:last_MNJ) = last_pts_tau(ii,:);
                pts_x(jj,1:last_MNJ,:) = last_pts_x(ii,:,:);
                pts_w(jj,1:last_MNJ,:) = last_pts_w(ii,:,:);
                pts_intx(jj,1:k-1,:) = last_pts_intx(ii,:,:);
                pts_lhood(jj,1:k-1,:) = last_pts_lhood(ii,:,:);
                
                [pred_lhood, pts_intx(jj,k,:)] = interp_and_lhood_tracking(flags, params, last_tau, t, last_w, last_x, observ(:,k));
                pts_lhood(jj,k,:) = pred_lhood;
                
                %             assert(all(eig(squeeze(pts_intP(jj,k,:,:)))>=-1E-6))
                
                % Update weight
                pts_weights(jj) = pts_weights(jj) + pred_lhood;
                
                % Run a RM step
                if (ch>1)%&&((t-last_tau)>(params.dt*20))%%&&(rand<0.05)
                    
                    kk = k;
                    
                    % Randomly choose one of two types of MH step
                    if (pts_Ns(jj)>1)&&(rand<0.5)
                        
                        tau_penult = pts_tau(jj,pts_Ns(jj)-1);
                        
                        % Propose a move in the last jump time
                        tau_replace = normrnd(last_tau, 0.2/(params.rate_shape*params.rate_scale));
                        tau_replace = min(tau_replace, t-params.dt);
                        
                        % Propose a resample-move step to adjust the last w
                        [w_replace, ppsl_prob, rev_ppsl_prob] = tracking_acceleration_proposal(flags, params, last_x, tau_replace, last_w, times(1:kk), observ(:,1:kk));
                        
                        ppsl_prob = ppsl_prob + log(normpdf(tau_replace, last_tau, 0.2/(params.rate_shape*params.rate_scale)));
                        rev_ppsl_prob = rev_ppsl_prob + log(normpdf(last_tau, tau_replace, 0.2/(params.rate_shape*params.rate_scale)));
                        
                        % Calculate likelihoods
                        start_idx = find(min(times(times>tau_replace))==times);
                        old_lhood = sum(pts_lhood(jj,start_idx:end));
                        [new_lhood, ppsl_intx] = tracking_calc_likelihood(flags, params, last_x, tau_replace, w_replace, times(1:kk), observ(:,1:kk));
                        old_accel_prob = log(mvnpdf(last_w', zeros(1,params.rnd_dim), params.Q));
                        new_accel_prob = log(mvnpdf(w_replace', zeros(1,params.rnd_dim), params.Q));
                        [ new_trans_prob ] = tracking_calc_jump_trans_prob( params, tau_penult, tau_replace ) + tracking_calc_jump_trans_prob( params, tau_replace );
                        [ old_trans_prob ] = tracking_calc_jump_trans_prob( params, tau_penult, last_tau ) + tracking_calc_jump_trans_prob( params, last_tau );
                        
                        % MH acceptance
                        acc_prob = (sum(new_lhood)+new_accel_prob+new_trans_prob)-(sum(old_lhood)+old_accel_prob+old_trans_prob)+(rev_ppsl_prob-ppsl_prob);
                        if log(rand)<acc_prob
                            pts_tau(jj,pts_Ns(jj)) = tau_replace;
                            pts_w(jj,pts_Ns(jj),:) = w_replace;
                            pts_intx(jj,start_idx:kk,:) = ppsl_intx(1,start_idx:kk,:);
                            pts_lhood(jj,start_idx:kk,:) = new_lhood;
                            MH_type_counts(2) = MH_type_counts(2)+1;
                        end
                        
                    else
                        
                        if (last_tau+params.dt<t)&&(rand<0.5)
                            
                            % Add a jump
                            
                            % First find the reverse probabilities
                            start_idx = find(min(times(times>last_tau))==times);
                            [~, ~, rev_ppsl_prob] = tracking_acceleration_proposal(flags, params, last_x, last_tau, last_w, times(1:kk), observ(:,1:kk));
                            old_accel_prob = log(mvnpdf(last_w', zeros(1,params.rnd_dim), params.Q));
                            old_trans_prob = tracking_calc_jump_trans_prob( params, last_tau );
                            old_lhood = sum(pts_lhood(jj,start_idx:end));
                            
                            % Now add a jump
                            extra_tau = unifrnd(last_tau+params.dt, t);
                            stop_idx = find(max(times(times<extra_tau))==times);
                            [w1, ppsl_prob_part1, ~] = tracking_acceleration_proposal(flags, params, last_x, last_tau, last_w, times(1:stop_idx), observ(:,1:stop_idx));
                            extra_x = tracking_calc_next_state(flags, last_x, extra_tau-last_tau, w1);
                            [w2, ppsl_prob_part2, ~] = tracking_acceleration_proposal(flags, params, extra_x, extra_tau, last_w, times(1:kk), observ(:,1:kk));
                            ppsl_prob = ppsl_prob_part1 + ppsl_prob_part2 + log(unifpdf(extra_tau, last_tau, t));
                            [new_lhood_part1, ppsl_intx_part1] = tracking_calc_likelihood(flags, params, last_x, last_tau, w1, times(1:stop_idx), observ(:,1:stop_idx));
                            [new_lhood_part2, ppsl_intx_part2] = tracking_calc_likelihood(flags, params, extra_x, extra_tau, w2, times(1:kk), observ(:,1:kk));
                            ppsl_intx = ppsl_intx_part2; ppsl_intx(1,start_idx:stop_idx,:) = ppsl_intx_part1(1,start_idx:stop_idx,:);
                            new_lhood = sum(new_lhood_part1) + sum(new_lhood_part2);
                            new_trans_prob = tracking_calc_jump_trans_prob( params, last_tau, extra_tau ) + tracking_calc_jump_trans_prob( params, extra_tau );
                            new_accel_prob = log(mvnpdf(w1', zeros(1,params.rnd_dim), params.Q)) + log(mvnpdf(w2', zeros(1,params.rnd_dim), params.Q));
                            
                            % MH acceptance
                            acc_prob = (sum(new_lhood)+new_accel_prob+new_trans_prob)-(sum(old_lhood)+old_accel_prob+old_trans_prob)+(rev_ppsl_prob-ppsl_prob);
                            if log(rand)<acc_prob
                                pts_Ns(jj) = pts_Ns(jj)+1;
                                pts_w(jj,pts_Ns(jj)-1,:) = w1;
                                pts_w(jj,pts_Ns(jj),:) = w2;
                                pts_x(jj,pts_Ns(jj),:) = extra_x;
                                pts_tau(jj,pts_Ns(jj),:) = extra_tau;
                                pts_intx(jj,start_idx:kk,:) = ppsl_intx(1,start_idx:kk,:);
                                pts_lhood(jj,start_idx:kk,:) = new_lhood;
                                last_x = extra_x;
                                last_tau = extra_tau;
                                last_w = w2;
                                MH_type_counts(3) = MH_type_counts(3)+1;
                            end
                            
                        elseif (pts_Ns(jj)>1)
                            
                            % Remove the last jump
                            
                            % First find the probabilities for the backwards move
                            tau_penult = pts_tau(jj,pts_Ns(jj)-1);
                            x_penult = squeeze(pts_x(jj,pts_Ns(jj)-1,:));
                            w_penult = squeeze(pts_w(jj,pts_Ns(jj)-1,:));
                            start_idx = find(min(times(times>tau_penult))==times);
                            stop_idx = find(max(times(times<last_tau))==times);
                            [~, ~, rev_ppsl_prob_part1] = tracking_acceleration_proposal(flags, params, x_penult, tau_penult, w_penult, times(1:stop_idx), observ(:,1:stop_idx));
                            [~, ~, rev_ppsl_prob_part2] = tracking_acceleration_proposal(flags, params, last_x, last_tau, last_w, times(1:kk), observ(:,1:kk));
                            rev_ppsl_prob = rev_ppsl_prob_part1 + rev_ppsl_prob_part2 + log(unifrnd(tau_penult,t));
                            old_lhood = sum(pts_lhood(jj,start_idx:end));
                            old_accel_prob = log(mvnpdf(w_penult', zeros(1,params.rnd_dim), params.Q)) + log(mvnpdf(last_w', zeros(1,params.rnd_dim), params.Q));
                            [ old_trans_prob ] = tracking_calc_jump_trans_prob( params, tau_penult, last_tau ) + tracking_calc_jump_trans_prob( params, last_tau );
                            
                            
                            % Now replace it
                            [w_replace, ppsl_prob, ~] = tracking_acceleration_proposal(flags, params, x_penult, tau_penult, w_penult, times(1:stop_idx), observ(:,1:stop_idx));
                            [new_lhood, ppsl_intx] = tracking_calc_likelihood(flags, params, x_penult, tau_penult, w_replace, times(1:kk), observ(:,1:kk));
                            new_accel_prob = log(mvnpdf(w_replace', zeros(1,params.rnd_dim), params.Q));
                            [ new_trans_prob ] = tracking_calc_jump_trans_prob( params, tau_penult);
                            
                            % MH acceptance
                            acc_prob = (sum(new_lhood)+new_accel_prob+new_trans_prob)-(sum(old_lhood)+old_accel_prob+old_trans_prob)+(rev_ppsl_prob-ppsl_prob);
                            if log(rand)<acc_prob
                                pts_Ns(jj) = pts_Ns(jj)-1;
                                pts_w(jj,pts_Ns(jj)+1,:) = zeros(dr,1);
                                pts_x(jj,pts_Ns(jj)+1,:) = zeros(ds,1);
                                pts_tau(jj,pts_Ns(jj)+1,:) = 0;
                                pts_w(jj,pts_Ns(jj),:) = w_replace;
                                pts_intx(jj,start_idx:kk,:) = ppsl_intx(1,start_idx:kk,:);
                                pts_lhood(jj,start_idx:kk,:) = new_lhood;
                                MH_type_counts(4) = MH_type_counts(4)+1;
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
    % Chop off excess columns;
    Np = jj;
    pts_weights(Np+1:end)=[];
    pts_Ns(Np+1:end)=[];
    pts_tau(Np+1:end,:)=[];
    pts_x(Np+1:end,:,:)=[];
    pts_w(Np+1:end,:,:)=[];
    pts_intx(Np+1:end,:,:)=[];
    pts_lhood(Np+1:end,:,:)=[];
    
    % Normalise weights
    lin_weights = exp(pts_weights-max(pts_weights)); assert(all(isreal(lin_weights)));
    pts_weights = log(lin_weights/sum(lin_weights));
    
    % Store particle output
    last_MNJ = MNJ;
    last_t = t;
    last_Np = Np;
    last_pts_weights = pts_weights;
    last_pts_Ns = pts_Ns;
    last_pts_tau = pts_tau;
    last_pts_x = pts_x;
    last_pts_w = pts_w;
    last_pts_intx = pts_intx;
    last_pts_lhood = pts_lhood;
    
    part_sets{k}.Np = Np;
    part_sets{k}.pts_weights = pts_weights;
    part_sets{k}.pts_Ns = pts_Ns;
    part_sets{k}.pts_tau = pts_tau;
    part_sets{k}.pts_x = pts_x;
    part_sets{k}.pts_w = pts_w;
    part_sets{k}.pts_intx = pts_intx;
    part_sets{k}.pts_lhood = pts_lhood;
    
    % Output
    fprintf('*** Completed frame %d, at time %4.3f, using %d particles.\n', k, t, Np);
    
end

MH_type_counts

end

