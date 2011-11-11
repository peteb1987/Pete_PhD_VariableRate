function [ part_sets ] = vr_filter( flags, params, times, observ )
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

% Count the obervations
K = numel(times);
T = times(K);
assert(size(observ,2)==K);

% Initialise the particle set
MNJ = 2;                                    % MAXIMUM NUMBER OF JUMPS - this is a parameter for array creation
last_pts_weights = log(ones(Np, 1)/Np);             % Particle weights
last_pts_Ns = ones(Np, 1);                          % Number of states per particle
last_pts_tau = zeros(Np, MNJ);                      % Particle jump times
last_pts_type = zeros(Np, MNJ);                     % Particle jump types
last_pts_intx = zeros(Np, 1, ds);                   % States interpolated at observation times
last_pts_intmu = zeros(Np, 1, ds);                  % Means interpolated at observation times (for rb case)
last_pts_intP = zeros(Np, 1, ds, ds);               % Covariances interpolated at observation times (for rb case)

% Initialise first points
[last_pts_x, last_pts_mu, last_pts_P, last_pts_w] = initialise_state(flags, params, Np, MNJ, observ);

% Initialise interpolated variables
last_pts_intx(:,1,:) = last_pts_x(:,1,:);
last_pts_intmu(:,1,:) = last_pts_mu(:,1,:);
last_pts_intP(:,1,:,:) = last_pts_P(:,1,:,:);
last_t = 0;
last_Np = Np;
last_MNJ = MNJ;

% Store
part_sets{1}.Np = last_Np;
part_sets{1}.pts_weights = last_pts_weights;
part_sets{1}.pts_Ns = last_pts_Ns;
part_sets{1}.pts_tau = last_pts_tau;
part_sets{1}.pts_type = last_pts_type;
part_sets{1}.pts_x = last_pts_x;
part_sets{1}.pts_w = last_pts_w;
part_sets{1}.pts_mu = last_pts_mu;
part_sets{1}.pts_P = last_pts_P;
part_sets{1}.pts_intx = last_pts_intx;
part_sets{1}.pts_intmu = last_pts_intmu;
part_sets{1}.pts_intP = last_pts_intP;

% Loop through observations
for k = 2:K
    
    t = times(k);
    
    % Create particle arrays
    MNJ = max(last_MNJ, 2 + max(last_pts_Ns));
    pts_weights = log(ones(Np, 1)/Np);             % Particle weights
    pts_Ns = ones(Np, 1);                          % Number of states per particle
    pts_tau = zeros(Np, MNJ);                      % Particle jump times
    pts_type = NaN(Np, MNJ);                       % Particle jump types
    pts_x = zeros(Np, MNJ, ds);                    % Particle states
    pts_w = zeros(Np, MNJ, dr);                    % Random variables
    pts_mu = zeros(Np, MNJ, ds);                   % Particle means (for rb case)
    pts_P = zeros(Np, MNJ, ds, ds);                % Particle covariances (for rb case)
    pts_intx = zeros(Np, k, ds);                   % States interpolated at observation times
    pts_intmu = zeros(Np, k, ds);                  % Means interpolated at observation times (for rb case)
    pts_intP = zeros(Np, k, ds, ds);               % Covariances interpolated at observation times (for rb case)
    
    % New particle counter
    jj = 0;
    
    % Calculate number of children for each particle
    round_weights = max(1/last_Np, exp(last_pts_weights));
    round_weights(last_pts_weights-max(last_pts_weights)<-5000)=0;
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
        
        last_intx = squeeze(last_pts_intx(ii,k-1,:));
        last_intmu = squeeze(last_pts_intmu(ii,k-1,:));
        last_intP = squeeze(last_pts_intP(ii,k-1,:,:));
        
        non_jumping_kids = 0;
        
        % Loop through children
        for ch = 1:Ni
            
            % Sample next jump time to see if one has happened since t-1
            [tau, type, w, x, mu, P ] = sample_next_state( flags, params, last_t, last_tau, last_x, last_w, last_intmu, last_intP );
            
            % If a jump has occured, add it to the state and update weight
            if tau < t
                
                jj = jj + 1;
                
                % Copy old particle bits accross
                pts_Ns(jj) = last_pts_Ns(ii);
                pts_tau(jj,1:last_MNJ) = last_pts_tau(ii,:);
                pts_type(jj,1:last_MNJ) = last_pts_type(ii,:);
                pts_x(jj,1:last_MNJ,:) = last_pts_x(ii,:,:);
                pts_w(jj,1:last_MNJ,:) = last_pts_w(ii,:,:);
                pts_mu(jj,1:last_MNJ,:) = last_pts_mu(ii,:,:);
                pts_P(jj,1:last_MNJ,:,:) = last_pts_P(ii,:,:,:);
                pts_intx(jj,1:k-1,:) = last_pts_intx(ii,:,:);
                pts_intmu(jj,1:k-1,:) = last_pts_intmu(ii,:,:);
                pts_intP(jj,1:k-1,:,:) = last_pts_intP(ii,:,:,:);
                
                % Add new bits
                pts_weights(jj) = last_pts_weights(ii)-log(Ni);
                pts_Ns(jj) = pts_Ns(jj) + 1;
                pts_tau(jj,pts_Ns(jj)) = tau;
                pts_type(jj,pts_Ns(jj)) = type;
                pts_x(jj,pts_Ns(jj),:) = x;
                pts_w(jj,pts_Ns(jj),:) = w;
                pts_mu(jj,pts_Ns(jj),:) = mu;
                pts_P(jj,pts_Ns(jj),:,:) = P;
                
                [pred_lhood, pts_intx(jj,k,:), pts_intmu(jj,k,:), pts_intP(jj,k,:,:)] = interp_and_lhood(flags, params, tau, t, last_w, x, mu, P, observ(:,k));
                
%                 assert(all(eig(squeeze(pts_intP(jj,k,:,:)))>=-1E-6))
                
                % Update weight
                pts_weights(jj) = pts_weights(jj) + pred_lhood;
                
            else
                
                % Count it as a non-jumping particle
                non_jumping_kids = non_jumping_kids + 1;
                
            end
            
        end
        
        
        if non_jumping_kids>0
            
            % Copy a non-jumping copy of the particle and update weight
            jj = jj + 1;
            pts_weights(jj) = log(non_jumping_kids)+last_pts_weights(ii)-log(Ni);
            pts_Ns(jj) = last_pts_Ns(ii);
            pts_tau(jj,1:last_MNJ) = last_pts_tau(ii,:);
            pts_type(jj,1:last_MNJ) = last_pts_type(ii,:);
            pts_x(jj,1:last_MNJ,:) = last_pts_x(ii,:,:);
            pts_w(jj,1:last_MNJ,:) = last_pts_w(ii,:,:);
            pts_mu(jj,1:last_MNJ,:) = last_pts_mu(ii,:,:);
            pts_P(jj,1:last_MNJ,:,:) = last_pts_P(ii,:,:,:);
            pts_intx(jj,1:k-1,:) = last_pts_intx(ii,:,:);
            pts_intmu(jj,1:k-1,:) = last_pts_intmu(ii,:,:);
            pts_intP(jj,1:k-1,:,:) = last_pts_intP(ii,:,:,:);
            
            [pred_lhood, pts_intx(jj,k,:), pts_intmu(jj,k,:), pts_intP(jj,k,:,:)] = interp_and_lhood(flags, params, last_tau, t, last_w, last_x, last_intmu, last_intP, observ(:,k));

%             assert(all(eig(squeeze(pts_intP(jj,k,:,:)))>=-1E-6))
            
            % Update weight
            pts_weights(jj) = pts_weights(jj) + pred_lhood;
            
            % See if its a while since there was last a state
            if (rand<1)&&((t-last_tau)>(params.dt*20))
                
                % Propose a resample-move step to adjust the last w
                [w, ppsl_prob, rev_ppsl_prob] = tracking_acceleration_proposal(flags, params, last_x, last_tau, last_w, times(1:k), observ(:,1:k));
                
                % Calculate likelihoods
                [old_lhood, ~] = tracking_calc_likelihood(flags, params, last_x, last_tau, last_w, times(1:k), observ(:,1:k));
                [new_lhood, ppsl_intx] = tracking_calc_likelihood(flags, params, last_x, last_tau, w, times(1:k), observ(:,1:k));
                old_accel_prob = log(mvnpdf(last_w', zeros(1,params.rnd_dim), params.Q));
                new_accel_prob = log(mvnpdf(w', zeros(1,params.rnd_dim), params.Q));
                
                % MH acceptance
                acc_prob = (new_lhood+new_accel_prob)-(old_lhood+old_accel_prob)+(rev_ppsl_prob-ppsl_prob);
                if rand<acc_prob
                    start_idx = find(min(times(times>last_tau))==times);
                    pts_w(jj,pts_Ns(jj),:) = w;
                    pts_intx(jj,start_idx:k,:) = ppsl_intx(1,start_idx:k,:);
                end
                
            end
            
        end
        
    end
    
    % Chop off excess columns;
    Np = jj;
    pts_weights(Np+1:end)=[];
    pts_Ns(Np+1:end)=[];
    pts_tau(Np+1:end,:)=[];
    pts_type(Np+1:end,:)=[];
    pts_x(Np+1:end,:,:)=[];
    pts_w(Np+1:end,:,:)=[];
    pts_mu(Np+1:end,:,:)=[];
    pts_P(Np+1:end,:,:,:)=[];
    pts_intx(Np+1:end,:,:)=[];
    pts_intmu(Np+1:end,:,:)=[];
    pts_intP(Np+1:end,:,:,:)=[];
    
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
    last_pts_type = pts_type;
    last_pts_x = pts_x;
    last_pts_w = pts_w;
    last_pts_mu = pts_mu;
    last_pts_P = pts_P;
    last_pts_intx = pts_intx;
    last_pts_intmu = pts_intmu;
    last_pts_intP = pts_intP;
    
    part_sets{k}.Np = Np;
    part_sets{k}.pts_weights = pts_weights;
    part_sets{k}.pts_Ns = pts_Ns;
    part_sets{k}.pts_tau = pts_tau;
    part_sets{k}.pts_type = pts_type;
    part_sets{k}.pts_x = pts_x;
    part_sets{k}.pts_w = pts_w;
    part_sets{k}.pts_mu = pts_mu;
    part_sets{k}.pts_P = pts_P;
    part_sets{k}.pts_intx = pts_intx;
    part_sets{k}.pts_intmu = pts_intmu;
    part_sets{k}.pts_intP = pts_intP;
    
    % Output
    fprintf('*** Completed frame %d, at time %4.3f, using %d particles.\n', k, t, Np);
    
end

end

