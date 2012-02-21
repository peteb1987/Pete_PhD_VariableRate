function [ part_sets, weight_sets ] = rb_vr_filter( flags, params, times, observs )
%RB_VR_FILTER Run RB variable rate particle filter with model-switching for tracking

% flags is a structure with all the model and algorithm flags
% params is a structure with all the model and algorithm parameters
% times is a vector specifying the time grid
% observs is a vector/matrix of observations, first index is time

% Set local variables for commonly-used parameters
Np = params.Np;             % Number of particles (this is overwritten in each iteration of the loop)

% Create an array to store particle sets for each time frame
part_sets = cell(params.K, 1);
weight_sets = cell(params.K, 1);

% Count the obervations
K = numel(times);
T = times(K);
assert(size(observs,2)==K);

% Initialise the particle set and weights
pts = initialise_particles(flags, params, Np);
weights = log(ones(Np, 1)/Np);

% Store the initial particles
part_sets{1} = pts;

last_pts = pts;
last_weights = weights;
last_t = 0;
last_Np = params.Np;

% Loop through observations
for k = 1:K
    
    t = times(k);
    
    % New particle counter
    jj = 0;
    
    % Calculate number of children for each particle
    round_weights = max(1/last_Np, exp(last_weights));
%     round_weights(last_pts_weights-max(last_pts_weights)<-2000)=0;
    round_weights = round_weights/sum(round_weights);
    Nchild = systematic_resample(round_weights, params.Np);
%     Nchild = max(1, floor(last_Np*exp(last_pts_weights)));
    
    % Loop through particles
    for ii = 1:last_Np
        
        % How many children for this particle?
        Ni = Nchild(ii);
        
        % Get state of particle from last time frame
        if k == 1
            last_kin = last_pts(ii).kin0;
        else
            last_kin = last_pts(ii).kin;
        end
        
        last_cp = last_pts(ii).cp;
        last_Ns = last_cp.Ns;
        
        non_jumping_kids = 0;
        
        % Loop through children
        for ch = 1:Ni
            
            % If using RM, need to reset "last_" values
            
            % Sample next jump time to see if one has happened since t-1
            [tau, m, u, trans_prob] = sample_jump_time( flags, params, last_cp.tau(:,last_Ns), last_t);
            
            % If a jump has occured, add it to the state and update weight
            if tau < t
                
                jj = jj + 1;
                
                % Copy old particle accross
                pts(jj) = last_pts(ii);
                
                % Auxiliary weight
                weights(jj) = last_weights(ii)-log(Ni);
                
%                 % Resample move - propose change to accelerations
%                 if flags.resam_move
%                     pts(jj) = resample_move(flags, params, k-1, pts(jj), times, observs, 1);
%                 end
%                 last_w = pts(jj).w(:,last_Ns);
%                 x = next_state(flags, params, last_x, last_w, tau-last_tau);
%                 [w, w_prob] = sample_acceleration(flags, params);
                
                % Add new bits
                Ns = last_Ns + 1;
                pts(jj).cp.Ns = Ns;
                pts(jj).cp.tau(Ns) = tau;
                pts(jj).cp.m(:,Ns) = m;
                pts(jj).cp.u(:,Ns) = u;
                
                % KF and calculate likelihood
                [pts(jj).kin.mu(:,k), pts(jj).kin.P(:,k), pred_lhood] = KF_kinematic_state(flags, params, pts(jj).cp, last_kin, observs(:,k));
                
                % Update weight
                weights(jj) = weights(jj) + pred_lhood;
                
            else
                
                % Count it as a non-jumping particle
                non_jumping_kids = non_jumping_kids + 1;
                
            end
            
        end
        
%         for ch = 1:non_jumping_kids
        if non_jumping_kids > 0
            
            % If using RM, need to reset "last_" values
            
            jj = jj + 1;
            
            % Copy old particle accross
            pts(jj) = last_pts(ii);
            
            % Auxiliary weight
%             weights(jj) = last_weights(ii)-log(Ni);
            weights(jj) = log(non_jumping_kids)+last_weights(ii)-log(Ni);

            % Interpolate state and calculate likelihood
            [pts(jj).mu(:,k), pts(jj).P(:,k), pred_lhood] = KF_kinematic_state(flags, params, k, pts(jj), times, observs);

            % Update weight
            weights(jj) = weights(jj) + pred_lhood;
            
%             % If more than one non-jumping child, rejeuvenate with resample move.
%             if (ch>1)&&(pts(jj).Ns>1)
%                 
%                 % Run a RM step
%                 if flags.resam_move
%                     pts(jj) = resample_move( flags, params, k, pts(jj), times, observs, 2 );
%                 end
%                 
%             end
            
        end
        
    end
    
    % Normalise weights
    lin_weights = exp(weights-max(weights)); assert(all(isreal(lin_weights)));
    weights = log(lin_weights/sum(lin_weights));
    
    % Store particles
    part_sets{k} = pts;
    weight_sets{k} = weights;
    
    % Store for next time
    last_pts = pts;
    last_weights = weights;
    last_t = times(k-1);
    last_Np = Np;
    
    % Output
    fprintf('*** Completed frame %d, at time %4.3f, using %d particles.\n', k, t, Np);
    
end

end

