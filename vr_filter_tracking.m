function [ part_sets ] = vr_filter_tracking( flags, params, times, observ )
%VR_INFERENCE Run variable rate particle filter and smoother algorithms

% flags is a structure with all the model and algorithm flags
% params is a structure with all the model and algorithm parameters
% times is a vector specifying the time grid
% observ is a vector/matrix of observations, first index is time

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

% Initialise the particle set and weights
pts = initialise_particles_tracking(flags, params, observ);
weights = log(ones(Np, 1)/Np);

% Store the initial particles
part_sets{1} = pts;

% Loop through observations
for k = 2:K
    
    t = times(k);
    
    % Initialise previous frame values
    last_pts = pts;
    last_weights = weights;
    last_t = times(k-1);
    last_Np = Np;
    
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
        last_Ns = last_pts(ii).Ns;
        last_tau = last_pts(ii).tau(last_Ns);
        last_x = last_pts(ii).x(:,last_Ns);
        last_w = last_pts(ii).w(:,last_Ns);
        
        non_jumping_kids = 0;
        
        % Loop through children
        for ch = 1:Ni
            
            % Sample next jump time to see if one has happened since t-1
            [tau, w, x ] = sample_next_state_tracking( flags, params, last_t, last_tau, last_x, last_w );
            
            % If a jump has occured, add it to the state and update weight
            if tau < t
                
                jj = jj + 1;
                
                % Copy old particle accross
                pts(jj) = last_pts(ii);
                
                % Auxiliary weight
                weights(jj) = last_weights(ii)-log(Ni);
                
                % Resample move - propose change to accelerations
                pts(jj) = resam_move_tracking(flags, params, k-1, pts(jj), times, observ, 1);
                last_w = pts(jj).w(:,last_Ns);
                x = tracking_calc_next_state(flags, last_x, tau-last_tau, last_w);
                
                % Add new bits
                Ns = last_Ns + 1;
                pts(jj).Ns = Ns;
                pts(jj).tau(Ns) = tau;
                pts(jj).x(:,Ns) = x;
                pts(jj).w(:,Ns) = w;
                
                % Interpolate state and calculate likelihood
                [pred_lhood, pts(jj).intx(:,k)] = interp_and_lhood_tracking(flags, params, tau, t, w, x, observ(:,k));
                pts(jj).lhood(k,1) = pred_lhood;
                
                % Update weight
                weights(jj) = weights(jj) + pred_lhood;
                
            else
                
                % Count it as a non-jumping particle
                non_jumping_kids = non_jumping_kids + 1;
                
            end
            
        end
        
        for ch = 1:non_jumping_kids
            
            jj = jj + 1;
            
            % Copy old particle accross
            pts(jj) = last_pts(ii);
            
            % Auxiliary weight
            weights(jj) = last_weights(ii)-log(Ni);
            % weights(jj) = log(non_jumping_kids)+last_weights(ii)-log(Ni);

            % Interpolate state and calculate likelihood
            [pred_lhood, pts(jj).intx(:,k)] = interp_and_lhood_tracking(flags, params, last_tau, t, last_w, last_x, observ(:,k));
            pts(jj).lhood(k,1) = pred_lhood;

            % Update weight
            weights(jj) = weights(jj) + pred_lhood;
            
            % If more than one non-jumping child, rejeuvenate with resample move.
            if (ch>1)&&(pts(jj).Ns>1)
                
                    % Run a RM step
                    pts(jj) = resam_move_tracking( flags, params, k, pts(jj), times, observ, 2 );
                
            end
            
        end
        
    end
    
    % Normalise weights
    lin_weights = exp(weights-max(weights)); assert(all(isreal(lin_weights)));
    weights = log(lin_weights/sum(lin_weights));
    
    % Store particles
    part_sets{k} = pts;
    
    % Output
    fprintf('*** Completed frame %d, at time %4.3f, using %d particles.\n', k, t, Np);
    
end

end

