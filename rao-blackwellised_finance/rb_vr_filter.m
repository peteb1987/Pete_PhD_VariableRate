function [ part_sets, weight_sets ] = rb_vr_filter( flags, params, times, observ )
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
Np = params.Np;             % Number of particles (this is overwritten in each iteration of the loop)

% Create an array to store particle sets for each time frame
part_sets = cell(params.K, 1);
weight_sets = cell(params.K, 1);

% Count the obervations
K = numel(times);
assert(size(observ,2)==K);

% Initialise particles
[pts, weights] = initialise_particles(flags, params, params.Np, observ);
last_pts = pts;
last_weights = weights;
last_Np = Np;

% Store initial particles
part_sets{1} = pts;
weight_sets{1} = weights;

% Loop through observations
for k = 2:K
    
    t = times(k);
    last_t = times(k-1);
    
    % New particle counter
    jj = 0;
    weights = zeros(last_Np,1);
    
    % Calculate number of children for each particle
    round_weights = max(1/last_Np, exp(last_weights));
    round_weights(last_weights-max(last_weights)<-2000)=0;
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

        last_tau = last_pts(ii).tau(last_pts(ii).Ns);
        last_intmu = last_pts(ii).intmu(:,k-1);
        last_intP = squeeze(last_pts(ii).intP(:,:,k-1));
        
        non_jumping_kids = 0;
        
        % Loop through children
        for ch = 1:Ni
            
            % Sample next jump time to see if one has happened since t-1
            [tau, type] = sample_next_state( flags, params, last_t, last_tau);
            
            % If a jump has occured, add it to the state and update weight
            if tau < t
                
                jj = jj + 1;
                
                % Copy particle
                pts(jj) = last_pts(ii);
                
                % Add new bits
                weights(jj) = last_weights(ii)-log(Ni);
                pts(jj).Ns = pts(jj).Ns + 1;
                pts(jj).tau(pts(jj).Ns) = tau;
                pts(jj).type(pts(jj).Ns) = type;
                
                [pred_lhood, pts(jj).intmu(:,k), pts(jj).intP(:,:,k)] = interp_and_lhood(flags, params, last_intmu, last_intP, observ(:,k), t-last_t, type);
                
%                 assert(all(eig(squeeze(pts_intP(jj,k,:,:)))>=-1E-6))
                
                % Update weight
                weights(jj) = weights(jj) + pred_lhood;
                
            else
                
                % Count it as a non-jumping particle
                non_jumping_kids = non_jumping_kids + 1;
                
            end
            
        end
        
        
        if non_jumping_kids>0
            
            % Collapse to a single non-jump particle
            jj = jj + 1;
            
            % Copy particle
            pts(jj) = last_pts(ii);
            
            weights(jj) = log(non_jumping_kids)+last_weights(ii)-log(Ni);
            
            [pred_lhood, pts(jj).intmu(:,k), pts(jj).intP(:,:,k)] = interp_and_lhood(flags, params, last_intmu, last_intP, observ(:,k), t-last_t, 0);
            
            %             assert(all(eig(squeeze(pts_intP(jj,k,:,:)))>=-1E-6))
            
            % Update weight
            weights(jj) = weights(jj) + pred_lhood;
            
        end
        
    end
    
    % Update Np
    Np = jj;
    weights(Np+1:end) = [];
    pts(Np+1:end) = [];
    
    % Normalise weights
    lin_weights = exp(weights-max(weights)); assert(all(isreal(lin_weights)));
    weights = log(lin_weights/sum(lin_weights));
    
    % Store particle output
    part_sets{k} = pts;
    weight_sets{k} = weights;
    last_pts = pts;
    last_Np = Np;
    last_weights = weights;
    
    % Output
    fprintf('*** Completed frame %d, at time %4.3f, using %d particles.\n', k, t, Np);
    
end

% Final frame resample
Nchild = systematic_resample(lin_weights, Np);
parent = zeros(Np,1); cnt=0;
for ii = 1:Np
    parent(cnt+1:cnt+Nchild(ii)) = ii;
    cnt = cnt + Nchild(ii);
end
weight_sets{K} = log(ones(Np,1)/Np);
part_sets{K} = part_sets{K}(parent);

end

