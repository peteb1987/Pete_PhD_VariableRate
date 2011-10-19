function [ part_sets ] = vr_inference( flags, params, times, observ )
%VR_INFERENCE Run variable rate particle filter and smoother algorithms

% params is a structure with all the model and algorithm parameters
% times is a vector specifying the time grid
% observ is a vector/matrix of observations, first index is time

% The first index in any particle array is the particle index.

% Offset time so first observation occurs at t=0
times = times - times(1);

% Set local variables for commonly-used parameters
ds = params.state_dim;      % State dimensionality
do = params.obs_dim;        % Observation dimensionality
Np = params.Np;             % Number of particles (this is overwritten in each iteration of the loop)

% Create an array to store particle sets for each time frame
part_sets = cell(params.T, 1);

% Count the obervations
K = numel(times);
T = times(K);
assert(numel(observ)==K);

% Initialise the particle set
last_pts_weights = log(ones(Np, 1)/Np);             % Particle weights
last_pts_Ns = ones(Np, 1);                         % Number of states per particle
last_pts_tau = zeros(Np, 1);                        % Particle jump times
last_pts_x = zeros(Np, 1, ds);                      % Particle states
last_pts_mu = zeros(Np, 1, ds);                     % Particle means (for rb case)
last_pts_sigma = zeros(Np, 1, ds, ds);              % Particle covariances (for rb case)
last_pts_intx = zeros(Np, K, ds);                   % States interpolated at observation times
last_pts_intmu = zeros(Np, K, ds);                  % Means interpolated at observation times (for rb case)
last_pts_intsigma = zeros(Np, K, ds, ds);           % Covariances interpolated at observation times (for rb case)

% Initialise first points - PUT THIS AWAY IN A FUNCTION - ITS NOT GENERAL
x_init = [observ(1), zeros(1,ds-do)];                      % Initialise value as first observation padded with zeros (i.e. partial linear observation)
last_pts_x(:,1,:) = repmat(x_init, Np, 1);
last_pts_mu(:,1,:) = repmat(x_init, Np, 1);
last_pts_sigma(:,1,:,:) = permute(repmat([params.x_start_sigma, 0; 0, params.xdot_start_var], [1, 1, Np]), [3,1,2]);

% Loop through observations
for k = 1:K
    
    % Create particle arrays
    pts_Ns = last_pts_Ns;
    pts_tau = last_pts_tau;
    pts_x = last_pts_x;
    pts_mu = last_pts_mu;
    pts_sigma = last_pts_sigma;
    pts_intx = last_pts_intx;
    pts_intmu = last_pts_intmu;
    pts_intsigma = last_pts_intsigma;
    
    % Count the particles
    Np = size(pts_tau, 1);
    
    % Resample/Jump proposal loop: Loop through particles
    for ii = 1:Np
        
        % Calculate number of children of this particle
        Ni = max(1, floor(Np*last_pts_weights(ii)));
        
        % Loop through children
        for jj = 1:Ni
            
            % Sample next jump time to see if one has happened since t-1
            
            %%% CONTINUE FROM HERE
            
            next_tau = sample_next_state(flags, params, pts_tau(ii, pts_Ns(ii)));
            
            
            % If a jump has occured, add it to the state and update weight
            
        end
            
        % Update weight for non-jumping particles
        
    end
    
    % Weighting loop
    for ii = 1:Np
        
        % Calculate predictive likelihood
        
        % Update weights
        
    end
    
end

% Store particle output


end

