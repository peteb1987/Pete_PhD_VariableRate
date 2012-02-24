function [ pts ] = initialise_particles(flags, params, Np)
%INITIALISE_PARTICLES Set up particles for particle filtering

ds = params.state_dim;

cp = struct('tau', 0,...                    % Changepoint times
             'm', 1,...                     % Model following each changepoint
             'u', 0,...                     % Model parameters following each changepoint
             'Ns', 1);                      % Number of changepoints

pts = struct('mu', zeros(ds,0),...           % State mean at each observation time
             'P', zeros(ds,ds,0),...          % State covariance at each obsrvation time
             'mu0', params.start_state,...  % Prior state mean at time 0
             'P0', params.start_var,...     % Prior state covariance at time 0
             'cp', repmat({cp}, [Np,1]));

end