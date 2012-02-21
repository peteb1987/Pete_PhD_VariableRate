function [ pts ] = initialise_particles(flags, params, Np)
%INITIALISE_PARTICLES Set up particles for particle filtering

% Generate a random set of starting points and accelerations
mu = repmat({params.start_state}, [Np,1]);
P = repmat({params.start_var}, [Np,1]);

cp = struct('tau', 0,...                    % Changepoint times
             'm', 1,...                     % Model following each changepoint
             'u', 1,...                     % Model parameters following each changepoint
             'Ns', 1);                      % Number of changepoints

kin = struct('mu', zeros(6,0),...           % State mean at each observation time
             'P', zeros(6,6,0));            % State covariance at each obsrvation time

kin0 = struct('mu0', mu,...                 % State mean at each observation time
             'P0', P);                      % State covariance at each obsrvation time

% Put the values into cell array
pts = struct('kin0', kin0,...
             'kin', kin,...
             'cp', cp);

end