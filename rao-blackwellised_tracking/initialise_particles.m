function [ pts, weights ] = initialise_particles( flags, params, Np, observ )
%INITIALISE_PARTICLES Initialise particles for rbvrpf

% Starting mean and covariance
mu_init = [observ(1); 0];
P_init = [params.x_start_sd^2, 0; 0 params.xdot_start_sd^2];

% Initialise the particle set
pts = struct('Ns', 1, ...
             'tau', 0, ...
             'type', 0, ...
             'intmu', num2cell(repmat(mu_init,[1,Np]), 1)', ...
             'intP', squeeze(num2cell(repmat(P_init,[1,1,Np]), [1,2])) );

weights = log(ones(Np,1)/Np);

end

