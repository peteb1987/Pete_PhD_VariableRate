function [ mu ] = observation_mean( flags, params, x )
%OBSERVATION_MEAN Deterministically calculate the expected observation
%given the current state

% Can handle multiple state in columns of x, or multiple values of random
% variables in columns of w, but not both.

Ns = size(x,2);

if flags.obs_mod == 1
    % Direct observation
    mu = x(1:params.obs_dim,:);
    
elseif flags.obs_mod == 2
    
    % Radar observation
    mu = zeros(params.obs_dim,Ns);
    
    % Position - bearing and range
    [mu(1,:), mu(2,:)] = cart2pol(x(1,:), x(2,:));
    
    if flags.obs_vel
        % Calculate polar unit vectors
        er = unit(x(1:2,:),1);
        etheta = [-er(2,:); er(1,:)];
        % Velocity - bearing rate and range rage
        v = x(3:4,:);
        mu(4,:) = dot(v, er);
        mu(3,:) = dot(v, etheta)./mu(2,:);
    end

else
    error('unhandled option');
end

end

