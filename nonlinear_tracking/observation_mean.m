function [ mu ] = observation_mean( flags, params, x, w )
%OBSERVATION_MEAN Deterministically calculate the expected observation
%given the current state, the random variables using the various
%observation models.

% Can handle multiple state in columns of x, or multiple values of random
% variables in columns of w, but not both.

Ns = size(x,2);

if flags.obs_mod == 1
    % Direct observation
    mu = x(1:params.obs_dim,:);
    
elseif flags.obs_mod == 2
    % Radar observation
    mu = zeros(params.obs_dim,Ns);
    
    if flags.space_dim == 2
        % Position - bearing and range
        [mu(1,:), mu(2,:)] = cart2pol(x(1,:), x(2,:));
        if flags.obs_vel
            % Calculate polar unit vectors
            er = unit(x(1:2,:));
            etheta = [-er(2,:); er(1,:)];
            % Velocity - bearing rate and range rage
            if flags.dyn_mod == 1
                v = x(3:4,:);
            elseif flags.dyn_mod == 2
                v = bsxfun(@plus, x(3:4,:), w(3:4,:));
            else
                error('unhandled option');
            end
            mu(4,:) = dot(v, er);
            mu(3,:) = dot(v, etheta)./mu(2,:);
            
        end
    elseif flags.space_dim == 3
        % Position - bearing, elevation and range
        [mu(1,:), mu(2,:), mu(3,:)] = cart2sph(x(1,:), x(2,:), x(3,:));
        if flags.obs_vel
            % Calculate polar unit vectors
            er = unit(x(1:3,:));
            epsi = cross2(er, [0;0;1]);
            epsi = unit(epsi);
            etheta = cross(epsi, er);

            % Velocity - bearing rate, elevation rate and range rage
            if flags.dyn_mod == 1
                v = x(4:6,:);
            elseif flags.dyn_mod == 2
                v = bsxfun(@plus, x(4:6,:), w(4:6,:));
            else
                error('unhandled option');
            end
            mu(6,:) = dot(v, er);
            mu(5,:) = dot(v, etheta)./mu(3,:);
            mu(4,:) = dot(v, epsi)./(mu(3,:).*cos(mu(2,:)));

        end
    else
        error('unhandled option');
    end
else
    error('unhandled option');
end

end

