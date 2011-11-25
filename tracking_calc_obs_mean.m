function [ mu ] = tracking_calc_obs_mean( flags, params, x, w )
%TRACKING_CALC_OBS_MEAN Calculate the Gaussian mean of the observation
%distribution for various observation models

Ns = size(x,2);

if flags.obs_mod == 1
    % Direct observation
    mu = x(1:params.obs_dim,:);
elseif flags.obs_mod == 2
    % Radar observation
    mu = zeros(params.obs_dim,Ns);
    if flags.dyn_mod < 5
        [mu(1,:), mu(2,:)] = cart2pol(x(1,:), x(2,:));
        if params.obs_dim == 4
            % Bearing and range rate
            if flags.dyn_mod == 2
                xdot = x(4,:).*cos(x(3,:)) + w(3,:);
                ydot = x(4,:).*sin(x(3,:)) + w(4,:);
                [x(3,:),x(4,:)] = cart2pol(xdot, ydot);
            end
            mu(4,:) = x(4,:).*cos(x(3,:)-mu(1,:));
            mu(3,:) = x(4,:).*sin(x(3,:)-mu(1,:))./mu(2,:);
        end
    else
        [mu(1,:), mu(2,:), mu(3,:)] = cart2sph(x(1,:), x(2,:), x(3,:));
        if params.obs_dim == 6
            r = x(1:3,:);
            v = bsxfun(@plus, x(4:6,:), w(4:6,:));
            er = bsxfun(@rdivide, r, sqrt(sum(r.^2)));
            epsi = cross(er, repmat([0;0;1],1,size(r,2)));
            epsi = bsxfun(@rdivide, epsi, sqrt(sum(epsi.^2)));
            etheta = cross(epsi, er);
            rdot = dot(v, er);
            thetadot = dot(v, etheta)./mu(3,:);
            psidot =dot(v, epsi)./(mu(3,:).*sin(mu(2,:)));
            mu(6,:) = rdot;
            mu(5,:) = thetadot;
            mu(4,:) = psidot;
        end
    end
	
end

end

