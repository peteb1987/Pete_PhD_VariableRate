function [ rmse, MAP_rmse ] = filter_errors( flags, params, filt_part_sets, times, true_tau, true_intx )
%FILTER_ERRORS Calculate MMSE and MAP errors of filtering estimates

sd = flags.space_dim;
K = length(times);

rmse.pos_over_time = zeros(1, K);
rmse.vel_over_time = zeros(1, K);
MAP_rmse.pos_over_time = zeros(1, K);
MAP_rmse.vel_over_time = zeros(1, K);


for k = 1:K
    
    pts = filt_part_sets{k};
    
    % Count particles
    Np = length(pts);
    
    % Calculate RMSEs for MMSE estimate
    intx = mean(cat(3, pts.intx),3);
    error = abs(true_intx(:,k) - intx(:,k));
    rmse.pos_over_time(k) = sqrt(mean( sum(error(1:sd,:,:).^2,1), 3));
    rmse.vel_over_time(k) = sqrt(mean( sum(error(sd+1:2*sd,:,:).^2,1), 3));
    
    % Find MAP particle
    prob = zeros(Np,1);
    for ii = 1:Np
        prob(ii) = sum(pts(ii).tau_prob(1:pts(ii).Ns)) + sum(pts(ii).w_prob(1:pts(ii).Ns)) + sum(pts(ii).lhood(1:k)); % There should be P(x_0) term in this.
    end
    [~, MAP_ind] = max(prob);
    
    % Calculate RMSEs for MAP estimate
    intx = pts(MAP_ind).intx;
    error = abs(bsxfun(@minus, true_intx(:,k,:), intx(:,k,:)));
    MAP_rmse.pos_over_time = sqrt(mean( sum(error(1:sd,:,:).^2,1), 3));
    MAP_rmse.vel_over_time = sqrt(mean( sum(error(sd+1:2*sd,:,:).^2,1), 3));

end

rmse.pos = sqrt(mean(rmse.pos_over_time.^2));
rmse.vel = sqrt(mean(rmse.vel_over_time.^2));
MAP_rmse.pos = sqrt(mean(MAP_rmse.pos_over_time.^2));
MAP_rmse.vel = sqrt(mean(MAP_rmse.vel_over_time.^2));



% figure, plot(times, rmse_over_time)

end

