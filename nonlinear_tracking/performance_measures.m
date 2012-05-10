function [ mNs, mospa, rmse, MAP_rmse, corr_rmse ] = performance_measures( flags, params, pts, times, true_tau, true_w, true_intx )
%PERFORMANCE_MEASURES Calculates various performance measures for VRPF/S
%output

if isempty(pts)
    [mNs, mospa, rmse.pos, rmse.vel, rmse.pos_over_time, rmse.vel_over_time] = deal([]);
    [rmse, MAP_rmse, corr_rmse] = deal(rmse);
    return
end

sd = flags.space_dim;
K = params.K;

mospa = [];
rmse = [];
MAP_rmse = [];
corr_rmse = [];

% Count particles
Np = length(pts);

% Calculate mean number of jumps
mNs = mean(cat(1,pts.Ns));

% Find MAP particle
prob = zeros(Np,1);
for ii = 1:Np
    prob(ii) = sum(pts(ii).tau_prob) + sum(pts(ii).w_prob) + sum(pts(ii).lhood); % There should be P(x_0) term in this.
end
[~, MAP_ind] = max(prob);

% Calculate MAP estimate OSPA between particle jump sequences and real sequence
if ~isempty(true_tau)
    mospa = OSPA(true_tau, pts(MAP_ind).tau(1:pts(MAP_ind).Ns), 1, 2);
end

% Add drift to velocity and work out RMSE
if ~isempty(true_intx) && (flags.dyn_mod == 2)
    % Inferred
    pts_intx = cat(3, pts.intx);
    for ii = 1:Np
        cpi = 1;
        for kk = 1:K
            if (cpi<pts(ii).Ns)&&(times(kk)>pts(ii).tau(cpi+1))
                cpi = cpi + 1;
            end
            pts_intx(sd+1:2*sd, kk, ii) = pts_intx(sd+1:2*sd, kk, ii) + pts(ii).w(sd+1:2*sd,cpi);
        end
    end
    intx = mean(pts_intx, 3);
    
    % True
    mod_true_intx = true_intx;
    if ~isempty(true_w)
        cpi = 1;
        for kk = 1:K
            if (cpi<length(true_tau))&&(times(kk)>true_tau(cpi+1))
                cpi = cpi + 1;
            end
            mod_true_intx(sd+1:2*sd, kk) = mod_true_intx(sd+1:2*sd, kk) + true_w(sd+1:2*sd,cpi);
        end
    end
    
    error = abs(mod_true_intx - intx);
    corr_rmse.pos_over_time = sqrt( sum(error(1:sd,:).^2,1));
    corr_rmse.vel_over_time = sqrt( sum(error(sd+1:2*sd,:).^2,1));
    corr_rmse.pos = sqrt(mean(corr_rmse.pos_over_time.^2));
    corr_rmse.vel = sqrt(mean(corr_rmse.vel_over_time.^2));
end


% Calculate RMSEs for MMSE estimate
if ~isempty(true_intx)
    intx = mean(cat(3, pts.intx), 3);
    error = abs(true_intx - intx);
    rmse.pos_over_time = sqrt( sum(error(1:sd,:).^2,1));
    rmse.vel_over_time = sqrt( sum(error(sd+1:2*sd,:,:).^2,1));
    rmse.pos = sqrt(mean(rmse.pos_over_time.^2));
    rmse.vel = sqrt(mean(rmse.vel_over_time.^2));
end

% Calculate RMSEs for MAP estimate
if ~isempty(true_intx)
    intx = pts(MAP_ind).intx;
    error = abs(bsxfun(@minus, true_intx, intx));
    MAP_rmse.pos_over_time = sqrt(mean( sum(error(1:sd,:,:).^2,1), 3));
    MAP_rmse.vel_over_time = sqrt(mean( sum(error(sd+1:2*sd,:,:).^2,1), 3));
    MAP_rmse.pos = sqrt(mean(MAP_rmse.pos_over_time.^2));
    MAP_rmse.vel = sqrt(mean(MAP_rmse.vel_over_time.^2));
end

% figure, plot(times, rmse_over_time)

end

