function [ ENEES ] = calc_filter_ENEES( flags, params, times_array, true_tau, true_w, true_intx, filt_part_sets, filt_weight_sets )
%CALC_FILTER_ENEES Calculate ENEES for a set of filter particles

K = params.K;
sd = flags.space_dim;

ENEES = ones(1, K);

Np = length(filt_part_sets{end});

% Velocity correction
if (flags.dyn_mod == 2)
    mod_true_intx = true_intx;
    if ~isempty(true_w)
        cpi = 1;
        for kk = 1:K
            if (cpi<length(true_tau))&&(times_array(kk)>true_tau(cpi+1))
                cpi = cpi + 1;
            end
            mod_true_intx(sd+1:2*sd, kk) = mod_true_intx(sd+1:2*sd, kk) + true_w(sd+1:2*sd,cpi);
        end
    end
end

for kk = 2:K
    
    pts = filt_part_sets{kk};
    wts = filt_weight_sets{kk};
    
    pts_intx = cat(3, pts.intx);
    
    % Velocity correction
    if (flags.dyn_mod == 2)
        % Inferred
        pts_intx = cat(3, pts.intx);
        for ii = 1:Np
            cpi = find(pts(ii).tau==max(pts(ii).tau(pts(ii).tau<times_array(kk))));
            pts_intx(sd+1:2*sd, kk, ii) = pts_intx(sd+1:2*sd, kk, ii) + pts(ii).w(sd+1:2*sd,cpi);
        end
    end
    
    % Errors
    errors = bsxfun(@minus, squeeze(pts_intx(:,kk,:)), mod_true_intx(:,kk));
    
    % Calculate error in MMSE point estimate
    mmse_error = errors*exp(wts);
    
    % Calculate sample error covariance
    sample_error_cov = errors*diag(exp(wts))*errors';
    
    % Calculate ENEES
    if isposdef(sample_error_cov)
        ENEES(kk) = (mmse_error'/sample_error_cov)*mmse_error;
    else
        ENEES(kk) = 1;
    end
end

end

