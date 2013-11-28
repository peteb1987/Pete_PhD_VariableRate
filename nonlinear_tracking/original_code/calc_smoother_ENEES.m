function [ ENEES ] = calc_smoother_ENEES( flags, params, times, true_tau, true_w, true_intx, pts )
%CALC_SMOOTHER_ENEES Calculate ENEES for a set of smoother particles

K = params.K;
sd = flags.space_dim;

ENEES = zeros(1, K);

pts_intx = cat(3, pts.intx);
Np = length(pts);

% Velocity correction
if (flags.dyn_mod == 2)
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
end

for kk = 1:K
    
    % Velocity correction
    if (flags.dyn_mod == 2)
        for ii = 1:Np
            cpi = 1;
            if (cpi<pts(ii).Ns)&&(times(kk)>pts(ii).tau(cpi+1))
                cpi = cpi + 1;
            end
            pts_intx(sd+1:2*sd, kk, ii) = pts_intx(sd+1:2*sd, kk, ii) + pts(ii).w(sd+1:2*sd,cpi);
        end
    end
    
    % Errors
    errors = bsxfun(@minus, squeeze(pts_intx(:,kk,:)), mod_true_intx(:,kk));
    
    % Calculate error in MMSE point estimate
    mmse_error = mean(errors, 2);
    
    % Calculate sample error covariance
    sample_error_cov = errors*errors'/Np;
    
    % Calculate ENEES
    if isposdef(sample_error_cov)
        ENEES(kk) = (mmse_error'/sample_error_cov)*mmse_error;
    else
        ENEES(kk) = 1;
    end
end

end

