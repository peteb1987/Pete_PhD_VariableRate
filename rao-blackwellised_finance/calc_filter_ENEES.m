function [ ENEES ] = calc_filter_ENEES( flags, params, true_intx, filt_pts )
%CALC_FILTER_ENEES Calculate ENEES for a set of filter particles

K = params.K;
Np = length(filt_pts);

ENEES = ones(1, K);

for kk = 2:K
    
    pts_intx = cat(3, filt_pts.intmu);
    
    % Errors
    errors = bsxfun(@minus, squeeze(pts_intx(:,kk,:)), true_intx(:,kk));
    
    % Calculate error in MMSE point estimate
    mmse_error = errors/Np;
    
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

