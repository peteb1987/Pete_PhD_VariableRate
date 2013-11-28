function [ ENEES ] = calc_smoother_ENEES( flags, params, true_intx, pts )
%CALC_SMOOTHER_ENEES Calculate ENEES for a set of smoother particles

K = params.K;

ENEES = zeros(1, K);

pts_intx = cat(3, pts.intmu);
Np = length(pts);

for kk = 1:K
    
    % Errors
    errors = bsxfun(@minus, squeeze(pts_intx(:,kk,:)), true_intx(:,kk));
    
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

