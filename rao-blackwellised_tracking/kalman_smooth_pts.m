function [ pts ] = kalman_smooth_pts( flags, params, times, pts )
%KALMAN_SMOOTH_PTS Run a RTS smoother over each particle

ds = params.state_dim;
F = params.F; C = params.C; L = params.L;

K =length(times);
Np = length(pts);

for ii = 1:Np
    
    % Initialise arrays
    A_arr = zeros(ds,ds,K);
    Q_arr = zeros(ds,ds,K);
    last_t = 0;
    
    ji = 2;
    for k = 2:K
        
        t = times(k);
        
        % Create transition matrices
        [A, Q] = lti_disc(F, L, C, t-last_t);
        
        % See if a jump happened
        if (ji<length(pts(ii).tau))&&(t>pts(ii).tau(ji))
            if pts(ii).type(ji)==1
                Q = Q + [params.x_jump_sd^2, 0; 0, 0];
            elseif pts(ii).type(ji)==2
                Q = Q + [0, 0; 0, params.xdot_jump_sd^2];
            end
            ji = ji + 1;
        end
        
        % Store
        A_arr(:,:,k-1) = A;
        Q_arr(:,:,k-1) = Q;
        
        last_t = t;
        
    end
    
    % RTS smooth
    [intmu, intP] = rts_smooth(pts(ii).intmu, pts(ii).intP, A_arr, Q_arr);
    
    % Store
    pts(ii).intmu = intmu;
    pts(ii).intP = intP;
    
end

end

