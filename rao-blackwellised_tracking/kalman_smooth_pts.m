function [ pts ] = kalman_smooth_pts( flags, params, times, pts )
%KALMAN_SMOOTH_PTS Run a RTS smoother over each particle

ds = params.state_dim;

K =length(times);
Np = length(pts);

for ii = 1:Np
    
    % Initialise arrays
    A_arr = zeros(ds,ds,K);
    Q_arr = zeros(ds,ds,K);
    
    ji = 1;
    for k = 2:K
        
        t = times(k);
        interm_t = times(k-1);
        
        if (ji < size(pts(ii).cp.tau,2)) && (t > pts(ii).cp.tau(ji+1))
        
            [A, Q, ~] = construct_transmats(pts(ii).cp.tau(ji+1)-interm_t, pts(ii).cp.m(ji), pts(ii).cp.u(ji), params.proc_var);
            interm_t = pts(ii).cp.tau(ji);
            
            ji = ji + 1;
            [~, ~, Ajump] = construct_transmats(0, pts(ii).cp.m(ji), pts(ii).cp.u(ji), params.proc_var);
            A = Ajump*A;
            Q = Ajump*Q*Ajump';
           
        else
            
            A = eye(6);
            Q = zeros(6);
            
        end
        
        [Ainc, Qinc, ~] = construct_transmats(t-interm_t, pts(ii).cp.m(ji), pts(ii).cp.u(ji), params.proc_var);
        A = Ainc*A;
        Q = Ainc*Q*Ainc' + Qinc;
        
        % Store
        A_arr(:,:,k-1) = A;
        Q_arr(:,:,k-1) = Q;
        
    end
    
    % RTS smooth
    [mu, P] = rts_smooth(pts(ii).mu, pts(ii).P, A_arr, Q_arr);
    
    % Store
    pts(ii).mu = mu;
    pts(ii).P = P;
    
end

end

