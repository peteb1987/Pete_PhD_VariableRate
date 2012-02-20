function [ MAP_pt ] = pick_max_particle( params, pts, times, observs )
%PICK_MAX_PARTICLE Find the MAP particle

F = params.F; C = params.C; L = params.L;

Np = size(pts, 1);

prob = zeros(Np,1);

last_tau = 0;

for ii = 1:Np
    
    % First add up chamgepoint sequence prior
    for k = 1:length(pts(ii).tau)
        
        tau = pts(ii).tau(k);
        prob(ii) = prob(ii) + log( exppdf(tau-last_tau, params.x_jump_rate+params.xdot_jump_rate) );
        
    end
    
    prob(ii) = prob(ii) + log(1 - expcdf(params.T-last_tau, params.x_jump_rate+params.xdot_jump_rate) );
    
    % Next add up predictive likelihoods
    m = [observs(1); 0];
    P = [params.x_start_sd^2, 0; 0 params.xdot_start_sd^2];
    ji = 2;
    for k = 2:params.K
        
        t = times(k);
        
        % See if a jump happened
        [A, Q] = lti_disc(F, L, C,t-times(k-1));
        if (ji<=length(pts(ii).tau))&&(t>pts(ii).tau(ji))
            if pts(ii).type(ji)==1
                Q = Q + [params.x_jump_sd^2, 0; 0, 0];
            elseif pts(ii).type(ji)==2
                Q = Q + [0, 0; 0, params.xdot_jump_sd^2];
            end
            ji = ji + 1;
        end
        
        [m, P] = kf_predict(m,P,A,Q);
        [m, P, ~, ~, ~, lhood] = kf_update(m, P, observs(k), [1, 0], params.R);
        
        prob(ii) = prob(ii) + log(lhood);
        
    end
    
end

MAP_ind = find(prob==max(prob),1);

MAP_pt = pts(MAP_ind);

end

