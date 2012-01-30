function [ MAP_pt ] = pick_max_particle( params, pts )
%PICK_MAX_PARTICLE Find the MAP particle

Np = size(pts, 1);

prob = zeros(Np,1);

last_tau = 0;

for ii = 1:Np
    
    for k = 1:length(pts(ii).tau)
        
        tau = pts(ii).tau(k);
        prob(ii) = prob(ii) + log( exppdf(tau-last_tau, params.x_jump_rate+params.xdot_jump_rate) );
        
    end
    
    prob(ii) = prob(ii) + log(1 - expcdf(params.T-last_tau, params.x_jump_rate+params.xdot_jump_rate) );
    
end

MAP_ind = find(prob==max(prob),1);

MAP_pt = pts(MAP_ind);

end

