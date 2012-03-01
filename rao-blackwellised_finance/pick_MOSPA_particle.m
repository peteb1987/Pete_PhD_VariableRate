function [ MOSPA_pt ] = pick_MOSPA_particle( params, pts )
%MINIMUM_OSPA_SEQUENCE Find the particle which gives the minimum OSPA

Np = length(pts);
ospas = zeros(Np);

for ii = 1:Np
    
    fprintf(1, '%u. ', ii);
    
    for jj = 1:ii
        
        ospas(ii,jj) = OSPA(pts(ii).tau, pts(jj).tau, 1, 0.01);
        ospas(jj,ii) = ospas(ii,jj);
        
    end
end

sum_ospas = sum(ospas);

MOSPA_idx = find(sum_ospas==min(sum_ospas), 1);
MOSPA_pt = pts(MOSPA_idx);

end

