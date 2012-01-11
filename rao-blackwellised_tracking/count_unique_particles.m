function [ Nup ] = count_unique_particles( times, pts )
%COUNT_UNIQUE_PARTICLES Does what it says on the tin

K = length(times);
Np = length(pts);
Nup = Np*ones(1,K);

% Loop through time
for k = 1:K
    
    if rem(k, 100)==0
        disp(['Time step ' num2str(k)]);
    end
    
    unique = true(Np,1);
    
    % Loop through particles
    for ii = 1:Np
        
        tau_ii = pts(ii).tau(pts(ii).tau<times(k));
        
        % Loop through later particles
        for jj = ii+1:Np
            
            % Only look at particles which we still know to be unique
            if unique(jj)
                
                tau_jj = pts(jj).tau(pts(jj).tau<times(k));
                
                % Compare number of jumps
                if length(tau_ii)==length(tau_jj)
                    
                    % Compare jump sequences
                    if all(tau_ii==tau_jj)
                        
                        % Particle is not unique, cross it off
                        unique(jj) = false;
                        Nup(k) = Nup(k) - 1;
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

end

