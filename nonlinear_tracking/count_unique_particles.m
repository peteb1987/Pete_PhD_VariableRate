function [ Nup, Nut, Nujt ] = count_unique_particles( times, pts )
%COUNT_UNIQUE_PARTICLES Does what it says on the tin

K = length(times);
Np = length(pts);
Nup = Np*ones(1,K);

% Loop through time
for k = 1:K
    
    if rem(k, 100)==0
        disp(['Time step ' num2str(k)]);
    end
    
    uniq = true(Np,1);
    
    % Loop through particles
    for ii = 1:Np
        
        tau_ii = pts(ii).tau(pts(ii).tau<times(k));
        
        % Loop through later particles
        for jj = ii+1:Np
            
            % Only look at particles which we still know to be unique
            if uniq(jj)
                
                tau_jj = pts(jj).tau(pts(jj).tau<times(k));
                
                % Compare number of jumps
                if length(tau_ii)==length(tau_jj)
                    
                    % Compare jump sequences
                    if all(tau_ii==tau_jj)
                        
                        % Compare start point
                        if (pts(ii).x(:,1)==pts(jj).x(:,1))
                            
                            % Particle is not unique, cross it off
                            uniq(jj) = false;
                            Nup(k) = Nup(k) - 1;
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

Nut = Nup(end);

tau = [];
% Loop through particles
for ii = 1:Np
    
    tau = [tau pts(ii).tau];
    
end

Nujt = length(unique(tau));

end

