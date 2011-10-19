function [ part_sets ] = vr_inference( params, times, observs )
%VR_INFERENCE Run variable rate particle filter and smoother algorithms

% params is a structure with all the model and algorithm parameters
% times is a vector specifying the time grid
% observs is a vector of observations

% Create an array to store particle sets for each time frame
part_sets = cell(params.T, 1);

% Count the obervations
K = numel(times);
T = times(K);
assert(numel(observs)==K);



% Loop through observations
for k = 1:K
    
    
    
end









for tind=2:NumObservations
    
    tprev = observations(tind-1,1);
    t = observations(tind,1);

    NumParticles = numel(particles);
    NumParticlesNextGen = 0;
    NextGenParticles = particle.empty(1,0);
    wsum = 0;
    jumpingweight = 0;
    numjump    = 0;
    numnonjump = 0;
    Nis = [];
    
    for i = 1:NumParticles
    
        p = particles(i);
        
        N0 = 0;
        N1 = 0;
        
        Ni = max(1, floor(N*p.w(tind-1)));
         
        % if Ni==1 then only probabilistically propagate this particle
        if(Ni==1 && rand(1)<1-N*p.w(tind-1)), Ni=0; end
        
        Nis(i) = Ni;
        
        if(Ni>0) %if zero then don't sample this particle
            
            % For this particle, N0 is the number of offspring that do not
            % jump and N1 is the number of offspring that do
            
            newparticles = particle.empty(1,0);

            for j=1:Ni
                taunew = p.samplenewtau(tprev);
                hasjumped = taunew < t;
                if(max(hasjumped) > 0)
                    % In the (unlikely) event that both jump we approximate
                    % by saying that both jumps occur at the first jump
                    % time.  This is an approximation and could be improved
                    % (but probably doesn't matter much).
                    jumptime = min(taunew);
                    q = p.AddState(jumptime, hasjumped);
                    N1 = N1 + 1;
                    newparticles(N1) = q;  % Copy the particle to this group
                else
                    N0 = N0 + 1;
                end
            end
            
            numjump = numjump+N1;
            numnonjump = numnonjump+N0;

            Nitilde = N1 + 1;

            for j=1:N1
                newparticles(j).w(tind-1) = p.w(tind-1) / Nitilde;
            end

            p.w(tind-1) = N0*p.w(tind-1) / Ni;
            newparticles(Nitilde) = p;

            for j=1:Nitilde
                %Nwprev(j) = newparticles(j).w(tind-1);
                
                
                % Use the following line with particleA
                %[obsprob, newparticles(j)] = newparticles(j).observationprobability(observations(1:tind,:));
                
                % Use the following line with particle
                [obsprob, newparticles(j)] = newparticles(j).addobservation(observations(tind,:));
                
                newparticles(j).w(tind) = newparticles(j).w(tind-1)*obsprob;
                %Nw(j) = newparticles(j).w(tind);
                wsum = wsum + newparticles(j).w(tind);
                if(j<=N1)
                    jumpingweight = jumpingweight+newparticles(j).w(tind);
                end
            end
            
            % add these particles to the total set
            for j=1:Nitilde
                NextGenParticles(NumParticlesNextGen+j) = newparticles(j);
            end
            NumParticlesNextGen = NumParticlesNextGen + Nitilde;
        end
        
    end
    
    % Normalize weights
    for i=1:NumParticlesNextGen
        NextGenParticles(i).w(tind) = NextGenParticles(i).w(tind)/wsum;
    end
    
    particles = NextGenParticles;
    
   
    % Predict ahead
    % predict slightly behind next observation to prevent look-ahead using next
    % observation
    
    if(tind<NumObservations), tpredict = observations(tind+1, 1) - 0.00001;
    else tpredict = 2*observations(tind,1) - observations(tind-1,1);
    end
    
    %tpredict = observations(tind,1) + 2;
    %tpredict = t;
    
    %tind-1 because want predict to start in first row
    predict(tind-1,1) = t;
    predict(tind-1,2) = tpredict;
    
    predict(tind-1,3) = Predict(particles, tind, tpredict, params);
    
    %particle storage
    particlecollection = [particlecollection {particles}];
end

for i=1:NumParticlesNextGen
    NS(i) = particles(i).NumStates;
end



end

