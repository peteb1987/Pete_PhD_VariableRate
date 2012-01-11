function [ filt_pts ] = create_filter_set( Np, filt_part_sets, filt_weight_sets )
%CREATE_FILTER_SET Take a cell array with the Kitigawa smoother output from
%each frame and create a single set of particles representing the filter
%output.

filt_pts = filt_part_sets{1};
K = length(filt_part_sets);

for k = 1:K
    
    % Resample to uniform weight and correct number of particles
    [~, parent] = systematic_resample(exp(filt_weight_sets{k}), Np);
    
    % Add to array
    for ii = 1:Np
        
        filt_pts(ii).intmu(:,k) = filt_part_sets{k}(parent(ii)).intmu(:,k);
        
    end
    
end

end

