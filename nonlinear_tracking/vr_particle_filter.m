function [ pf, pf_diagnostics ] = vr_particle_filter( display, algo, model, time, observ )
%VR_PARTICLE_FILTER Run block variable rate particle filter for
%conditionally deterministic tracking.

% Local copies
N = algo.N;
B = algo.B;
L = algo.L;

% Initialise particles
for ii = 1:N
    pf(ii).cp_time = 0;
    pf(ii).cp_state = model.x0;
    pf(ii).cp_param = tracking_paramtrans(model, []);
    pf(ii).cp_pqr = 0;
    pf_weight = 0;
end

% Initialise loop variables
stop_time = 0;
kk = -B;

% Main loop
while 1
    
    % Times and indexes
    kk = kk + B;
    start_time = stop_time;
    stop_time = time(kk+L);
    
    % Store the old particle set
    last_pf = pf;
    last_pf_weight = pf_weight;
    
    % Resample
    ancestor = sample_weights(pf_weight, N);
    pf = pf(ancestor);
    
    % Loop through particles
    for ii = 1:N
        
        % Get the preceding changepoint
        cp_ind = most_recent_changepoint(pf(ii).cp_time, start_time);
        pre_cp_time = pf(ii).cp_time(cp_ind);
        pre_cp_param = pf(ii).cp_param(:,cp_ind);
        pre_cp_state = pf(ii).cp_state(:,cp_ind);
        
        % Propose an intermediate sequence
        [int_cp_time, int_cp_param, int_cp_state] = tracking_sequencetrans( model, start_time, stop_time, pre_cp_time, pre_cp_param, pre_cp_state);
        
        % Propose better parameters values, and calculate pqr
        [cp_time, cp_param, cp_state, cp_pqr] = tracking_sequenceimprovement(model, start_time, stop_time, int_cp_time, int_cp_param, int_cp_state, time(kk+1:kk+L), observ(:,kk+1:kk+L));
        
        % Calculate probabilities
        
        % Calculate weight
        
    end
    
end


end

