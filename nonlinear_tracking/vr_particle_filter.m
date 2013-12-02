function [ pf, pf_diagnostics ] = vr_particle_filter( display, algo, model, time, observ )
%VR_PARTICLE_FILTER Run block variable rate particle filter for
%conditionally deterministic tracking.

pf_diagnostics = [];

% Local copies
N = algo.N;
B = algo.B;
L = algo.L;
K = model.K;

fprintf(1, 'Initialising particles.\n');

% First frame
for ii = 1:N
    
    Lpre = L-B;
    start_time = 0;
    stop_time = time(Lpre);
    
    % First changepoint
    cp_time = 0;
    cp_param = tracking_paramtrans(model, []);
    cp_state = model.x0;
    
    % Sample first window from prior
    [int_cp_time, int_cp_param, int_cp_state] = tracking_sequencetrans( model, 0, time(1:Lpre), cp_time, cp_param, cp_state);
    int_cp_time = [cp_time, int_cp_time];
    int_cp_param = [cp_param, int_cp_param];
    int_cp_state = [cp_state, int_cp_state];
    
    % Propose better parameters values, and calculate pqr
    [cp_time, cp_param, cp_state, cp_pqr] = ...
                tracking_sequenceimprovement(...
                    model, start_time, stop_time, ...
                    [], [], [], ...
                    int_cp_time, int_cp_param, int_cp_state, time(1:Lpre), observ(:,1:Lpre));
    
    % Merge sequences
    pf(ii).cp_time  = cp_time;
    pf(ii).cp_param = cp_param;
    pf(ii).cp_state = cp_state;
    pf(ii).cp_pqr =   cp_pqr;
    
    % Likelihood
    lhood = tracking_likelihood(model, pf(ii).cp_time, pf(ii).cp_param, pf(ii).cp_state, time(1:Lpre), observ(:,1:Lpre));
    pf(ii).lhood = zeros(1,K);
    pf(ii).lhood(1:Lpre) = lhood;
    
    % Weight
    pf_weight(ii) = sum(lhood);
    
end
fprintf(1, '     ESS is %f.\n', calc_ESS(pf_weight));

% Initialise loop variables
stop_time = 0;
kk = -B;
Lt = L;
mm = 0;

% Main loop
while 1
    
    % Store the old particle set
    last_pf = pf;
    last_pf_weight = pf_weight;
    last_Lt = Lt;
    
    % Times and indexes
    kk = kk + B;
    if kk > K, break; end
    Lt = min(L,K-kk);
    if kk > 0
        start_time = time(kk);
    else
        start_time = 0;
    end
    stop_time = time(min(K,kk+Lt));
    
    % 
    mm = mm + 1;
    tic
    fprintf(1, 'Beginning processing step %u: indexes %u to %u of %u.\n', mm, kk+1, kk+Lt, K);
    
    % Resample
    ancestor = sample_weights(pf_weight, N, 2);
    pf = pf(ancestor);
    
    % Loop through particles
    for ii = 1:N
        
        % Get the preceding changepoint
        cp_ind = most_recent_changepoint(pf(ii).cp_time, start_time);
        pre_cp_time = pf(ii).cp_time(cp_ind);
        pre_cp_param = pf(ii).cp_param(:,cp_ind);
        pre_cp_state = pf(ii).cp_state(:,cp_ind);
                
%         % Propose an intermediate sequence
%         [int_cp_time, int_cp_param, int_cp_state] = tracking_sequencetrans( model, start_time, stop_time, pre_cp_time, pre_cp_param, pre_cp_state);
        
        % Find the first changepoint preceding the new stretch of time
        ext_cp_ind = most_recent_changepoint(pf(ii).cp_time, time(kk-B+last_Lt));
        pre_ext_cp_time = pf(ii).cp_time(cp_ind);
        pre_ext_cp_param = pf(ii).cp_param(:,cp_ind);
        pre_ext_cp_state = pf(ii).cp_state(:,cp_ind);
        
        % Build an intermediate sequence
        [ext_cp_time, ext_cp_param, ext_cp_state] = tracking_sequencetrans( model, time(kk-B+last_Lt), stop_time, pre_ext_cp_time, pre_ext_cp_param, pre_ext_cp_state);
        int_cp_time  = [pf(ii).cp_time(cp_ind+1:ext_cp_ind)    ext_cp_time];
        int_cp_param = [pf(ii).cp_param(:,cp_ind+1:ext_cp_ind) ext_cp_param];
        int_cp_state = [pf(ii).cp_state(:,cp_ind+1:ext_cp_ind) ext_cp_state];
        
        % Propose better parameters values, and calculate pqr
        [cp_time, cp_param, cp_state, cp_pqr] = ...
            tracking_sequenceimprovement(...
                            model, start_time, stop_time, ...
                            pre_cp_time, pre_cp_param, pre_cp_state, ...
                            int_cp_time, int_cp_param, int_cp_state, time(kk+1:kk+Lt), observ(:,kk+1:kk+Lt));
        
        % Merge sequences
        new_pt.cp_time  = [pf(ii).cp_time(1:cp_ind)    cp_time];
        new_pt.cp_param = [pf(ii).cp_param(:,1:cp_ind) cp_param];
        new_pt.cp_state = [pf(ii).cp_state(:,1:cp_ind) cp_state];
        new_pt.cp_pqr =   [pf(ii).cp_pqr(1:cp_ind)     cp_pqr];
        
        % Calculate probabilities
        lhood = tracking_likelihood(model, new_pt.cp_time, new_pt.cp_param, new_pt.cp_state, time(kk+1:kk+Lt), observ(:,kk+1:kk+Lt));
        new_pt.lhood = zeros(1,K);
        new_pt.lhood(1:kk) = pf(ii).lhood(1:kk);
        new_pt.lhood(kk+1:kk+Lt) = lhood;
        
        % Calculate weight
        %%% THERE SHOULD BE A SURVIVOR FUNCTION RATIO IN HERE %%%
        pf_weight(ii) = sum(new_pt.lhood(kk+1:kk+Lt)) - sum(pf(ii).lhood(kk+1:kk-B+last_Lt)) ...
                        + sum(new_pt.cp_pqr(cp_ind+1:end)) - sum(pf(ii).cp_pqr(cp_ind+1:end));
        
        % Assert that cp_time is increasing
        assert(all(diff(new_pt.cp_time)>0));
        
        % Replace the old particle
        pf(ii) = new_pt;
        
    end
    
    fprintf(1, '     Frame complete.\n');
    fprintf(1, '     ESS is %f.\n', calc_ESS(pf_weight));
    fprintf(1, '     Time taken was %fs.\n', toc);
    
end


end

