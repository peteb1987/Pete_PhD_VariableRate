function [w, ppsl_prob, rev_ppsl_prob] = acceleration_proposal(flags, params, tau, x, w, times, observs)
%ACCELERATION_PROPOSAL Propose a change in the latest random
%variables (accelerations) in a variable rate track using a UKF
%approximation to the optimal importance distribution

% The last state is specified by x, tau and w. This function proposes a
% replacement for w, using the oberservations in observs at times. The
% forward and reverse proposal probabilities are also calculated.

% The proposal mechanism uses an unscented transform scheme, akin to a UKF
% but without the prediction step.

K = length(times);

% Set local variables
ds = params.state_dim;
dr = params.rnd_dim;
do = params.obs_dim;
Q = params.Q;
R = params.R;

% Find sigma point weights
[WM,WC,c] = ut_weights(dr+do, 0.5, 2, dr+do+1);
N_sigs = length(WM);

% Set initial mean and covariance
w_mn = zeros(dr, 1);
w_var = Q;

% Work out where to start
start_idx = find(min(times(times>tau))==times);

% Work out which frames to use
total_num_frames = K-start_idx+1;
used_num_frames = min(params.min_num_ppsl_frames,total_num_frames);
used_num_frames = used_num_frames + round((total_num_frames - used_num_frames)*params.prop_ppsl_frames);
frame_list = round( (1:used_num_frames)*total_num_frames/used_num_frames ) + start_idx-1;

% Loop through time
for k = frame_list
    
    % Calculate sigma points
    W = ut_sigmas([w_mn; zeros(do,1)],[w_var zeros(dr,do); zeros(do,dr) R],c);
    
    % Propagate SPs through transition and observation functions
    u = next_state(flags, params, x, W(1:dr,:), times(k)-tau);
    Y = observation_mean(flags, params,u,W);

    % Add observation noise
    Y = Y + W(dr+1:end,:);
    
    % Calculate observation mean, covariance and cross-covariance
    S  = zeros(do);
    C  = zeros(dr,do);
    mu = Y*WM;
    for ii=1:N_sigs
      S = S + WC(ii) * (Y(:,ii) - mu) * (Y(:,ii) - mu)';
      C = C + WC(ii) * (W(1:dr,ii) - w_mn) * (Y(:,ii) - mu)';
    end
    
    % Update recursions
    if times(k)~=tau
        K = C / S;
        w_mn = w_mn + K * (bsxfun(@minus, observs(:,k), mu));
        w_var = w_var - K * S * K';
    end
    
end

% Make sure its exactly symmetric (tolerances on mvnrnd are pretty tight)
w_var = (w_var+w_var')/2;

% Calculate probability of reverse move
rev_ppsl_prob = log_mvnpdf_fast_batch(w,w_mn,w_var);

% Propose a value for w
w = mvnrnd(w_mn', w_var)';
ppsl_prob = log_mvnpdf_fast_batch(w,w_mn,w_var);

end