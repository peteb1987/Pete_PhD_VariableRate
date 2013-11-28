function [ time, cp_time, cp_param, cp_state, state, observ ] = tracking_generatedata( model )
%tracking_generatedata Generate a data set for the variable rate tracking
%model.

% This version uses a regular time grid.

% Times
time = (1:model.K)*model.T;

%%% Simulate changepoint times %%%

% First one
cp_time = 0;
cp_param = tracking_paramtrans(model, []);
cp_state = model.x0;

% Loop through the rest
while cp_time(end) < time(end)
    
    % Sample next changepoint
    next_cp_time = cp_time(end) + tracking_periodtrans(model, cp_param(:,end));
    next_cp_param = tracking_paramtrans(model, cp_param(:,end));
    next_cp_state = tracking_evaluatestate(model, cp_time(end), cp_param(:,end), cp_state(:,end), next_cp_time);
    
    % Append it
    cp_time = [cp_time, next_cp_time];
    cp_param = [cp_param, next_cp_param];
    cp_state = [cp_state, next_cp_state];
    
end

% Remove the last one that overshot
cp_time(end) = [];
cp_param(:,end) = [];
cp_state(:,end) = [];

% Initialise state and observation arrays
state = zeros(model.ds, model.K);
observ = zeros(model.do, model.K);

% Loop through observation times
for kk = 1:model.K
    
    % Find the most recent changepoint in the list
    cp_ind = most_recent_changepoint(cp_time, time(kk));
    
    % Evaluate the state
    state(:,kk) = tracking_evaluatestate(model, cp_time(cp_ind), cp_param(:,cp_ind), cp_state(:,cp_ind), time(kk));
    
    % Simulate an observation
    observ(:,kk) = tracking_observation(model, state(:,kk));
    
end

end
