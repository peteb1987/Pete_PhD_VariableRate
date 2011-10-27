function [trans_prob, prev_ji, next_ji] = rb_transition_probability(params, t, past_tau, future_tau, past_type, future_type)
%RB_TRANSITION_PROBABILITY Calculate transition probability of a future set
% of jump times and types given a past set of jump times and types

% Get next jump times
[next_1_tau, next_1_ji] = min(future_tau(((future_type==0)|(future_type==1))&(future_tau>t)));
[next_2_tau, next_2_ji] = min(future_tau(((future_type==0)|(future_type==2))&(future_tau>t)));
if ~isempty(next_1_tau), next_1_ji = find(next_1_tau==future_tau); end
if ~isempty(next_2_tau), next_2_ji = find(next_2_tau==future_tau); end

% Find previous jump times
temp = past_tau; temp(temp>=t)=0;
temp1 = temp; temp1(~((past_type==0)|(past_type==1))) = NaN;
temp2 = temp; temp2(~((past_type==0)|(past_type==2))) = NaN;

[prev_1_tau, prev_1_ji] = max(temp1,[],2);
[prev_2_tau, prev_2_ji] = max(temp2,[],2);

% Create an array of transition probabilities
trans_prob = zeros(size(future_tau,1),1);

% Type 1 jumps
if ~isempty(next_1_tau)
    trans_prob = trans_prob + log(exppdf(next_1_tau-prev_1_tau, 1/params.x_jump_rate));
else
    trans_prob = trans_prob + log(1-expcdf(params.T-prev_1_tau, 1/params.x_jump_rate));
end

% Type 2 jumps
if ~isempty(next_2_tau)
    trans_prob = trans_prob + log(exppdf(next_2_tau-prev_2_tau, 1/params.xdot_jump_rate));
else
    trans_prob = trans_prob + log(1-expcdf(params.T-prev_2_tau, 1/params.xdot_jump_rate));
end

% Find indexes of next and previous jumps for gluing
prev_ji = max(prev_1_ji,prev_2_ji);
next_ji = [];
if ~isempty(next_1_ji)
    next_ji = next_1_ji;
    if ~isempty(next_2_ji)
        next_ji = min(next_1_ji,next_2_ji);
    end
else
    if ~isempty(next_2_ji)
        next_ji = next_2_ji;
    end
end

end

