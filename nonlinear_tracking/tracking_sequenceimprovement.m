function [cp_time, cp_param, cp_state, cp_pqr] = tracking_sequenceimprovement(model, start_time, stop_time, pre_cp_time, pre_cp_param, pre_cp_state, int_cp_time, int_cp_param, int_cp_state, time, observ)
%TRACKING_SEQUENCEIMPROVEMENT Optimise changepoint parameter and propose a
%new value. Calculate p/q ratio for each.

num_cp = length(int_cp_time);

cp_time = int_cp_time;
cp_state = int_cp_state;

if num_cp == 0
    % No changepoints in the window, so nothing to do. Hooray!
    cp_param = int_cp_param;
    cp_pqr = zeros(1,0);
    
    return
    
% elseif num_cp == 1
%     % One changepoint
%     
% %     cp_param = int_cp_param;
% %     cp_pqr = 0;
%     
%     % Optimise it
%     tol = 0.01;
%     h_of = @(u) neg_log_lhood(model, time, observ, pre_cp_time, pre_cp_param, pre_cp_state, int_cp_time, u, int_cp_state);
%     options = optimset('Display','notify-detailed', 'GradObj','on', 'TolX',tol);
% %     options = optimset('Display','notify-detailed', 'TolX',tol);
%     [opt_cp_param, ~, ~, ~, ~, hess] = fminunc(h_of, int_cp_param, options);
%     
%     % Propose a new one
%     ppsl_vr = 0.001*eye(2); inv(full(hess));
%     cp_param = mvnrnd(opt_cp_param', ppsl_vr)';
%     
%     % p/q ratio
%     q_prob = loggausspdf(cp_param, opt_cp_param, ppsl_vr);
%     p_prob = loggausspdf(cp_param, model.param_trans_mn, model.param_trans_vr);
%     cp_pqr = p_prob - q_prob;

else
    % Optimise changepoint parameters
    
%     cp_param = int_cp_param;
%     cp_pqr = zeros(size(cp_time));
    
    cp_param = zeros(size(int_cp_param));

    % Loop through
    for cc = 1:num_cp
        
        % Find the first and last observations for this segment
        first_ob_ind = earliest_observation(time, int_cp_time(cc));
        if cc < num_cp
            last_ob_ind = earliest_observation(time, int_cp_time(cc+1))-1;
        else
            last_ob_ind = length(time);
        end
        
        % Optimise parameters
        h_of = @(u) neg_log_lhood(model, time(first_ob_ind:last_ob_ind), observ(:,first_ob_ind:last_ob_ind), [], [], [], cp_time(cc), u, cp_state(:,cc));
        options = optimset('Display','notify-detailed', 'GradObj','on', 'TolX',0.01, 'TolFun',0.1);
        [opt_cp_param, ~, ~, ~, ~, hess] = fminunc(h_of, int_cp_param(:,cc), options);
        
        % Propose a new one
        ppsl_vr = 0.1*eye(2); inv(full(hess));
        cp_param(:,cc) = mvnrnd(opt_cp_param', ppsl_vr)';
        
        % p/q ratio
        q_prob = loggausspdf(cp_param(:,cc), opt_cp_param, ppsl_vr);
        p_prob = loggausspdf(cp_param(:,cc), model.param_trans_mn, model.param_trans_vr);
        cp_pqr(cc) = p_prob - q_prob;
        
        % Update next state
        if cc < num_cp
            cp_state(:,cc+1) = tracking_evaluatestate(model, cp_time(:,cc), cp_param(:,cc), cp_state(:,cc), cp_time(cc+1));
        end
        
    end
    
end

end

function [func, grad] = neg_log_lhood(model, time, observ, pre_cp_time, pre_cp_param, pre_cp_state, cp_time, cp_param, cp_state)

% Find first observation after changepoint
ob_ind = earliest_observation(time, cp_time);

func = 0;
grad = 0;
for kk = ob_ind:length(time)
    x = tracking_evaluatestate(model, cp_time, cp_param, cp_state, time(kk));
    du = cp_param-model.param_trans_mn;
    dy = observ(:,kk)-tracking_h(model, x);
    if dy(1) > pi
        dy(1) = dy(1) - 2*pi;
    elseif dy(1) < -pi
        dy(1) = dy(1) + 2*pi;
    end
    func = func + dy'*(model.R\dy) + du'*(model.param_trans_vr\du);
    
    if nargout > 1
        H = tracking_obsjacobian(model,x);
        J = tracking_statejacobian(model, cp_time, cp_param, cp_state, time(kk), x);
        grad = grad - 2*( dy'/model.R ) * H * J + 2*du'/model.param_trans_vr;
    end
    
end

grad = grad';

end