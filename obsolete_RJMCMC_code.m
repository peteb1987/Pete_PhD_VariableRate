if (last_tau+params.dt<t)&&(rand<0.5)
    
    % Add a jump
    
    % First find the reverse probabilities
    start_idx = find(min(times(times>last_tau))==times);
    [~, ~, rev_ppsl_prob] = tracking_acceleration_proposal(flags, params, last_x, last_tau, last_w, times(1:kk), observ(:,1:kk));
    old_accel_prob = log(mvnpdf(last_w', zeros(1,params.rnd_dim), params.Q));
    old_trans_prob = tracking_calc_jump_trans_prob( params, last_tau );
    old_lhood = sum(pts_lhood(jj,start_idx:end));
    
    % Now add a jump
    extra_tau = unifrnd(last_tau+params.dt, t);
    stop_idx = find(max(times(times<extra_tau))==times);
    [w1, ppsl_prob_part1, ~] = tracking_acceleration_proposal(flags, params, last_x, last_tau, last_w, times(1:stop_idx), observ(:,1:stop_idx));
    extra_x = tracking_calc_next_state(flags, last_x, extra_tau-last_tau, w1);
    [w2, ppsl_prob_part2, ~] = tracking_acceleration_proposal(flags, params, extra_x, extra_tau, last_w, times(1:kk), observ(:,1:kk));
    ppsl_prob = ppsl_prob_part1 + ppsl_prob_part2 + log(unifpdf(extra_tau, last_tau, t));
    [new_lhood_part1, ppsl_intx_part1] = tracking_calc_likelihood(flags, params, last_x, last_tau, w1, times(1:stop_idx), observ(:,1:stop_idx));
    [new_lhood_part2, ppsl_intx_part2] = tracking_calc_likelihood(flags, params, extra_x, extra_tau, w2, times(1:kk), observ(:,1:kk));
    ppsl_intx = ppsl_intx_part2; ppsl_intx(1,start_idx:stop_idx,:) = ppsl_intx_part1(1,start_idx:stop_idx,:);
    new_lhood = sum(new_lhood_part1) + sum(new_lhood_part2);
    new_trans_prob = tracking_calc_jump_trans_prob( params, last_tau, extra_tau ) + tracking_calc_jump_trans_prob( params, extra_tau );
    new_accel_prob = log(mvnpdf(w1', zeros(1,params.rnd_dim), params.Q)) + log(mvnpdf(w2', zeros(1,params.rnd_dim), params.Q));
    
    % MH acceptance
    acc_prob = (sum(new_lhood)+new_accel_prob+new_trans_prob)-(sum(old_lhood)+old_accel_prob+old_trans_prob)+(rev_ppsl_prob-ppsl_prob);
    if log(rand)<acc_prob
        pts_Ns(jj) = pts_Ns(jj)+1;
        pts_w(jj,pts_Ns(jj)-1,:) = w1;
        pts_w(jj,pts_Ns(jj),:) = w2;
        pts_x(jj,pts_Ns(jj),:) = extra_x;
        pts_tau(jj,pts_Ns(jj),:) = extra_tau;
        pts_intx(jj,start_idx:kk,:) = ppsl_intx(1,start_idx:kk,:);
        pts_lhood(jj,start_idx:kk,:) = new_lhood;
        last_x = extra_x;
        last_tau = extra_tau;
        last_w = w2;
        MH_type_counts(3) = MH_type_counts(3)+1;
    end
    
elseif (pts_Ns(jj)>1)
    
    % Remove the last jump
    
    % First find the probabilities for the backwards move
    tau_penult = pts_tau(jj,pts_Ns(jj)-1);
    x_penult = squeeze(pts_x(jj,pts_Ns(jj)-1,:));
    w_penult = squeeze(pts_w(jj,pts_Ns(jj)-1,:));
    start_idx = find(min(times(times>tau_penult))==times);
    stop_idx = find(max(times(times<last_tau))==times);
    [~, ~, rev_ppsl_prob_part1] = tracking_acceleration_proposal(flags, params, x_penult, tau_penult, w_penult, times(1:stop_idx), observ(:,1:stop_idx));
    [~, ~, rev_ppsl_prob_part2] = tracking_acceleration_proposal(flags, params, last_x, last_tau, last_w, times(1:kk), observ(:,1:kk));
    rev_ppsl_prob = rev_ppsl_prob_part1 + rev_ppsl_prob_part2 + log(unifrnd(tau_penult,t));
    old_lhood = sum(pts_lhood(jj,start_idx:end));
    old_accel_prob = log(mvnpdf(w_penult', zeros(1,params.rnd_dim), params.Q)) + log(mvnpdf(last_w', zeros(1,params.rnd_dim), params.Q));
    [ old_trans_prob ] = tracking_calc_jump_trans_prob( params, tau_penult, last_tau ) + tracking_calc_jump_trans_prob( params, last_tau );
    
    
    % Now replace it
    [w_replace, ppsl_prob, ~] = tracking_acceleration_proposal(flags, params, x_penult, tau_penult, w_penult, times(1:stop_idx), observ(:,1:stop_idx));
    [new_lhood, ppsl_intx] = tracking_calc_likelihood(flags, params, x_penult, tau_penult, w_replace, times(1:kk), observ(:,1:kk));
    new_accel_prob = log(mvnpdf(w_replace', zeros(1,params.rnd_dim), params.Q));
    [ new_trans_prob ] = tracking_calc_jump_trans_prob( params, tau_penult);
    
    % MH acceptance
    acc_prob = (sum(new_lhood)+new_accel_prob+new_trans_prob)-(sum(old_lhood)+old_accel_prob+old_trans_prob)+(rev_ppsl_prob-ppsl_prob);
    if log(rand)<acc_prob
        pts_Ns(jj) = pts_Ns(jj)-1;
        pts_w(jj,pts_Ns(jj)+1,:) = zeros(dr,1);
        pts_x(jj,pts_Ns(jj)+1,:) = zeros(ds,1);
        pts_tau(jj,pts_Ns(jj)+1,:) = 0;
        pts_w(jj,pts_Ns(jj),:) = w_replace;
        pts_intx(jj,start_idx:kk,:) = ppsl_intx(1,start_idx:kk,:);
        pts_lhood(jj,start_idx:kk,:) = new_lhood;
        MH_type_counts(4) = MH_type_counts(4)+1;
    end
    
end