function [ lhood ] = tracking_likelihood( model, cp_time, cp_param, cp_state, time, observ )
%TRACKING_LIKELIHOOD Calculate the likelihood of a set of observations

lhood = zeros(size(time));

for kk = 1:length(time)
    
    cp_ind = most_recent_changepoint(cp_time, time(kk));
    
    x = tracking_evaluatestate(model, cp_time(cp_ind), cp_param(:,cp_ind), cp_state(:,cp_ind), time(kk));
    dy = observ(:,kk) - tracking_h(model, x);
    if dy(1) > pi
        dy(1) = dy(1) - 2*pi;
    elseif dy(1) < -pi
        dy(1) = dy(1) + 2*pi;
    end
    
    lhood(kk) = loggausspdf(dy, zeros(size(dy)), model.R);

end

