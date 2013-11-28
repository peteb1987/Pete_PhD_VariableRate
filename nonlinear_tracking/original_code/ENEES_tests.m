bias = [3 0]';
Sigma = eye(2);

N = 100;


rpts = 1000;
MENEES = 0;
for ii = 1:rpts
    x = mvnrnd(bias', Sigma, N)';
    
    x_est = mean(x,2);
    cov_est = cov(x');% + x_est*x_est';
    
    if isposdef(cov_est)
        ENEES = (x_est'/cov_est)*x_est;
    else
        ENEES = 1;
    end
    
    MENEES = MENEES + ENEES;
end

MENEES = MENEES/rpts