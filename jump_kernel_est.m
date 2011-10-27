function [t_grid, jump_kdest] = jump_kernel_est(T, pts_tau)
%JUMP_KERNEL_EST Make a kernel density estimate of the jump times from the
% particle output

Np = size(pts_tau,1);

% Make a finely spaced time grid
t_grid = linspace(0,T,1001)';

% Set the standard deviation equal to the spacing
sd = 2*(t_grid(2)-t_grid(1));

% Make a vector list of all the jump times
jump_list = pts_tau(:);
jump_list(jump_list==0)=[];

% Loop through and make the kd estimate
jump_kdest = zeros(size(t_grid));
for tt = 1:length(jump_list)
    jump_kdest = jump_kdest + exp(-0.5*(t_grid-jump_list(tt)).^2/(sd^2))/Np;
%     jump_kdest = jump_kdest + normpdf(t_grid,jump_list(tt),sd)/Np;
end

end

