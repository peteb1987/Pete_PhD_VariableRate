function [ params, times, observs, true_x ] = load_benchmark( flags, params )
%LOAD_BENCHMARK Load in an aeroplane benchmark trajectory

ds = params.state_dim; do = params.obs_dim; dr = params.rnd_dim;

load('TARGET6_legacy.mat')

if flags.space_dim == 2
    true_x = xyz_targ([1,4,2,5],:);
elseif flags.space_dim == 3
    true_x = xyz_targ([1,4,7,2,5,8],:);
end
times = xyz_targ(10,:)';

% Downsample
factor = 20;

times = times(factor:factor:end);
params.K = length(times);
params.dt = times(1);
params.T = times(end);
true_x = true_x(:,factor:factor:end);

% Start point
params.start_state = true_x(:,1);

% Generate some observations
observs = zeros(do, params.K);
for k = 1:params.K
    
    % Sample observsation
    mu = observation_mean(flags, params, true_x(:,k), zeros(dr,1));
    observs(:,k) = mvnrnd(mu', params.R)';
    
end

end

