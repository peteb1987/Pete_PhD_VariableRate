function [ Nchild, parent ] = systematic_resample( weights, Np_out )
%SYSTEMATIC_RESAMPLE Calculates number of children for each particle
%according to a systematic resampling scheme, and generates an array of
%indexes to the parent of each child particle.

% weights are assumed to be linear

% Number of particles in
Np_in = length(weights);


% Generate random index array
u = (1/Np_out)*ones(Np_out, 1);
u(1) = rand/Np_out;
u = cumsum(u);

% Generate cumulative weight array
w_sum = cumsum(weights);

% Create array of parent indexes
parent = zeros(Np_out, 1);

% Enumerate offspring
Nchild = zeros(Np_in,1);
Nchild(1) = sum(u < w_sum(1));
parent(1:Nchild(1)) = 1;
cnt = Nchild(1);
for ii = 2:Np_in
    Nchild(ii) = sum((u < w_sum(ii))&(u > w_sum(ii-1)));
    parent(cnt+1:cnt+Nchild(ii)) = ii;
    cnt = cnt + Nchild(ii);
end

end

