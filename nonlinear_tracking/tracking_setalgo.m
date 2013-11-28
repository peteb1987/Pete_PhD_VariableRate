function algo = tracking_setalgo(test)
% Algorithm parameters

% Number of particles
algo.N = 1000;

% Smoothing trajectories
algo.S = 100;

% MH chain length / max rejections
algo.M = 10;

% Block and window lengths
algo.B = 5;
algo.L = 10;

end