function [ prob ] = expcdf( t, mu )
%EXPPDF Returns the cdf value for the exponential distribution. mu is the
% mean.

prob = 1 - exp(-t/mu);

end

