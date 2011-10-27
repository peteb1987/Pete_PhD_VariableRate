function [ prob ] = exppdf( t, mu )
%EXPPDF Returns the pdf value for the exponential distribution. mu is the
% mean.

prob = exp(-t/mu)/mu;

end

