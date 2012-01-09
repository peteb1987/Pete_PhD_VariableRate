function [ r ] = rande( lambda, n )
%RANDE Generate random numbers from an exponential distribution with mean
% lambda. (optional) n defines dimensions of array, just as for 'rand'.

if nargin==1
    n=1;
end

r = -log(rand(n))*lambda;

end

