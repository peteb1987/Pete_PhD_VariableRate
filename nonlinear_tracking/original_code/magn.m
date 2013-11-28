function [ m ] = magn( x )
%MAGN magnitude of a vector

m = sqrt( sum(x.^2, 1) );

end

