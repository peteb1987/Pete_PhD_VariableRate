function [ C ] = arraymatprod( A, B )
%ARRAYMATPROD Fast multiple matrix-vector multiplications

% The first input is a 3D array where each page is a matrix. Size MxNxP
% The second input is a 2D array where each column is a vector. Size NxP
% The output is given by C(:,:,ii) = A(:,:,ii) * B(:,ii)

% NO INPUT CHECKING. BE CAREFUL.

B = permute(B, [3 1 2]);
C = squeeze(sum( bsxfun(@times, A, B), 2));

end

