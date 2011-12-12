function [ c ] = cross( a, b )
%CROSS Fast implementation of the cross product.
% Input matrices a and b must be equal size with vectors as columns.
% Output given by c(:,ii) = a(:,ii) x b(:,ii)

% NO INPUT CHECKING. BE CAREFUL.

c = [a(2,:).*b(3,:)-a(3,:).*b(2,:);
     a(3,:).*b(1,:)-a(1,:).*b(3,:);
     a(1,:).*b(2,:)-a(2,:).*b(1,:)];

end

