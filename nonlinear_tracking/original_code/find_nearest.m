function [ idx ] = find_nearest( times, tau, higher )
%FIND_NEAREST Find the index of the element in times nearest to tau.
%Specify higher or lower

% If there is no appropriate element or no tau specified, return empty matrix

if isempty(tau)
    idx = [];
    return;
end

if higher
    el = min(times(times>tau));
else
    el = max(times(times<tau));
end

if ~isempty(el)
    idx = find(el==times);
else
    idx = [];
end

end

