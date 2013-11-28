function [ ospa ] = OSPA( s1, s2, p, c )
%OSPA Calculate a 1D OPSA between the list of points in s1 those in s2 (which are row vectors).

min_bid = 1E-5;

if size(s1,2)>size(s2,2)
    temp = s1;
    s1 = s2;
    s2 = temp;
    clear temp;
end

% Number of points
M = size(s1,2);
N = size(s2,2);

% Initialise arrays
assoc = zeros(M,M+N);           % Matix of associations. assoc(i,j)=1 indicates truth(i) assigned to estimate(j)
prices = zeros(1, M+N);         % Price vector. prices(j) is the cost paid for the assignment to estimate(j)
payoffs = zeros(M,M+N);         % Payoff matrix. payoffs(i,j) is the positive effect of assigning truth(i) to estimate(j)
payoffs(1:M,N+1:N+M) = -inf(M);
payoffs(logical([zeros(M,N) eye(M)])) = -c^p;
for ii = 1:M;
    for jj = 1:N
        payoffs(ii,jj) = -dist(s1(ii), s2(jj))^p;
    end
end

% Auction
while ~all(sum(assoc, 2))
    
    % Loop over list 1
    for ii = 1:M
        
        if sum(assoc(ii, :))==0
            
            % Find the target or clutter ("object") that gives the biggest reward
            rewards = payoffs(ii, :) - prices;
            k = find(rewards==max(rewards));
            if length(k)>1
                k = k(payoffs(ii, k)==max(payoffs(ii, k)));
            end
            
            % Calculate the bidding increment
            best_reward = rewards(k);            
            second_best_reward = max(rewards([1:k-1, k+1:end]));
            bid = best_reward - second_best_reward;
            if bid < min_bid
                bid = min_bid;
            end
            
            if isnan(bid)
                error('Bid is NaN!');
            end
            
            % Increase the price and assign this target
            prices(k) = prices(k) + bid;
            assoc(:, k) = 0;
            assoc(ii, k) = 1;
            
        end
        
    end
    
end

ospa = (mean(-payoffs(logical(assoc))))^(1/p);

end

function d = dist(p1, p2)
d = abs(p1-p2);
end
