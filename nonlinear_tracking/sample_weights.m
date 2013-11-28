function [ ancestor, selected_weight ] = sample_weights( weight, N, resam_type )
%SAMPLE_WEIGHTS Samples an array of N ancestor indexes from a set of
%weights. Also (optionally), the weights of the ancesors.

% Default resampling type
if (nargin<3)||isempty(resam_type)
    resam_type = 1;
end

% Row vectors only here, please
weight = weight(:)';

% Convert weights to linear domain and normalise
weight = weight - max(weight);
weight = exp(weight);
weight = weight/sum(weight);

% Catch NaNs
assert(~any(isnan(weight)));

% Create bin boundaries
edges = min([0 cumsum(weight)],1);
edges(end) = 1;

% Sample indexes
if resam_type == 1
    % Multinomial
    idx = rand(N,1);
elseif resam_type == 2
    % Systematic
    idx = (1/N)*ones(N, 1);
    idx(1) = rand/N;
    idx = cumsum(idx);
else
    error('That''s not a valid choice for resampling type');
end

% Draw particles
[~, ancestor] = histc(idx,edges);
ancestor = ancestor';

% Selected weights
if nargout > 1
    selected_weight = log(weight(ancestor));
end

end

