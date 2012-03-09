function [ w, prob ] = sample_acceleration( flags, params, w )
%SAMPLE_ACCELERATION Sample accelerations from prior

if nargin==2
    w = mvnrnd(zeros(1, params.rnd_dim), params.Q)';
end
prob = log_mvnpdf(w', zeros(1, params.rnd_dim), params.Q);

end

