function [low_i, high_i] = beta_confidence_interval(n_i, N, low, high, alpha);
% Compute confidence intervals of the beta distribution.
%
% ARGUMENTS
%
% OPTIONAL ARGUMENTS
%  alpha - prior (defaults to 1/2 for Jeffreys prior, uniform otherwise)
% 
% RETURNS 


% Set prior to Jeffreys if not specified.
if (nargin < 5)
  alpha = 1/2; 
end

% Determine number of points to evaluate
npoints = length(n_i);
low_i = zeros(size(n_i));
high_i = zeros(size(n_i));
for i = 1:npoints
  n = n_i(i);
  low_i(i) = betainv(low, alpha+n, alpha+(N-n));
  high_i(i) = betainv(high, alpha+n, alpha+(N-n));
end
