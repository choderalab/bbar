function [f_mean, f_lower, f_upper] =  BBAR_test(WF, WR, CI, permsigns)
% ---
% Bayesian BAR
% ---
%
% Calculate the free energy difference between two states using the
% Bayesian extension of the Bennett acceptance ratio, described in
%
% [1] Maragakis P, Ritort F, Bustamante C, Karplus M, and Crooks GE.
% Bayesian estimates of free energies from nonequilibrium work data in
% the presence of instrument noise. JCP 129:024102, 2008.
%
% This variant was suggested by Gavin Crooks.
%
% USAGE
% 
% [f_lower, f_mean, f_upper] = BBAR(WF, WR, CI, prior)
%
% PARAMETERS
%
% WF:       The total work done in forward processes (in dimensionless units)
% WR:       The total work done in reverse processes (in dimensionless units)
% CI:       Confidence interval(s) of posterior to (0 < CI < 1). If CI is a vector, f_lower and f_upper will be vectors.
% permsigns (optional) - each row is a unique sign permutation with NF +1 and NR -1
%
% RETURNS
%
% f_mean:     The mean of the posterior (in dimensionless units)
% f_lower:    The lower free energy bound for the specified confidence interval (in dimensionless units)
% f_upper:    The upper free energy bound for the specified confidence interval (in dimensionless units)
%
% NOTE
%
% At least one forward and one reverse work value must be provided in order for the PDF to be normalizable
% and confidence bounds to be computed.
%
% This is currently very expensive to compute!  Use it with small numbers of samples (< 10 total) only!

% Determine number of forward and backward samples
NF = length(WF);
NR = length(WR);

% Compute total number of samples.
N = NF+NR;

% Generate list of sign permutations.
truesigns = [ones(1,NF) -ones(1,NR)];
if (nargin == 3)
  nperms = nchoosek(N, NF);
  disp(sprintf('Generating %d sign permutations...', nperms));
  tic; permsigns = uniqueperms(truesigns); toc
else
  nperms = size(permsigns,1); % get size from input arguments
end

% Concatenate rectified work values.
W = [WF; -WR];

% Determine minimum and maximum work values, bracketing the free energy difference estimate.
Wmin = min(W);
Wmax = max(W);
DW = Wmax - Wmin; % interval size

% Define the logarithm of the unnormalized probability density function.
% This takes a row vector of free energy differences 'df' as input.
% This should be equivalent to
% log [ exp[\sum_n sign_n (W_n - df)] / \permsum_p  exp[\sum_n permsign_pn (W_n - df)] ]
log_unnormalized_posterior = @(df) - log(sum(exp(0.5*permsigns*(W*ones(size(df))-ones(size(W))*df) - 0.5*repmat(truesigns,[nperms,1])*(W*ones(size(df))-ones(size(W))*df)),1));

% Determine the maximum log likelihood estimate starting from the center of this interval.
df_ML = fminsearch(@(x) -log_unnormalized_posterior(x), (Wmin+Wmax)/2);

% Redefine log unnormalized posterior so that maximum value is zero.
max_log = log_unnormalized_posterior(df_ML);
log_unnormalized_posterior = @(df) - log(sum(exp(0.5*permsigns*(W*ones(size(df))-ones(size(W))*df) - 0.5*repmat(truesigns,[nperms,1])*(W*ones(size(df))-ones(size(W))*df)),1)) - max_log;

% Determine the support of the posterior by searching to the left and right for when the log unnormalized posterior has decayed sufficiently.
log_ratio = 10; % ten log units should be enough
Fmin = fminbnd(@(x) (log_unnormalized_posterior(x) + log_ratio).^2, df_ML-DW, df_ML); % search from left bound to find where log_unnormalized_posterior = log_ratio
Fmax = fminbnd(@(x) (log_unnormalized_posterior(x) + log_ratio).^2, df_ML, df_ML+DW); % search from right bound to find where log_unnormalized_posterior = log_ratio
%disp(sprintf('F bounds: [%8.3f (%8.3f) %8.3f]', Fmin, df_ML, Fmax));

% Compute maximum of log unnormalized posterior.

% Compute the normalization constant over this range.
Z = quad(@(x) exp(log_unnormalized_posterior(x)), Fmin, Fmax);

% Define the normalized posterior PDF.
pdf = @(df) exp(log_unnormalized_posterior(df)) / Z;

% Compute the posterior mean.
f_mean = quad(@(x) x.*pdf(x), Fmin, Fmax);
%disp(sprintf('F mean: %8.3f', f_mean));

% Compute the cumulative distribution function over a set of finite points.
% This limits the resolution of the confidence bounds, but it's the most efficient solution for now.
npoints = 1000; % number of points
Fpoints = linspace(Fmin, Fmax, npoints);
cdf = cumsum(pdf(Fpoints)); % estimate CDF on sampled points
cdf = cdf / cdf(npoints); % normalize

% Estimate the confidence interval(s).
f_lower = zeros(size(CI));
f_upper = zeros(size(CI));
for i = 1:length(CI)
  f_lower(i) = Fpoints(find(1/2 - CI(i)/2 < cdf, 1));
  f_upper(i) = Fpoints(find(1/2 + CI(i)/2 < cdf, 1));
%  disp(sprintf('CI %8.6f : [%8.3f, %8.3f]', CI(i), f_lower(i), f_upper(i)));  
end

return
