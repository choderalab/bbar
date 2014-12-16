function [f_mean, f_lower, f_upper] =  BBAR_mrsfix(WF, WR, CI)
% ---
% Bayesian BAR for fixed-number experiment.
% ---
%
% Calculate the free energy difference between two states using the
% Bayesian extension of the Bennett acceptance ratio, described in
%
% [1] Maragakis P, Ritort F, Bustamante C, Karplus M, and Crooks GE.
% Bayesian estimates of free energies from nonequilibrium work data in
% the presence of instrument noise. JCP 129:024102, 2008.
%
% In this variant, the M factor is treated as a nuisance parameter, where
%
% M = PF / (1 - PF)
%
% PF ~ Beta(NF+alpha, NR+alpha)
%
% where alpha is a prior hyperparameter, and marginalized out.
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
% PF (optional) probability of forward trajectories to be used in computing M
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

% Determine number of forward and backward samples
NF = length(WF);
NR = length(WR);

% We must have at least one forward and one reverse work value in order for the PDF to be normalizable.
if (NF < 1) | (NR < 1)
  error('There must be at least one forward and one reverse value.');
end

% Compute total number of samples.
N = NF+NR;

% Set the hyperparameter alpha, which sets prior pseudocounts to use for PF ~ Beta(NF+alpha, NR+alpha)
% alpha = 1 (Laplace) or alpha = 1/2 (Jeffreys)
alpha = 1;

% Determine the mean log ratio of samples M.
% Note that the "Laplace prior" form is used for M.
M = log( (NF+alpha) / (NR+alpha) );

% Determine minimum and maximum work values, bracketing the free energy difference estimate.
Wmin = min([WF; -WR]);
Wmax = max([WF; -WR]);
DW = Wmax - Wmin; % interval size

% Define the logarithm of the logistic function in such a way that overflow of the exponential is avoided.
log_logistic = @(x) -log(1 + exp(-x));
%log_logistic = inline('-log(1+exp(-max(x,0))) + min(x,0) - log(exp(+min(x,0))+1) + log(2)', 'x'); 

% Define the logarithm of the unnormalized probability density function.
% This takes a row vector of free energy differences as input.
log_unnormalized_posterior = @(df) sum(log_logistic([WF*ones(size(df)) - ones(size(WF))*df + M; WR*ones(size(df)) + ones(size(WR))*df - M]),1);

% Determine the maximum log likelihood estimate starting from the center of this interval.
df_ML = fminsearch(@(x) -log_unnormalized_posterior(x), (Wmin+Wmax)/2);

% Redefine log unnormalized posterior so that maximum value is zero.
max_log = log_unnormalized_posterior(df_ML);
log_unnormalized_posterior = @(df) sum(log_logistic([WF*ones(size(df)) - ones(size(WF))*df + M; WR*ones(size(df)) + ones(size(WR))*df - M]),1) - max_log;

% Determine the support of the posterior by searching to the left and right for when the log unnormalized posterior has decayed sufficiently.
log_ratio = 10; % ten log units should be enough
Fmin = fminbnd(@(x) (log_unnormalized_posterior(x) + log_ratio).^2, df_ML-DW, df_ML); % search from left bound to find where log_unnormalized_posterior = log_ratio
Fmax = fminbnd(@(x) (log_unnormalized_posterior(x) + log_ratio).^2, df_ML, df_ML+DW); % search from right bound to find where log_unnormalized_posterior = log_ratio
%disp(sprintf('F bounds: [%8.3f (%8.3f) %8.3f]', Fmin, df_ML, Fmax));

% Determine coefficients for conversion from fixed-number to fixed-probability.
A = zeros(N+1,N+1);
for i = 0:N
  for j = 0:N
    A(i+1,j+1) = binopdf(j,N,i/N);
  end 
end

% Determine coefficients for conversion from fixed-probability to fixed-number.
Ainv = inv(A);
disp(Ainv(NF+1,:));

% Redefine log unnormalized posterior for fixed-number experiment to plot.
log_unnormalized_posterior_variable_M = @(df,M) sum(log_logistic([WF*ones(size(df)) - ones(size(WF))*df + M; WR*ones(size(df)) + ones(size(WR))*df - M]),1) - max_log;
clf;
hold on;
Fmin = -4;
Fmax = +5;
nbins = 1000;
df = linspace(Fmin, Fmax, nbins);
avgpdf = zeros(size(df));
for k = 0:N
  M = log(k+1) - log(N-k+1); % compute M-factor
  log_pdf = log_unnormalized_posterior_variable_M(df, M); % compute log pdf
  [max_log_pdf, max_index] = max(log_pdf);
  log_pdf = log_pdf - max_log_pdf;
  Z = sum(exp(log_pdf)) * (Fmax - Fmin) / nbins;
  pdf = exp(log_pdf) / Z;
  weight = Ainv(NF+1,k+1); % weight from coefficient
  avgpdf = avgpdf + weight * pdf;
  plot(df, pdf, 'k-');

  % label
  text(df(max_index), 1, sprintf('%d', k));
end
avgpdf = avgpdf / max(avgpdf);
plot(df, avgpdf, 'r-', 'LineWidth', 3);

M = log(NF) - log(NR);
log_pdf = log_unnormalized_posterior_variable_M(df, M);
pdf = exp(log_pdf - max(log_pdf));
plot(df, pdf, 'k-', 'LineWidth', 3);

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
