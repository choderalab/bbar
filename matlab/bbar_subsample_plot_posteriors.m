% Plot the posteriors of subsampled and full BBAR.

% Test the Bayesian and asymptotic error estimates for BAR by performing many replications of an experiment 
% and determining fraction of time the true error is within the compute confidence bounds from the sample.
%
% This function produces a 1D "confidence plot" depicting the fraction of time each estimate falls within 
% various different confidence levels.  A correct posterior should produce a diagonal line.
% Because a finite number of replications are conducted, 95% confidence intervals are plotted to show
% whether discrepancies from x = y are significant.
% 
% This variant uses a fixed PROBABILITY of forward and backward measurements.

clear;

% PARAMETERS

NF = 10; % number of forward realizations per experiment for fixed-number experiment
NR = 10; % number of reverse realizations per experiment for fixed-number experiment

% convert number of forward and reverse realizations to probability for fixed-probability experiment
N = NF + NR;  % total number of samples/experiment
PF = NF / N;
PR = NR / N;

% DEFINE THE EXPERIMENT HERE
%
% Forward and reverse work measurements are given by Gaussian distributions:
%   WF ~ N(mu, sigma^2)
%   WR ~ N(-(mu + beta sigma^2), sigma^2)
%
% True free energy difference given by
%
% dF = mu - sigma^2 / 2

mu = 1.0; % mean of forward work distribution
sigma = 1.0; % std dev of forward work distribution
beta = 1.0; % inverse temperature

% Compute true free energy difference.
true_df = mu - sigma^2 / 2;

% Conduct fixed-number experiment.
% Compute forward and reverse work values.
WF = mu + sigma * randn([NF, 1]);
WR = -(mu - beta*sigma^2) + sigma * randn([NR, 1]);

% Make figure
clf;
hold on;

% DEBUG
BBAR_beta(WF, WR, 0.95);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute BBAR posterior using all data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the log ratio of samples M.
% Note that the "Laplace prior" form is used for M.
M = log( (NF+1) / (NR+1) );

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
disp(sprintf('F bounds: [%8.3f (%8.3f) %8.3f]', Fmin, df_ML, Fmax));

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

% plot pdf
plot(Fpoints, pdf(Fpoints), 'r-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute BBAR posterior using subsampling without replacement.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subsample without replacement, limiting ourselves to min(NF, NR);
Nmin = min(NF,NR);
this_N_f = binornd(Nmin, PF);
this_N_r = Nmin - this_N_f;
index_f = randperm(NF);
index_r = randperm(NR);
WF = WF(index_f(1:this_N_f));
WR = WR(index_r(1:this_N_r));

% Determine the log ratio of samples M.
% Note that the "Laplace prior" form is used for M.
M = log( (NF+1) / (NR+1) );

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
disp(sprintf('F bounds: [%8.3f (%8.3f) %8.3f]', Fmin, df_ML, Fmax));

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

% plot pdf
plot(Fpoints, pdf(Fpoints), 'g-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% label plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Label plot
legend('all data', 'subsample');
xlabel('\Delta F')
ylabel('p(\Delta F | Data)')

% plot true free energy difference
oldaxis = axis;
plot(true_df * [1 1], oldaxis(3:4), 'k-', 'LineWidth', 2);

return
