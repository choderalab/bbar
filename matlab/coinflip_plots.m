% Demonstrate Bayesian posterior estimate for inference of biased coin probability from a number of coin flips.
%
% This function produces a 1D "confidence plot" depicting the fraction of time each estimate falls within 
% various different confidence levels.  A correct posterior should produce a diagonal line.
% Because a finite number of replications are conducted, 95% confidence intervals are plotted to show
% whether discrepancies from x = y are significant.

%clear;

% PARAMETERS

nreplicates = 1000; % number of replications of the experiment to perform (the larger, the smaller the error bars in the plot)

cis = linspace(0.01,0.99,40); % confidence intervals at which to evaluate error

% DEFINE THE EXPERIMENT HERE

N = 1000; % total number of coin flips
PH_true = 0.3; % true probability of heads

alpha = 1/2; % hyperparameters for beta prior -- 1/2 corresponds to Jeffreys prior
beta = alpha;

fontsize = 8;

clf;

set(gcf, 'Units', 'inches', 'position', [0 0 3.25 3.25]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate plot of example posteriors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot('position', [0.1 0.6 0.8 0.4]);

% Run a coin flip experiment.
NH = binornd(N, PH_true); NT = N - NH; % coin flip experiment

% Set extents of pdf plot to encompass 99.9% confidence interval.
ci = 0.999;
x_low = betainv(0.5 - ci/2, NH+alpha, NT+beta);
x_high = betainv(0.5 + ci/2, NH+alpha, NT+beta);
x = linspace(x_low, x_high, 100);

hold on;

y = betapdf(x, NH+alpha, NT+beta);

% Shade 95% confidence interval.
ci = 0.95; % confidence interval
gray = 0.7; % gray level
PH_low = betainv(0.5 - ci/2, NH+alpha, NT+beta);
PH_high = betainv(0.5 + ci/2, NH+alpha, NT+beta);
indices = find((x >= PH_low) & (x <= PH_high));
fill([x(indices) x(fliplr(indices))], [y(indices) 0*y(indices)], gray*[1 1 1], 'LineStyle', 'none');

% Plot posterior pdf.
plot(x, y, 'k-', 'LineWidth', 2);

% Plot posteriors for several of additional realizations.
nadditional = 5;
for index = 1:nadditional 
  NH = binornd(N, PH_true); NT = N - NH;
  y = betapdf(x, NH+alpha, NT+beta);
  plot(x, y, 'k:', 'LineWidth', 1);
end

% Adjust and label plot.
%axis square;
oldaxis = axis;
axis([x_low x_high 0 oldaxis(4)]);
xlabel('\theta');
ylabel('p(\theta | D)');
set(gca, 'YTick', [0]);

% Show true PH.
plot(PH_true * [1 1], [0 oldaxis(4)], 'k-', 'LineWidth', 2);

% Turn on enclosing box.
set(gca, 'box', 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute observed confidence levels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine number of confidence intervals to evaluate
ncis = length(cis);

% Perform a number of replicates of the experiment, and tally the fraction of time we find the true error is within
% predicted confidence bounds.
NB = zeros([ncis,1]); % NB(c) is the number of realizations for which PH_low <= PH <= PH_high for the Bayesian posterior
dfB = zeros([nreplicates,1]); % dfB(r) is the posterior mean estimate of PH_true
for replicate = 1:nreplicates	
  % Perform a replicate of the experiment.
  if (mod(replicate, 10) == 0)
    disp(sprintf('Replicate %d / %d', replicate, nreplicates));
  end
 
  % Perform coin flips.
  NH = binornd(N, PH_true); % number of heads
  NT = N - NH; % number of tails

  % Compute maximum-likelihood estimate of PH.
  PH_mle = NH / N;
  
  % Compute Bayesian posterior mean.
  PH_mean = (NH + alpha) / ((NH + alpha) + (NT + beta));

  % Compute Bayesian posterior confidence intervals.
  PH_low = betainv(0.5 - cis/2, NH+alpha, NT+beta);
  PH_high = betainv(0.5 + cis/2, NH+alpha, NT+beta);

  % Accumulate counts whether true probability of heads PH lies in confidence intervals.
  indices = find((PH_low <= PH_true) & (PH_true <= PH_high));
  NB(indices) = NB(indices) + 1;
end
  
% Estimate fraction of replicates where true free energy was within various confidence intervals
% and estimate 95% confidence intervals for this estimator.
PB = NB / nreplicates;
[PB_lower, PB_upper] = beta_confidence_interval(NB, nreplicates, 0.025, 0.975);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate confidence level verification plot for coinflip experiment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot('position', [0.5 0.1 0.4 0.4]);

hold on;

% Plot X = Y line.
plot([0 1], [0 1], 'k-', 'LineWidth', 2); % x = y line

% Plot dots with error bars.
errorbar(cis, PB, PB - PB_lower, PB_upper - PB, 'k.', 'MarkerSize', 10);
plot(cis, PB, 'w.', 'MarkerSize', 5);

% Plot filled 95% confidence region -- looks poor here.
%plot(cis, PB, 'k:', 'LineWidth', 2);
%fill([cis fliplr(cis)], [PB_upper' fliplr(PB_lower')], gray * [1 1 1]);

% Label plot
%xlabel('desired confidence level');
%ylabel('actual confidence level');
%title(sprintf('coinflip (N = %d)', N));
%axis square;
axis([0 1 0 1]);
set(gca, 'box', 'on');
set(gca, 'XTick', [1]);
set(gca, 'YTick', []);
%set(gca, 'FontSize', fontsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate reference confidence level verification plot for normal distribution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot('position', [0.1 0.1 0.4 0.4]);

% Generate plot box.
axis([0 1 0 1]);
%axis square;
hold on;

% Generate normal deviations.
%sigmas = [1, 2, 10, 100]; % standard deviations to plot
sigmas = [100, 10, 2, 1, 1/2, 1/10, 1/100]; % standard deviations to plot
npoints = 1000; % number of points in each line
desired_cis = linspace(0, 1, npoints); % desired confidence intervals
handles = [];
for sigma = sigmas
  % Generate plot for under- and over-estimates of normal posterior by 'sigma' standard deviations.

  thetas_lower = norminv(0.5 - desired_cis/2);
  thetas_upper = norminv(0.5 + desired_cis/2);
  actual_cis = normcdf(thetas_upper, 0, sigma) - normcdf(thetas_lower, 0, sigma);

  handle = plot(desired_cis, actual_cis, 'k-');
  handles = [handles handle];

  if sigma == 1
    set(handle, 'LineWidth', 2);
  end
end

%legend('100 \sigma', '10 \sigma', '2 \sigma', '\sigma', '1/2 \sigma', '1/10 \sigma', '1/100 \sigma');

% Label curves.
fontsize = 6.5;
text(0.0250, 0.9701, '100 \sigma', 'FontSize', fontsize);
text(0.1506, 0.8562, '10 \sigma', 'FontSize', fontsize);
text(0.3520, 0.6152, '2 \sigma', 'FontSize', fontsize);
text(0.4370, 0.5044, '\sigma', 'FontSize', fontsize);
text(0.5096, 0.3974, '1/2 \sigma', 'FontSize', fontsize);
text(0.6646, 0.1418, '1/10 \sigma', 'FontSize', fontsize);
text(0.7319, 0.0395, '1/100 \sigma', 'FontSize', fontsize);

% Label plot
%xlabel('desired confidence level');
text(0.65, -0.15, 'desired confidence level');
ylabel('actual confidence level');
set(gca, 'XTick', [0 1]);
set(gca, 'YTick', [0 1]);
%axis square;
axis([0 1 0 1]);
set(gca, 'box', 'on');
%set(gca, 'FontSize', fontsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE THE PLOT AS PDF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = sprintf('../plots/coinflip.eps');
print('-depsc', filename);
exportfig(gcf, filename, 'width', 3.25, 'FontMode', 'fixed', 'FontSize', 6.5)
unix(sprintf('epstopdf %s', filename));
