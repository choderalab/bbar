% Test the Bayesian and asymptotic error estimates for BAR by performing many replications of an experiment 
% and determining fraction of time the true error is within the compute confidence bounds from the sample.
%
% This function produces a 1D "confidence plot" depicting the fraction of time each estimate falls within 
% various different confidence levels.  A correct posterior should produce a diagonal line.
% Because a finite number of replications are conducted, 95% confidence intervals are plotted to show
% whether discrepancies from x = y are significant.


% PARAMETERS

clear;
N_f = 4; % number of forward realizations per experiment
N_r = 4; % number of reverse realizations per experiment

% Generate list of sign permutations.
N = N_f + N_r;
nperms = nchoosek(N, N_f);
disp(sprintf('Generating %d sign permutations...', nperms));
truesigns = [ones(1,N_f) -ones(1,N_r)];
tic; permsigns = enhancedperms([ones(1,N_f) ones(1,N_r)], 'signs'); toc % generate all 2^N sign permutations

nreplicates = 1000; % number of replications of the experiment to perform (the larger, the smaller the error bars in the plot)

cis = linspace(0.01,0.99,99); % confidence intervals at which to evaluate error

% DEFINE THE EXPERIMENT HERE
%
% Harmonic oscillator is given by potential
%   U(x) = (K / 2) (x - x_0)^2
% probability density function
%   p(x) = (1/Z) exp[-\beta U(x)]
%
% Z = \int dx \exp[-\beta U(x)]
%   = \int dx \exp[-(x - x_0)^2 / (2 sigma^2)]
%   = sqrt(2 pi) sigma
% F = - ln Z = - (1/2) ln (2 pi) - ln sigma
% sigma^2 = (\beta K)^{-1}

x_0 = 0.0; % equilibrium spring position for initial state
K_0 = 1.0; % spring constant of initial state
x_1 = 1.5; % equilibrium spring position for final state
K_1 = 4.0; % spring constant of final state
beta = 1.0; % inverse temperature

% Compute Gaussian widths for harmonic oscillators
sigma_0 = 1 / sqrt(beta * K_0);
sigma_1 = 1 / sqrt(beta * K_1);

% Define instantaneous work functions
WF = @(x) (K_1/2)*(x-x_1).^2 - (K_0/2)*(x-x_0).^2; % work from state 0 -> 1
WR = @(x) (K_0/2)*(x-x_0).^2 - (K_1/2)*(x-x_1).^2; % work from state 1 -> 0

% Compute true free energy difference.
true_df = - log(sigma_1) + log(sigma_0);

% Determine number of confidence intervals to evaluate
ncis = length(cis);

% Perform a number of replicates of the experiment, and tally the fraction of time we find the true error is within
% predicted confidence bounds.
NA = zeros([ncis,1]); % NA(c) is the number of realizations for which df_lower <= true_df <= df_upper for the asymptotic BAR (ABAR)
NB = zeros([ncis,1]); % NB(c) is the number of realizations for which df_lower <= true_df <= df_upper for the Bayesian BAR (BBAR)
dfA = zeros([nreplicates,1]); % dfA(r) is the maximum likelihood free energy estimate of replicate r
dfB = zeros([nreplicates,1]); % dfB(r) is the posterior mean free energy estimate of replicate r
for replicate = 1:nreplicates	
  % Perform a replicate of the experiment with fixed number of forward and reverse trajectories

  disp(sprintf('Replicate %d / %d', replicate, nreplicates));
 
  % Draw samples from stationary distributions at inverse temperature beta.
  x_f = sigma_0 * randn([N_f, 1]) + x_0; % samples from state 0
  x_r = sigma_1 * randn([N_r, 1]) + x_1; % samples from state 1
  
  % Compute forward and reverse work values.
  w_f = WF(x_f);
  w_r = WR(x_r);
  
  % Compute BAR estimate and asymptotic covariance estimate.
  [df, ddf] = ABAR(w_f, w_r); % Shirts BAR with asymptotic covariance estimate
  % [df, ddf] = MBAR(w_f, w_r); % two-state version of multistate BAR can treat cases where N_f or N_f == 0.

  % Store maximum likelihood estimate.
  dfA(replicate) = df;
  
  % Compute posterior mean and confidence intervals by Bayesian BAR.
  %[df_mean, dfB_lower, dfB_upper] = BBAR(w_f, w_r, cis);
  [df_mean, dfB_lower, dfB_upper] = BBAR_test(w_f, w_r, cis, permsigns);

  % Store Bayesian estimate.
  dfB(replicate) = df_mean;

  % Determine fraction of time this estimate falls within various confidence intervals.
  for c = 1:ncis
    % Get confidence interval.
    ci = cis(c);
    
    % Determine whether this estimate falls within confidence interval for normal distribution estimate.
    df_ci_lower = df - sqrt(2)*ddf*erfinv(ci);
    df_ci_upper = df + sqrt(2)*ddf*erfinv(ci);
    if (df_ci_lower <= true_df) && (true_df <= df_ci_upper)
      NA(c) = NA(c) + 1;
    end	  	  

    % Determine whether this estimate falls within confidence interval for Bayesian posterior.
    % record whether this was in the allowed region
    if (dfB_lower(c) <= true_df) && (true_df <= dfB_upper(c))
      NB(c) = NB(c) + 1;
    end
    
  end
end
  
% Compute fraction of time true free energy was within various confidence intervals and estimate 95% confidence intervals on this estimate.
PA = NA / nreplicates;
[PA_lower, PA_upper] = beta_confidence_interval(NA, nreplicates, 0.025, 0.975);
PB = NB / nreplicates;
[PB_lower, PB_upper] = beta_confidence_interval(NB, nreplicates, 0.025, 0.975);

% Estimate bias and one-standard-deviation uncertainty.
bias_ABAR = mean(dfA) - true_df;		
dbias_ABAR = std(dfA) / sqrt(nreplicates);
bias_BBAR = mean(dfB) - true_df;
dbias_BBAR = std(dfB) / sqrt(nreplicates);

disp(sprintf('ABAR bias: %f +- %f', bias_ABAR, dbias_ABAR));
disp(sprintf('BBAR bias: %f +- %f', bias_BBAR, dbias_BBAR));

% Save all data for later analysis.
filename = sprintf('confidence-%d-%d.mat', N_f, N_r);
save(filename);

% Generate a plot for visualization.
clf;
plot([0 1], [0 1], 'k-');
hold on;
errorbar(cis, PA, PA - PA_lower, PA_upper - PA, 'r.');
errorbar(cis, PB, PB - PB_lower, PB_upper - PB, 'g.');
xlabel('confidence interval');
ylabel('observed fraction');
title(sprintf('N_f = %d, N_r = %d', N_f, N_r));
legend('x = y', 'ABAR', 'BBAR');
axis square;
axis([0 1 0 1]);

% Save the plot as PDF.
filename = sprintf('confidence-allperms-%d-%d.eps', N_f, N_r);
print('-depsc', filename);
unix(sprintf('epstopdf %s', filename));

