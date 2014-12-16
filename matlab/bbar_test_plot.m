% Test the Bayesian and asymptotic error estimates for BAR by performing many replications of an experiment 
% and determining fraction of time the true error is within the compute confidence bounds from the sample.
%
% This function produces a 1D "confidence plot" depicting the fraction of time each estimate falls within 
% various different confidence levels.  A correct posterior should produce a diagonal line.
% Because a finite number of replications are conducted, 95% confidence intervals are plotted to show
% whether discrepancies from x = y are significant.
% 
% This variant uses a fixed PROBABILITY of forward and backward measurements.

%clear;

% PARAMETERS

N_f = 10; % number of forward realizations per experiment for fixed-number experiment
N_r = 10; % number of reverse realizations per experiment for fixed-number experiment

nreplicates = 1000; % number of replications of the experiment to perform (the larger, the smaller the error bars in the plot)

cis = linspace(0.01,0.99,40); % confidence intervals at which to evaluate error

% convert number of forward and reverse realizations to probability for fixed-probability experiment
N_tot = N_f + N_r;  % total number of samples/experiment
P_f = N_f / N_tot;
P_r = N_r / N_tot;

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
NA = zeros([ncis,2]); % NA(c) is the number of realizations for which df_lower <= true_df <= df_upper for the asymptotic BAR (ABAR)
NB = zeros([ncis,2]); % NB(c) is the number of realizations for which df_lower <= true_df <= df_upper for the Bayesian BAR (BBAR)
dfA = zeros([nreplicates,2]); % dfA(r) is the maximum likelihood free energy estimate of replicate r
dfB = zeros([nreplicates,2]); % dfB(r) is the posterior mean free energy estimate of replicate r

for phase = 1:2 % FN, FP
  for replicate = 1:nreplicates	
    % Perform a replicate of the experiment with fixed number of forward and reverse trajectories
    
    if (mod(replicate,10) == 0)
      disp(sprintf('Phase %d/2 replicate %d / %d', phase, replicate, nreplicates));
    end
    
    if (phase == 1)
      % Conduct fixed-number experiment.
      this_N_f = N_f;
      this_N_r = N_r;
    else
      % Conduct fixed-probability experiment.
      this_N_f = binornd(N_tot, P_f);
      this_N_r = N_tot - this_N_f; 
    end
 
    % Draw samples from stationary distributions at inverse temperature beta.
    x_f = sigma_0 * randn([this_N_f, 1]) + x_0; % samples from state 0
    x_r = sigma_1 * randn([this_N_r, 1]) + x_1; % samples from state 1
  
    % Compute forward and reverse work values.
    w_f = WF(x_f);
    w_r = WR(x_r);
  
    % Compute BAR estimate and asymptotic covariance estimate.
    if (phase == 1)
      % Analyze with fixed-number.
      [df, ddf] = ABAR(w_f, w_r); % M-factor estimated from number of observed forward/reverse work measurements
    else
      % Analyze with fixed-probability correction.
      [df, ddf] = ABAR(w_f, w_r, P_f); % M-factor computed from given fixed probability of forward switching events
    end

    % Store maximum likelihood estimate.
    dfA(replicate,phase) = df;
    
    % Compute posterior mean and confidence intervals by Bayesian BAR.
    [df_mean, dfB_lower, dfB_upper] = BBAR(w_f, w_r, cis, P_f); % M-factor determined from fixed, known probability of forward/reverse work measurements
    
    % Store Bayesian estimate.
    dfB(replicate,phase) = df_mean;
    
    % Determine fraction of time this estimate falls within various confidence intervals.
    for c = 1:ncis
      % Get confidence interval.
      ci = cis(c);
      
      % Determine whether this estimate falls within confidence interval for normal distribution estimate.
      df_ci_lower = df - sqrt(2)*ddf*erfinv(ci);
      df_ci_upper = df + sqrt(2)*ddf*erfinv(ci);
      if (df_ci_lower <= true_df) && (true_df <= df_ci_upper)
	NA(c,phase) = NA(c,phase) + 1;
      end	  	  

      % Determine whether this estimate falls within confidence interval for Bayesian posterior.
      % record whether this was in the allowed region
      if (dfB_lower(c) <= true_df) && (true_df <= dfB_upper(c))
	NB(c,phase) = NB(c,phase) + 1;
      end
      
    end
  end
end
  
% Compute fraction of time true free energy was within various confidence intervals and estimate 95% confidence intervals on this estimate.
PA = NA;
PB = NB;
PA_lower = NA;
PA_upper = NA;
PB_lower = NB;
PB_upper = NB;
for phase = 1:2
  PA(:,phase) = NA(:,phase) / nreplicates;
  [PA_lower(:,phase), PA_upper(:,phase)] = beta_confidence_interval(NA(:,phase), nreplicates, 0.025, 0.975);
  PB(:,phase) = NB(:,phase) / nreplicates;
  [PB_lower(:,phase), PB_upper(:,phase)] = beta_confidence_interval(NB(:,phase), nreplicates, 0.025, 0.975);
end

% Save all data for later analysis.
filename = sprintf('confidence-%d-%d.mat', N_f, N_r);
save(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate confidence level verification plots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set global figure properties.
clf;
set(gcf, 'position', [0 0 3.25 3.25], 'units', 'inches');

% Choose fixed-number plot.
subplot('position', [0.1 0.1 0.4 0.4]);

hold on;

% Plot X = Y line.
plot([0 1], [0 1], 'k-', 'LineWidth', 2); % x = y line

% Plot BAR-FN as open circles with error bars.
errorbar(cis, PA(:,1), PA(:,1) - PA_lower(:,1), PA_upper(:,1) - PA(:,1), 'k.', 'MarkerSize', 10);
plot(cis, PA(:,1), 'w.', 'MarkerSize', 5);

% Plot BBAR as filled circles with error bars.
errorbar(cis, PB(:,1), PB(:,1) - PB_lower(:,1), PB_upper(:,1) - PB(:,1), 'k.', 'MarkerSize', 10);

% Label plot
ylabel('actual confidence level');
axis([0 1 0 1]);
set(gca, 'box', 'on');
set(gca, 'XTick', [0 1]);
set(gca, 'YTick', [0 1]);

title('fixed-number');

% Choose fixed-probability plot.
subplot('position', [0.5 0.1 0.4 0.4]);

hold on;

% Plot X = Y line.
plot([0 1], [0 1], 'k-', 'LineWidth', 2); % x = y line

% Plot BAR-FP as open circles with error bars.
errorbar(cis, PA(:,2), PA(:,2) - PA_lower(:,2), PA_upper(:,2) - PA(:,2), 'k.', 'MarkerSize', 10);
plot(cis, PA(:,2), 'w.', 'MarkerSize', 5);

% Plot BBAR as filled circles with error bars.
errorbar(cis, PB(:,2), PB(:,2) - PB_lower(:,2), PB_upper(:,2) - PB(:,2), 'k.', 'MarkerSize', 10);

% Label plot
axis([0 1 0 1]);
set(gca, 'box', 'on');
set(gca, 'XTick', [1]);
set(gca, 'YTick', [0 1]);

title('fixed-probability');

text(0.65-1, -0.15, 'desired confidence level');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE THE PLOT AS PDF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = sprintf('../plots/bbar-%d-%d.eps', N_f, N_r);
print('-depsc', filename);
exportfig(gcf, filename, 'width', 3.25, 'FontMode', 'fixed', 'FontSize', 6.5)
unix(sprintf('epstopdf %s', filename));
