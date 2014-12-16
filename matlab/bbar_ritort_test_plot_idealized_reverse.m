% Test the Bayesian and asymptotic error estimates for BAR by testing many bootstrap subsamples of 
% Ritort dataset appearing in
%
% [1] Collin D, Ritort F, Jarzynski C, Smith SB, Tinoco Jr I, and Bustamante C. Verification of the Crooks fluctuation theorem and recovery of RNA folding free energies. Nature 437:231, 2005.
%
% where forward work samples are drawn from a kernel density estimate of the forward work distribution and
% reverse work samples are drawn from the idealized reverse distribution given the forward distribution and
% the Crooks fluctuation theorem.  This is meant to mimic the effect of eliminating instrument noise.
%
% This function produces a 1D "confidence plot" depicting the fraction of time each estimate falls within 
% various different confidence levels.  A correct posterior should produce a diagonal line.
% Because a finite number of replications are conducted, 95% confidence intervals are plotted to show
% whether discrepancies from x = y are significant.

clear;

% PARAMETERS

% CD4 dataset with 20 pN/s pulling rate (about 500 measurements each).
% Note that the data here doesn't seem to exactly the data appearing in Figure 2 of Ref. [1].
forward_work_datafile = '../datasets/ritort/CD4-20pN_per_s-forward.dat'; % location of datafile containing forward work measurements, comma-delimited, in units of kT
reverse_work_datafile = '../datasets/ritort/CD4-20pN_per_s-reverse.dat'; % location of datafile containing reverse work measurements, comma-delimited, in units of kT

NF = 100; % number of forward realizations per experiment for fixed-number experiment
NR = 100; % number of reverse realizations per experiment for fixed-number experiment

nreplicates = 1000; % number of replications of the experiment to perform (the larger, the smaller the error bars in the plot)

cis = linspace(0.01,0.99,40); % confidence intervals at which to evaluate error

% READ DATA

% Read work measurements from comma-delimited files.
WF = dlmread(forward_work_datafile)';
WR = dlmread(reverse_work_datafile)';

% Do shifts we observe in plots and fix sign of WR.
%shift = 24; % this is estimated from looking at the plots!
shift = 0; % this is estimated from looking at the plots!
WF = WF - shift; % shift forward work measurements down to correspond to plots
WR = -(WR - shift); % shift the reverse work measurement and then negate it

% Compute best estimate of true free energy difference using all data.
[true_df, true_ddf] = BAR(WF, WR);
disp(sprintf('Best estimate from all data using BAR is %f +- %.1f kT', true_df, true_ddf));

% GENERATE KERNEL DENSITY ESTIMATES

% Get sampled range of work values.
Wmin = min([WF' -WR'])-10;
Wmax = max([WF' -WR'])+10;
nx = 1000;
xf = linspace(Wmin, Wmax, nx);
dx = ((Wmax-Wmin)/nx);

% Compute kernel density estimate of WF to get log pdf for forward distribution.
[log_PWF, bandwidth] = kde(WF, xf');
log_PWF = log_PWF';
disp(sprintf('kde bandwidth is %f', bandwidth));

% Compute pdf for forward distribution.
PWF = exp(log_PWF - max(log_PWF)); 
PWF = PWF / sum(PWF);

% Compute cdf for forward distribution.
PWF_cdf = cumsum(PWF);

% Compute pdf for reverse distribution using CFT relation.
% PWR(W) = PWF(-W) exp[beta(W + DF)]
xr = -xf;
log_PWR = log_PWF + xr + true_df;
PWR = exp(log_PWR - max(log_PWR));
PWR = PWR / sum(PWR);

% Compute cdf for reverse distribution.
PWR_cdf = cumsum(PWR);

% find true intersection
index = find(PWR(1:(nx-1)) > PWF(1:(nx-1)) & PWR(2:nx) < PWF(2:nx), 1, 'first');
cross_df = (xf(index) + xf(index+1)) / 2;
disp(sprintf('true_df revised from %f to %f', true_df, cross_df));
true_df = cross_df;

% DEBUG
clf;
hold on;
nbins = 10;
[HN,HX] = hist(WF,nbins);
bar(HX,HN,'r');
[HN,HX] = hist(-WR,nbins);
bar(HX,HN,'b');
plot(xf, PWF / dx * length(WF) * nbins, 'r-', 'LineWidth', 2);
plot(-xr, PWR / dx * length(WR) * nbins, 'b-', 'LineWidth', 2);

% convert number of forward and reverse realizations to probability for fixed-probability experiment
N = NF + NR;  % total number of samples/experiment
PF = NF / N;
PR = NR / N;

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
      % Create fixed-number subsample.
      this_NF = NF;
      this_NR = NR;
    else
      % Create fixed-probability subsample.
      this_NF = binornd(N, PF);
      this_NR = N - this_NF;
    end

    % Generate sample from idealized distributions.
    this_WF = zeros([this_NF,1]);
    this_WR = zeros([this_NR,1]);
    for n = 1:this_NF
      index = find(PWF_cdf < rand, 1, 'last');
      this_WF(n) = xf(index);
    end
    for n = 1:this_NR
      index = find(PWR_cdf < rand, 1, 'last');
      this_WR(n) = xr(index);
    end
  
    % Compute BAR estimate and asymptotic covariance estimate.
    if (phase == 1)
      % Analyze with fixed-number.
      [df, ddf] = ABAR(this_WF, this_WR); % M-factor estimated from number of observed forward/reverse work measurements
    else
      % Analyze with fixed-probability correction.
      [df, ddf] = ABAR(this_WF, this_WR, PF); % M-factor computed from given fixed probability of forward switching events
    end

    % Store maximum likelihood estimate.
    dfA(replicate,phase) = df;
    
    % Compute posterior mean and confidence intervals by Bayesian BAR.
    [df_mean, dfB_lower, dfB_upper] = BBAR(this_WF, this_WR, cis, PF); % M-factor determined from fixed, known probability of forward/reverse work measurements
    
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
filename = sprintf('confidence-gaussian-%d-%d.mat', NF, NR);
save(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY STATISTICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('bias for ABAR: fixed-number %f fixed-probability %f', mean(dfA(:,1) - true_df), mean(dfA(:,2) - true_df)));
disp(sprintf('bias for BBAR: fixed-number %f fixed-probability %f', mean(dfB(:,1) - true_df), mean(dfB(:,2) - true_df)));

disp('');

disp(sprintf('mean unsigned error for ABAR: fixed-number %f fixed-probability %f', mean(abs(dfA(:,1) - true_df)), mean(abs(dfA(:,2) - true_df))));
disp(sprintf('mean unsigned error for BBAR: fixed-number %f fixed-probability %f', mean(abs(dfB(:,1) - true_df)), mean(abs(dfB(:,2) - true_df))));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate confidence level verification plots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set global figure properties.
clf;
%set(gcf, 'position', [0 0 3.25 3.25], 'units', 'inches');

%% Plot work distributions.
%subplot('position', [0.1 0.6 0.8 0.4]);
%hold on;
%wmin = (mu - beta*sigma^2) - 3*sigma;
%wmax = mu + 3*sigma;
%w = linspace(wmin, wmax, 100);
%PWF = normpdf(w, mu, sigma);
%PWR = normpdf(-w, -(mu - beta*sigma^2), sigma);
%plot(w, PWF, 'k-');
%plot(w, PWR, 'k-.');
%fill([w fliplr(w)], [min(PWF,PWR) 0*PWF], 0.7*[1 1 1]);
%xlabel('W');
%ylabel('p(W)');
%set(gca,'YTick', []);
%legend('P_F(W)', 'P_R(-W)');

% DEBUG
subplot('position', [0.1 0.6 0.8 0.4]);
hold on;
nbins = 10;
[HN,HX] = hist(this_WF,nbins);
bar(HX,HN,'r');
plot(xf, PWF * max(HN) / max(PWF), 'r-', 'LineWidth', 2);
[HN,HX] = hist(-this_WR,nbins);
bar(HX,HN,'b');
plot(-xr, PWR * max(HN) / max(PWR), 'b-', 'LineWidth', 2);

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

% Plot BBAR as filled circles with error bars.
errorbar(cis, PB(:,2), PB(:,2) - PB_lower(:,2), PB_upper(:,2) - PB(:,2), 'k.', 'MarkerSize', 10);

% Plot BAR-FP as open circles with error bars.
errorbar(cis, PA(:,2), PA(:,2) - PA_lower(:,2), PA_upper(:,2) - PA(:,2), 'k.', 'MarkerSize', 10);
plot(cis, PA(:,2), 'w.', 'MarkerSize', 5);

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

filename = sprintf('../plots/bbar-ritort-idealized-%d-%d.eps', NF, NR);
print('-depsc', filename);
exportfig(gcf, filename, 'width', 3.25, 'FontMode', 'fixed', 'FontSize', 6.5)
unix(sprintf('epstopdf %s', filename));
