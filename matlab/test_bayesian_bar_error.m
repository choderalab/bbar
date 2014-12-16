% Test the Bayesian and asymptotic error estimates for BAR by performing many replications of an experiment and determining fraction of time the true error is within the anticipated error.

clear;

% PARAMETERS

Nmax = 10; % maximum number of samples N_f + N_r
nreplicates = 1000; % number of replications of the experiment to perform for each combination of (N_f,N_r)
cis = [0.68, 0.95]; % confidence intervals at which to evaluate error

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

% Determine number of confidence intervals to evaluate
ncis = length(cis);

% Compute Gaussian widths for harmonic oscillators
sigma_0 = 1 / sqrt(beta * K_0);
sigma_1 = 1 / sqrt(beta * K_1);

% Define instantaneous work functions
WF = @(x) (K_1/2)*(x-x_1).^2 - (K_0/2)*(x-x_0).^2; % work from state 0 -> 1
WR = @(x) (K_0/2)*(x-x_0).^2 - (K_1/2)*(x-x_1).^2; % work from state 1 -> 0

% Compute true free energy difference.
true_df = - log(sigma_1) + log(sigma_0);

% Plot system.
plot_system;
return

% Estimate error distributions for each combination of (N_f,N_r)
Pfrs_bayesian = zeros([Nmax+1, Nmax+1, ncis]); % Pfrs_bayesian(N_f,N_r,c) is the probability that the true error is within Bayesian confidence interval of cis(c) for N_f forward and N_r reverse samples
Pfrs_asymptotic = zeros([Nmax+1, Nmax+1, ncis]);
mask = zeros([Nmax+1, Nmax+1]);
BAR_ML_bias_fr = zeros([Nmax+1, Nmax+1]);
BAR_mean_bias_fr = zeros([Nmax+1, Nmax+1]);
df_ML = zeros([nreplicates,1]);
df_mean = zeros([nreplicates,1]);
for N_f = 0:Nmax
  for N_r = 0:Nmax     
%    if (N_f > 0) && (N_r > 0) && (N_f + N_r <= Nmax) % only evaluate those cases where there are some, but not too many, samples
%    if (N_f > 0) && (N_r > 0) && (N_f + N_r <= Nmax) % only evaluate those cases where there are some, but not too many, samples
    if (N_f + N_r > 1)
      mask(N_f+1, N_r+1) = 1;
      disp(sprintf('N_f = %d, N_r = %d', N_f, N_r));
      % perform a number of replicates
      Nasymptotic = zeros([nreplicates, ncis]);
      Nbayesian = zeros([nreplicates, ncis]);
      for replicate = 1:nreplicates	
	% perform a replicate of the experiment with fixed number of forward and reverse trajectories
		
	% Draw samples from stationary distributions at inverse temperature beta.
	x_f = sigma_0 * randn([N_f, 1]) + x_0; % samples from state 0
	x_r = sigma_1 * randn([N_r, 1]) + x_1; % samples from state 1
	
	% Compute forward and reverse work values.
	w_f = WF(x_f);
	w_r = WR(x_r);
	
	% compute confidence interval for asymptotic covariance estimate
	[df, ddf] = MBAR(w_f, w_r);
	% store statistics
	df_ML(replicate) = df;
	ddf_ML(replicate) = ddf;

	% Compute mean by Bayesian BAR.
	[df_low, df_this_mean, df_high] = BBAR(w_f, w_r, cis(1), df, ddf);
	df_mean(replicate) = df_this_mean;
	
	% compute free energy confidence intervals
	for c = 1:ncis
	  % get confidence interval cutoff
	  ci = cis(c);

	  % compute confidence interval for Bayesian method
	  [df_low, df_this_mean, df_high] = BBAR(w_f, w_r, ci, df, ddf);
	  % record whether this was in the allowed region
	  if (df_low <= true_df) && (true_df <= df_high)
	    Nbayesian(replicate, c) = Nbayesian(replicate, c) + 1;
	  end
	  
	  % compute normal confidence interval
	  df_ci_lower = df - sqrt(2)*ddf*erfinv(ci);
	  df_ci_upper = df + sqrt(2)*ddf*erfinv(ci);
	  % DEBUG
	  %disp(sprintf('MBAR: df = %8.3f +- %8.3f : [ %8.3f (%8.3f) %8.3f ]', df, ddf, df_ci_lower, true_df, df_ci_upper));
	  
	  % record whether this was in the allowed region
	  if (df_ci_lower <= true_df) && (true_df <= df_ci_upper)
	    Nasymptotic(replicate, c) = Nasymptotic(replicate, c) + 1;
	  end	  	  
	end
	
	% Compute fraction of time true free energy was within various confidence intervals.
	for c = 1:ncis
	  P = sum(Nasymptotic(:,c)) / nreplicates;
	  Pfrs_asymptotic(N_f+1, N_r+1, c) = P;
	  dPfrs_asymptotic(N_f+1, N_r+1, c) =  P * (1-P) / nreplicates;

	  P = sum(Nbayesian(:,c)) / nreplicates;
	  Pfrs_bayesian(N_f+1, N_r+1, c) = P;
	  dPfrs_bayesian(N_f+1, N_r+1, c) =  P * (1-P) / nreplicates;
	end

	% Estimate bias.
	BAR_ML_bias_fr(N_f+1, N_r+1) = mean(df_ML) - true_df;		
	BAR_mean_bias_fr(N_f+1, N_r+1) = mean(df_mean) - true_df;
      end      
      
    end
  end
end

% save all data 
save bar-test.mat

% plot data
generate_plot
