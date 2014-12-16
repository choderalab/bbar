function [df, ddf] =  MBAR(WF, WR)
% ---
% MBAR (for two states only)
% ---
%
% Calculates the free energy difference between two states using the multistate Bennett acceptance ratio
% estimator (MBAR) restricted to two states.  This method correctly handles cases where there are no work
% measurements in one direction.
%
% References
%
% [1] Shirts MR and Chodera JD. Statistically optimal analysis of data from multiple
% equilibrium states. Submitted, 2008.
%
% [deltaF,ddeltaF] = MBAR(WF, WR)
%
% Parameters
% WF:       The total work done in forward processes (in dimensionless units)
% WR:       The total work done in reverse processes (in dimensionless units)
%
% Output
% df:       The free energy difference estimate (in dimensionless units)
% ddf:      The uncertainty (estimated standard deviation) of the estimate (in dimensionless units)

% Determine number of forward and backward work measurements.
NF = length(WF);
NR = length(WR);

% Total number of work measurements.
Ntot = NF+NR;

if (NF > 0) & (NR > 0)
  % Initialize with BAR
  [df, ddf] = BAR(WF, WR);
else
  % Set initial estimate of free energy difference to zero.
  df = 0.0;
end

% Define functions.
maxits = 5000;
tol = 1.0e-6;
WFR = [WF; -WR];
for iteration = 1:maxits
  df_old = df;
  f_0 = - log(sum((NF + NR*exp(df-WFR)).^(-1)));
  f_1 = - log(sum((NF*exp(+WFR) + NR*exp(df)).^(-1)));
  df = f_1 - f_0;
  delta = abs(df - df_old);
  if (delta < tol) 
    %disp(sprintf('terminated with estimate df = %f after %d iterations (delta = %e)', df, iteration, delta));
    break
  end
end

%% Estimate the free energies Ft for the time-dependent states
% Determine uncertainties in Ft using MBAR machinery.
W = zeros(Ntot, 2); % weight matrix
W(:,1) = (NF + NR*exp(df-WFR)).^(-1);
W(:,1) = W(:,1) / sum(W(:,1));
W(:,2) = (NF*exp(+WFR) + NR*exp(df)).^(-1);
W(:,2) = W(:,2) / sum(W(:,2));

N = diag([NF NR]);
Theta = W'*pinv(eye(Ntot,Ntot)-W*N*W')*W; % compute covariance matrix by MBAR Eq. 8
variance = Theta(1,1) + Theta(2,2) - 2*Theta(1,2);
ddf = sqrt(variance);

return

[df_bar, ddf_bar] = BAR(WF, WR);
if (variance > 0)
  disp(sprintf('MBAR %12f +- %12f BAR %12f +- %12f   (%8d iterations)', df, ddf, df_bar, ddf_bar, iteration));
else
  disp(sprintf('MBAR %12f +- %12f BAR %12f +- %12f * (%8d iterations)', df, ddf, df_bar, ddf_bar, iteration));
end

return


