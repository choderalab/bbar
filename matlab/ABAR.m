function [deltaF,ddeltaF] =  ABAR(WF, WR, PF)
% ---
% ABAR
% ---
%
% Calculates the free energy difference between two states using the
% classic Bennett Acceptance Ratio method, as extended by Crooks.  The
% equation is solved by finding the zero of the implicit function, as
% decribed by Shirts et al.
%
% References
%
% [1] C. Bennett. Efficient Estimation of Free Energy Differences from Monte
% Carlo Data. Journal of Computational Physics 22, 245-268 (1976).
% [2] G. Crooks. Path-ensemble averages in systems driven far from equilibrium.
% Physical Review E 61, 2361-2366 (2000).
% [3] M. Shirts, E. Bair, G. Hooker, and V. Pande.  Equilibrium Free Energies
% from Nonequilibrium Measurements Using Maximum-Likelihood Methods.
% Physical Review Letters 91, 140601 (2003).
%
% USAGE
% 
% [deltaF,ddeltaF] = ABAR(WF,WR)
%
% PARAMETERS
%
% WF:       The total work done in forward processes (in dimensionless units)
% WR:       The total work done in reverse processes (in dimensionless units)
%
% OPTIONAL PARAMETERS
%
% PF        If specified, the probability of a forward measurement to be used in computing M. Disables correction for fixed number.
% 
% RETURNS
%
% deltaF:   The free energy difference estimate (in dimensionless units)
% ddeltaF:  The uncertainty (estimated standard deviation) of the estimate (in dimensionless units)
%
% This function was written by David D. L. Minh, NIH
% and was downloaded from here: http://mccammon.ucsd.edu/~dminh/software/
% and modified by John D. Chodera.

NF = length(WF);
NR = length(WR);
N = NF+NR;
M = log(NF/NR);

% Override M if probability of forward trajectories is specified
if (nargin > 2)
  M = log ( PF / (1-PF) );
end  

%%% Use Jarzynski's equality in the forward direction as
%%% an initial estimate of deltaF
deltaF = -log(mean(exp(-WF)));

% use zero as initial estimate of deltaF.
deltaF = 0.0;

fermi = @(x) 1./(1+exp(x));
dlogL = @(dF) sum(fermi(+(M + WF*ones(size(dF)) - ones(size(WF))*dF))) - ...
              sum(fermi(-(M - WR*ones(size(dF)) - ones(size(WR))*dF)));

% Solve using implicit equation from Eq. 8 of Ref. [1].
deltaF = fzero(dlogL,deltaF);

% Compute variance.
% Warning: this code may be vulnerable to numerical issues.
expect = (1/N) * ( sum((2. + 2.*cosh(M + WF - deltaF)).^(-1)) + sum((2. + 2*cosh(M - WR - deltaF)).^(-1)) );
variance = (1/N) * (1/expect);
% If a fixed probability PF of forward/reverse measurements is not specify, add correction for fixed number of forward/reverse measurements.
if (nargin == 2)
  variance = variance - (1/NF + 1/NR);
end

if (variance < 0)
  disp(sprintf('variance = %f = (1/N) * (%f - %f)', variance, 1/expect, (N/NF + N/NR)));
  ddeltaF = 0.0;
else
  ddeltaF = sqrt(variance);
end

return


