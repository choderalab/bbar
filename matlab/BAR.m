function [deltaF,ddeltaF] =  BAR(WF, WR, type)
% ---
% BAR
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
% [deltaF,ddeltaF] = BAR(WF,WR)
%
% Parameters
% WF:       The total work done in forward processes (in dimensionless units)
% WR:       The total work done in reverse processes (in dimensionless units)
% type: (optional) either 'fixed-number' or 'fixed-probability'
%
% Output
% deltaF:   The free energy difference estimate (in dimensionless units)
% ddeltaF:  The uncertainty (estimated standard deviation) of the estimate (in dimensionless units)

% SET DEFAULTS
if (nargin < 3) | isempty(type) 
  % use correction for fixed number
  type = 'fixed-number';
end

NF = length(WF);
NR = length(WR);
N = NF+NR;
M = log(NF/NR);

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

% Compute estimate of free energy difference by iterating Eq. 1 from MRS BAR paper.
%maxits = 5000;
%tol = 1.0e-7;
%for iteration = 1:maxits
%  numerator = mean(fermi(M + WF - deltaF));
%  denominator = mean(exp(-WR) .* fermi(M - WR - deltaF));
%  deltaF_old = deltaF;
%  deltaF = - log(numerator) + log(denominator);
%  delta = abs(deltaF - deltaF_old) / abs(deltaF);
%  if (delta < tol) 
%    %disp(sprintf('terminated after %d iterations', iteration));
%    break
%  end
%end

% Compute one-standard-deviation uncertainty
%ddeltaF = sqrt((mean((1./(2+2*cosh((M + [WF; -WR] - deltaF)))))^-1 - (N/NF + N/NR))/N);

% WARNING: This could be running into numerical issues.
expect = (1/N) * ( sum((2. + 2.*cosh(M + WF - deltaF)).^(-1)) + sum((2. + 2*cosh(M - WR - deltaF)).^(-1)) );
variance = (1/N) * (1/expect);
if (type == 'fixed-number')
  % use correction for fixed number
  variance = variance - (1/NF + 1/NR);
end

if (variance < 0)
  disp(sprintf('variance = %f = (1/N) * (%f - %f)', variance, 1/expect, (N/NF + N/NR)));
  ddeltaF = 0.0;
else
  ddeltaF = sqrt(variance);
end

%disp(sprintf('deltaF = %f +- %f', deltaF, ddeltaF));

return


