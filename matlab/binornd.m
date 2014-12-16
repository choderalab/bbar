function r=binornd(n,p,mm,nn)
% BINORND Random matrices from a binomial distribution.
%   R = BINORND(N,P,MM,NN)  is an MM-by-NN matrix of random
%   numbers chosen from a binomial distribution with parameters N and P.
%
%   The size of R is the common size of N and P if both are matrices.
%   If either parameter is a scalar, the size of R is the size of the other
%   parameter. 
%   Alternatively, R = BINORND(N,P,MM,NN) returns an MM by NN matrix. 
%
%   The method is direct using the definition of the binomial
%   distribution as a sum of Bernoulli random variables.

%   References:
%      [1]  L. Devroye, "Non-Uniform Random Variate Generation", 
%      Springer-Verlag, 1986
%   See Lemma 4.1 on page 428, method on page 524.

%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 2.11 $  $Date: 2002/03/31 22:26:57 $

if nargin == 2
    [errorcode rows columns] = rndcheck(2,2,n,p);
elseif nargin == 3
    [errorcode rows columns] = rndcheck(3,2,n,p,mm);
elseif nargin == 4
    [errorcode rows columns] = rndcheck(4,2,n,p,mm,nn);
else
    error('Requires at least two input arguments.'); 
end

if errorcode > 0
    error('Size information is inconsistent.');
end

% Handle the scalar params case efficiently
if prod([size(p),size(n)]) == 1 % scalar params
    if (0 <= p & p <= 1) & (0 <= n & round(n) == n)
        r = sum(rand(rows,columns,n) < p, 3);
    else
        r = repmat(NaN, rows, columns);
    end

% Handle the scalar n case efficiently
elseif prod(size(n)) == 1
    if 0 <= n & round(n) == n
        r = sum(rand(rows,columns,n) < repmat(p, [1,1,n]), 3);
        r(p < 0 | 1 < p) = NaN;
    else
        r = repmat(NaN, rows, columns);
    end

% Handle the non-scalar params case
else 
    if prod(size(p)) == 1, p = repmat(p, rows, columns); end
    r = zeros(rows,columns);
    for i = 1:max(max(n))
        k = find(n >= i);
        r(k) = r(k) + (rand(size(k)) < p(k));
    end
    r(p < 0 | 1 < p | n < 0 | round(n) ~= n) = NaN;
end
