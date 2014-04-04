function [m,g,c,d,e,f] = lcmtpl(a,b,acc)

%LCMTPL    Least common multiple.
%
%   LCMTPL(A,B) is the least common multiple of corresponding elements of
%   A and B.  The arrays A and B must be the same size (or either can be scalar).
%
%   LCMTPL(A,B,Accuracy) gives Value of Accuracy to stop loop using GCDIV,
%                         default: EPS
%
%   [M,G,C,D] = LCMTPL(A,B) returns C and D so that G = A.*C + B.*D
%   These are useful for solving Diophantine equations and computing
%   Hermite transformations.
%
%   [M,G,C,D,E,F] = LCMTPL(A,B) returns E and F so that A = E.*G and B = F.*G
%
%   See also GCDIV.
%

%   Copyright (c) 1984-98 by The MathWorks, Inc.  &  cbegler 2001
%   $Revision: 5.6 $  $Date: 1997/11/21 23:45:41 $


if nargin < 3
   acc = eps;
end

[g,c,d,e,f] = gcdiv(a,b,acc);

m = abs(a.*f);
