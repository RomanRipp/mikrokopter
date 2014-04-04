function [g,c,d,e,f] = gcdiv(a,b,acc)

%GCDIV    Greatest common divisor.
%
%   GCDIV(A,B) is the greatest common divisor of corresponding
%   elements of A and B.  The arrays A and B must be the same size 
%    (or either can be scalar).
%
%   GCDIV(0,0) is 0 by convention; all other GCDs are positive integers.
%
%   GCDIV(A,B,Accuracy) gives Value of Accuracy to stop loop, default: EPS
%
%   [G,C,D] = GCDIV(A,B) returns C and D so that G = A.*C + B.*D
%   These are useful for solving Diophantine equations and computing
%   Hermite transformations.
%
%   [G,C,D,E,F] = GCDIV(A,B) returns E and F so that A = E.*G and B = F.*G
%
%   See also LCMTPL.
%

%   Algorithm: See Knuth Volume 2, Section 4.5.2, Algorithm X.
%   Author:    John Gilbert, Xerox PARC
%   Copyright (c) 1984-98 by The MathWorks, Inc.  &  cbegler 2001
%   $Revision: 5.10 $  $Date: 1998/02/17 18:40:23 $

 
% Do scalar expansion if necessary
if length(a) == 1
   a = a(ones(size(b)));
elseif length(b) == 1
   b = b(ones(size(a)));
end

if ~isequal(size(a),size(b))
    error('Inputs must be the same size.')
end;

if nargin < 3
   acc = eps;
end

n = prod(size(a));

e = 0*a;
f = e;
c = e;
d = e;
g = e;
 

for k = 1 : n

   u = [1 0 abs(a(k))];
   v = [0 1 abs(b(k))];

%  while abs(v(3)) > 100*eps

   dg = acc + 1;

   while ( dg > acc )  &  v(3)

       q = round( ( u(3) - mod(u(3),v(3)) ) / v(3) );

       t = u - v*q;

       u = v;
       v = t;
 
       dg = abs( abs( a(k) / ( v(2) + (v(2)==0) ) ) - ...
                 abs( b(k) / ( v(1) + (v(1)==0) ) )       );

   end

   e(k) = abs(v(2)) * sign(a(k));
   f(k) = abs(v(1)) * sign(b(k));

   c(k) = u(1) * sign(a(k));
   d(k) = u(2) * sign(b(k));

%  g(k) = u(3);

   g(k) = ~( v(2) == 0 ) * abs( a(k) / ( v(2) + (v(2)==0) ) ) + ...
          ~( v(1) == 0 ) * abs( b(k) / ( v(1) + (v(1)==0) ) );

   g(k) = g(k) / sum( ~( v([1 2]) == 0 ) );


end

