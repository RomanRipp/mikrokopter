function [p,k] = pairs(n)

% PAIRS  Returns all Combinations of Pairs
%
% P = PAIRS(N)  returns all pairs on N choices
%
% P = [ M by 2 ]  IndexVector,  M = (N-1) * N / 2
%
%------------------------------------------------------
%                  
% PAIRS works about 40 times faster then 
%
%       NCHOOSEK( 1:N , 2 ) (SpecializedFunctions)
%   or  COMBNTNS( 1:N , 2 ) (MappingToolbox)
%
%------------------------------------------------------
%                  
% [P,K] = PAIRS(N)  returns the StartIndex K
%                    for Groups in P(:,1)
%
% K = [ N-1 by 1 ]    
%
% K = cumsum([ 1 ( N-1 : -1 : 2 ) ])
%
% K(N-1) = M
%
% P(K,1) = ( 1 : N-1 )
% P(K,2) = P(K,1) + 1
% P(L,2) = N             L = K( 2 : N-1 ) - 1
%
%------------------------------------------------------
%                  
% see also:  COMBINE, NCHOOSEK, COMBNTNS
%

if nargin < 1
   n = [];
end

if isempty(n)
   n = 0;
end

ok = ( isnumeric(n) & ( prod(size(n)) == 1 ) );
if ok
   ok = ( ( n == round(n) ) & ( n >= 0 ) );
end

if ~ok
    error('N must be a positive Integer.');
end

if n == 0
   p = zeros(0,2);
   k = zeros(0,1);
   return
end

m = n - 1;

p = zeros(m*n/2,2);

  k    = cat( 2 , 1 , ( m : -1 : 2 ) );
  k    = cumsum(k,2);

p(k,1) = 1;
p(:,1) = cumsum(p(:,1),1);

p(:,2) = 1;

p(1,2) = 2;

  l    = k(2:m);

p(l,2) = p(l,1) + 1  - n;

p(:,2) = cumsum(p(:,2),1);
