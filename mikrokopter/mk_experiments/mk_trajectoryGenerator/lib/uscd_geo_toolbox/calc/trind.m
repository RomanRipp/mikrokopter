function p = trind(n)

% TRIND  returns Triangle Indize 
%
% P = TRIND( N ) 
%
% P = [ M by 2 ]  IndexVector,  
%       M = N! / ( K! * (N-K)! ) + N 
%
%------------------------------------------------------
%                  
% see also:  COMBINE, PAIRS, NCHOOSEK, COMBNTNS
%

%***************************************************************
% Check Inputs

Nin = nargin;

if Nin == 0
   error('Not enough InputArguments.');
end

ok = ( isnumeric(n) & ( prod(size(n)) == 1 ) );
if ok
   ok = ( ( n == round(n) ) & ( n >= 0 ) );
end
if ~ok
    error('N must be a positive Integer.');
end

%******************************************************

k = 2;

m = nok(n,k) + n;

p = zeros(m,k);

z = 0;

for ii = 1 : n
 
     m = nok(n-ii,k-1) + 1;
    jj = z + (1:m);

    p(jj,1) = ii;
    p(jj,2) = ( 1 : (n-ii+1) )';

     z = z + m;

end

%******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function n = nok(n,k)

n = prod(n-k+1:n) / prod(1:k);
