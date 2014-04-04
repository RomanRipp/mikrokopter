function c = bincoeff(n,k)

% BINCOEFF  returns binominal series or single coefficients
%
%--------------------------------------------------------------
%
% C = BINCOEFF( N )  % N'th Order Series, N = [ 0 1 2 3 4 ... ]
%
% C = [ 1 by N+1 ]   % Coefficients
%
% sum(C) = 2^N
%
%--------------------------------------------------------------
%
% C = BINCOEFF( N , K )  (K+1). Coefficient of N'th Order Series
%                     
% K = [ 0 .. N ]
%
% C = N! / ( K! * (N-K)! ) 
%
%--------------------------------------------------------------
% Algorithm
%
% (n,k+1) = (n-k)/(k+1) * (n,k) 
%
% Start with (n,0) = 1, use CUMPROD
%

Nin = nargin;

ok = ( isnumeric(n)  &  ( prod(size(n)) == 1 ) );
if ok
   ok = ( ( mod(n,1) == 0 ) &  ( n >= 0 ) );
end

if ~ok
   error('N must be an non-negative Integer.');
end


if Nin == 1

  k = n;

else

 ok = ( isnumeric(k)  &  ( prod(size(k)) == 1 ) );
 if ok
   ok = ( ( mod(k,1) == 0 ) &  ( 0 <= k ) & ( k <= n ) );
 end

 if ~ok
    error(sprintf('K must be an non-negative Integer between 0 and %.0f.',n));
 end

 % (n,k) == (n,n-k)

 k = k + ( n - 2*k ) * ( k > n/2 );

end


c = ones(2,k+1);

c(1,2:k+1) = n - ( 0 : k-1 );
c(2,2:k+1) = ( 1 : k );

c = round( cumprod( c(1,:) ./ c(2,:) ) );

if Nin == 2
   c = c(k+1);
end
