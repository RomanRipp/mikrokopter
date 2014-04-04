function [ii,kk] = findextr(x,int)

% FINDEXTR  Find Extrema of Vector
%
%  [ Index , Sign ] = FINDEXTR( X , WindowLength )
%
%  Index  IndexVector of Extrema of X
%  Sign    -1 == Minima
%           1 == Maxima
%
%  The Extrema will searched in a running Window.
%
%  Give a negative WindowLength to search Extrema from End.
%  The Index will returned in descending order.
%
% example:
%
%   x = linspace( 0 , 4*pi , 500 );  % 2 Periods == 2 Minima
%   x = sin(x) - 0.25*x;
%
%   i1 = findextr( x ,  10 );   % Small Window
%   i2 = findextr( x , 200 );   % Large Window
%
%   figure,hold on
%
%   plot(x);
%
%   plot(i1,x(i1),'r*');
%   plot(i2,x(i2),'ko');
%

Nin = nargin;

if Nin < 1
   error('Not enough InputArguments.');
end

if Nin < 2
   int = 1;
end

%--------------------------------------
% Check Inputs

if ~( isnumeric(x)  &  any( prod(size(x)) == size(x) ) )
   error('X must be a numeric Vector');
end

if ~( isnumeric(int)  &  ( prod(size(int)) == 1 ) )
   error('Intervall must be an Integer');
end

if ~isfinite(int)
   error('Intervall must be an Integer');
end

%-------------------------------------
% Analyze Syntax of Intervall
%

int = int + ( int == 0 );

swp = ( int < 0 );  % Swap X

int = floor(abs(int));

int = int + ( int == 0 );

%--------------------------------------

ii = [];
kk = [];

n = prod(size(x));

if n < 3
   return
end

x = x(:);

% Apply  Swap and Operator on X

x = x( (n+1)*swp + (1-2*swp)*(1:n) );

opp = [ -1  1 ];


%-------------------------------------
% Find 1. Extrema

ii = zeros([2 1]);

ii(1) = findmax( opp(1)*x , int );  % Minima
ii(2) = findmax( opp(2)*x , int );  % Maxima

ii = ii + ( n - ii ) .* ( ii == 1 );

[ii,kk] = min(ii);

z = 1;  % Counter

kk = opp(kk);

%-------------------------------------
% Find following Extrema

while ( 1 < ii(z) ) & ( ii(z) < n )
  
      z = z + 1;

      kk(z) = kk(z-1) - 2 * kk(z-1);

      ii(z) = findmax(  kk(z)*x( ii(z-1)+1 : n ) , int ) + ii(z-1);

end

nn = prod(size(ii));

jj = ( 1+(ii(1)==1) : nn-(ii(end)==n) );

ii = ii(jj);
kk = kk(jj);


%---------------------------------------
% Undo Swap in Index

ii = (n+1)*swp + (1-2*swp)*ii;

%---------------------------------------

%************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function i1 = findmax(x,int)

%-------------------------------------
% Initialisation of Window

n = size(x,1);

i1  = 1;
ind = i1 + ( 1 : min(int,n-i1) );

ok  = ( x(ind) >= x(i1) );

%-------------------------------------
% Search for Maxima

while any(ok) & ( i1 < n )

      [mm,jj] = max(x(ind));

      i1 = i1 + jj;

     ind = i1 + ( 1 : min(int,n-i1) );

     ok  = ( x(ind) >= x(i1) );

end


   