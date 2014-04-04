function [c,x] = window(n,m,p)

% WINDOW  returns an symetric, normalized window
%
% [C,X] = WINDOW( N , Mode , Power )
%
% N       = Length of window, oddinary Integer
% N < 0   = Right part of Window is ZERO
% 
% Mode    = Mode of window: 'linear' | 'triangle' | 'cosine' | ...
%                                      'gauss'    | 'binomial'
%
% Mode    = positive Numeric for Power of Window in Range ( -1 .. 1 )
%             Mode == 0    ==>  Impuls
%             Mode == 1    ==>  Triangle
%             Mode == Inf  ==>  Linear (constant)
% 
% Power   = Power for pre-normalized ( max == 1 ) window
%
%
% C       = normalized window, [ 1 by N ], sum(C) == 1 
%
% X       = X-Coordinates, [ 1 by N ], [ -1 .. 1 ]
%
%
%  defaults:  N = 1
%             Mode = 'linear'
%             Power = 1
%
% For Demonstration see WINDEMO, type: >> windemo 
%

%---------------------------------------------------
% Check Inputs

msg = '';

nl = char(10);

Nin = nargin;

%---------------------------------------------------
% Check Length

if Nin < 1

   n = 1;

else

  ok = ( isnumeric(n)  &  ( prod(size(n)) == 1 ) );
  if ok
     ok = ( mod(n,2) == 1 );
  end

  if ~ok
     msg = 'N must be an odd Integer.';
  end

end

%---------------------------------------------------
% Check Mode

e = 1;

if Nin < 2

  m = 'l';

else

  ok = ( ischar(m) & ( prod(size(m)) == size(m,2) ) & ~isempty(m) );

  if ~ok
      ok = ( isnumeric(m)  &  ( prod(size(m)) == 1 ) );
      if ok
         ok = ( ~isnan(m) & ( m >= 0 ) );
      end
      if ok
         e = m;
         m = 'p';
      end
  end

  if ~ok
      msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                 'Mode must be a String or positive Numeric.' );
  end

end

%---------------------------------------------------
% Check Potenz

if Nin < 3

  p = 1;

else

  ok = ( isnumeric(p)  &  ( prod(size(p)) == 1 ) );
  if ok
     ok = isfinite(p);
  end

  if ~ok
     msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                'Potenz must be a finite Numeric.'    );
  end

end


%---------------------------------------------------

if ~isempty(msg)
   error(msg)
end

%---------------------------------------------------
% Build Window

is2 = ( n < 0 );
n   = abs(n);

if n == 1
   x = 0;
   c = 1;
   return
end

x = ( 1 : n );

switch lower(m(1))

  %-------------------------------------------------
  % Binomial

  case 'b'

    nc = n - 1;

    c = ones(2,nc+1);

    c(1,2:nc+1) = nc - ( 0 : nc-1 );
    c(2,2:nc+1) = ( 1 : nc );

    c = round( cumprod( c(1,:) ./ c(2,:) ) );

  %-------------------------------------------------
  % Cosine

  case 'c'

    c = 1 - cos( 2*pi * x/(n+1) );

  %-------------------------------------------------
  % Gauss

  case 'g'
 
    c = pi;  % exp(1);

    c = c * ( 2 * x / (n+1) - 1 );  % linspace(-c,c,n+2)

    c = exp( (-1) * c.^2 / 2 );

  %-------------------------------------------------
  % Triangle | Potenz

  case { 't'  'p' }

    if isinf(e)

       c = ones(1,n);

    else

       c = 2 * x / (n+1) - 1;

       if e == 0
          c = double( c == 0 );
       elseif e == 1
          c = 1 - abs(c);
       else    
          c = 1 - abs(c) .^ e;
       end

    end

  %-------------------------------------------------
  % linear 
 
  otherwise


    c = ones(1,n);


end

x = 2 * (x-1) / (n-1) - 1;  % linspace(-1,1,n)

%---------------------------------------------------
% Half Window

if is2
   c( (n+1)/2+1 : n ) = 0;
end

%---------------------------------------------------
% Normalize

if ( p == 1 )
  c = c / sum(c);
  return
end

%---------------------------------------------------
% Potenz and Normalize

c = ( c / max(c) ) .^ p;

c = c / sum(c);


