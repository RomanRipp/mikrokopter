function [data,nn,dev,c] = meanind2(dat,n,m,p)

% MEANIND2  2-Dimensional Running Mean of Matrice
%
% [ Y , NN , STDEV ] = MEANIND2( X , N , Mode , P );
%
% X       = Matrice to filter 
%
% N       = Size of window, oddinary Integer(s): [ Ny Nx ]
%
% Mode    = Mode of window:  'linear' | 'triangle' | 'cosine' | ...
%                            'gauss'  | 'binomial'
%
% Mode    = positive Numeric for Power of Window in Range  ( -1 .. 1 )
%             Mode == 0    ==>  Impuls
%             Mode == 1    ==>  Triangle
%             Mode == Inf  ==>  Linear (constant)
%
% Mode    = Matriceor Vector  of Weights for window with Size N
%
% P       = Power for pre-normalized window
%
%  defaults:  N    = [ 1  1 ]
%             Mode = 'linear'
%             P    = 1
%
% NaN's in X make no Problems.
%
%------------------------------------------------------
%
% [ Y , NN , STDEV , Coeff ] = MEANIND2( X , N , ... )
%
% returns the Number NN of averaged Elements, the StandartDeviation 
%  per Element and the Coefficients of the selected Window.
%
%------------------------------------------------------
%
% See also: MEANIND1, MFILTER, WINDOW, NOISE
%

Nin  = nargin;
Nout = nargout;

%---------------------------------------------------
% Check Inputs

if Nin < 1
   error('Not enough InputArguments.');
end

si = size(dat);
ps = prod(si);

if ( Nin < 2 ) | ( ps <= 1 )
   data = dat;
   if Nout > 1
    dev =  zeros(si);
     nn = ~isnan(dat);
   end
   return
end

if Nin < 3
   m = 'linear';
end

if Nin < 4
   p = 1;
end

%--------------------------------------
% Check for old version with SET_NAN

if 0 %%% isnumeric(m)
   set_nan = m(1);
   m = 'linear';
else
   set_nan = NaN;
end

%--------------------------------------
% Check N

ok = ( isnumeric(n)  &  ( prod(size(n)) <= 2 ) );
if ok
   ok = all( ( mod(n,2) == 1 ) & ( n > 0 ) );
end
if ~ok
   error('N must be odd Integer(s).');
end

%--------------------------------------

n = n(:)';
n = cat( 2 , n , n );
n = n([1 2]);

sm = size(m); pm = prod(sm); vm = ( max(sm) == pm );

if isnumeric(m) & ( pm > 1 )
   ok = isequal(sm,n);
   if ~ok 
       if vm
          ok = isequal( n , [pm pm] );
          if ~ok
              ii = ( n == pm );
              ok = any(ii);
              if ok
                 ii = find(ii);
                 m = m(:);
                 if ii == 1
                    m = m * ones(1,n(2));
                 else
                    m = ones(n(1),1) * m';
                 end
              end
          else
              m = m(:);
              m = m * m';
          end
          if ~ok 
              error('Size of Vector for Mode correspond with N.')
          end
       else 
          error('Size of Matrice for Mode correspond with N.')
       end
   end
   try
      c = window( n(1) , 'l' , p );
   catch
      error(sprintf('Can''t build Window with Potenz.\n%s',lasterr));
   end
   c = m;
   m = max(m(:));
   c = ( c / (m+(m==0)) ) .^ p;
   m = sum(c(:));
   c = c / (m+(m==0));
else
   try
     c1 = window(n(1),m,p);
     c2 = window(n(2),m,p);
   catch
     error(sprintf('Can''t build Window.\n%s',lasterr));
   end
   c = c1(:) * c2;   
end

if isequal(n,[1 1])
   data = dat;
   if Nout > 1
      dev =  zeros(si);
      nn  = double(~isnan(dat));
   end
   return
end

%---------------------------------------------------
% Expand dat with -n2  ..  +n2 - Elements

n2  = ( n - 1 ) / 2;

s1 = size(dat,1);
s2 = size(dat,2);

i1 = n2(1);
i2 = n2(2);

o1 = ones(1,i1);
o2 = ones(1,i2);

dat  = dat( cat(2,o1,(1:s1),s1*o1) , cat(2,o2,(1:s2),s2*o2) );

dat( cat(2,(1:i1),(1:i1)+(s1+i1)) , : ) = NaN;
dat( : , cat(2,(1:i2),(1:i2)+(s2+i2)) ) = NaN;

%---------------------------------------------------
% Prepare Variables

quo = ~isnan(dat);

dat( find(~quo) ) = 0;

data = zeros(s1,s2);
quot = zeros(s1,s2);

if Nout >= 2
   nn = zeros(s1,s2);
end

i1 = ( 1 : s1 );
i2 = ( 1 : s2 );

for ii1 = 1 : n(1)
    for ii2 = 1 : n(2)
        data = data + c(ii1,ii2)*dat(i1+(ii1-1),i2+(ii2-1));
        quot = quot + c(ii1,ii2)*quo(i1+(ii1-1),i2+(ii2-1));
        if Nout > 1
           nn = nn + quo(i1+(ii1-1),i2+(ii2-1));
        end
    end
end

quo = find(~quo(i1+n2(1),i2+n2(2)));

data = data ./ ( quot + ( quot == 0 ) );

%------------------------------------------
if Nout > 2

   dev = zeros(s1,s2);

   for ii1 = 1 : n(1)
       for ii2 = 1 : n(2)
           dev = dev + (dat(i1+(ii1-1),i2+(ii2-1))-data).^2;
       end
   end

   qadd  = 2*(nn==0)+(nn==1);

   dev = sqrt( dev ./ (nn+qadd-1) );

   dev = dev .* (~(qadd==1));

   dev(quo) = NaN;

end
%------------------------------------------

data(quo) = set_nan;

if Nout > 1
   nn(quo) = 0;
end  

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
% Mode    = positive Numeric for Power of Window: ( -1 .. 1 )
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
  % Triangle | Power

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

