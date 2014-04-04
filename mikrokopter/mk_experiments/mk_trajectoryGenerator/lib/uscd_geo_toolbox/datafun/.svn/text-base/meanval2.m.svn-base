function data = meanval2(dat,n,w,m,v,p)

% MEANVAL2  Running Value-Mean along 1. and 2. Dimension
%
% C = MEANVAL2( Z , N , W , Mode , V , P );
%
% N = Length of window, oddinary Integer(s): [ Ny Nx ]
%
% W = Inputs for N-based window: { Mode  Potenz }
%
%    Mode = Mode of window: {'linear'} | 'triangle' | 'cosine' | ...
%                            'gauss'  | 'binomial'
%
%    Mode = positive Numeric for Power of Window: ( -1 .. 1 )
%             Mode == 0    ==>  Impuls
%             Mode == 1    ==>  Triangle
%             Mode == Inf  ==>  Linear (constant)
%
%    Potenz = Potenz for pre-normalized window, default: 1
%
%----------------------------------------------------------------------
%
% Mode = Mode of Value-based window: {'triangle'} | 'cosine' |  'gauss'   
%
% Mode = positive Numeric for Power of Window: ( -1 .. 1 )
%             Mode == 0    ==>  Impuls
%             Mode == 1    ==>  Triangle
%             Mode == Inf  ==>  Linear (constant)
% 
% V = Deviation + Offset*i
%
%     Deviation = shift of Value-based window:
%                  -1 left    ¯¯\_          1 ---> 0  focus on MIN
%         default:  0 center  _/\_   0 ---> 1 ---> 0  focus on MEAN
%                   1 right   _/¯¯   0 ---> 1         focus on MAX
%
%     Offset    = shift of ZERO-Value of window 
%                  0 < Offset <= 1; default: 1/N
%
% P = Potenz of Offset-shifted, Value-based window, default: 1
%
%----------------------------------------------------------------------
%
% NaN's in Z make no Problems.
%
 
Nin  = nargin;
Nout = nargout;

if Nin < 2
   data = dat;
   return
end

%***************************************************************
% Check Inputs

msg = cell(0,1);

%---------------------------------------------------
% Check N

  ok = ( isnumeric(n)  &  ( prod(size(n)) <= 2 ) & ~isempty(n) );
  if ok
      n = abs(real(n(:)'));
     ok = all( mod(n,2) == 1 );
  end
  if ~ok
      msg = cat( 1 , msg , {'N must be an odd Integer with max. 2 Elements.'} );
  end

%---------------------------------------------------
% WindowInput

if Nin < 3
   w = [];
end

%---------------------------------------------------
% Check Mode

e = 1;

if Nin < 4
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
      msg = cat( 1 , msg , {'Mode must be a String or positive Numeric.'} );
  end
end

%---------------------------------------------------
% Check V

if Nin < 5
   v = 0;
else
  ok = ( isnumeric(v)  &  ( prod(size(v)) == 1 ) );
  if ok
     ok = isfinite(v);
  end
  if ~ok
     msg = cat( 1 , msg , {'V must be a single finite Numeric.'} );
  end
end

%---------------------------------------------------
% Check Potenz

if Nin < 6
   p = 1;
else
  ok = ( isnumeric(p)  &  ( prod(size(p)) == 1 ) );
  if ok
     ok = isfinite(p);
  end
  if ~ok
     msg = cat( 1 , msg , {'Potenz must be a single finite Numeric.'} );
  end
end

%---------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Invalid Inputs.\n%s',msg));
end

%*****************************************************************

%--------------------------------------

if prod(size(dat)) <= 1
   data = dat;
   return
end

%--------------------------------------

n = cat( 2 , n , n );
n = n([1 2]);

%---------------------------------------------------
% Check for Deviation and Offset

off = abs(imag(v));

v   = sign(real(v));

off = min(max(off,0),1);

if off == 0
   off = 1 / prod(n);
end

%--------------------------------------
% Check for Window-Inputs

if isempty(w) 
   w = {};
elseif ~iscell(w)
   w = {w};
end

%--------------------------------------

try
  c1 = window(n(1),w{:});
  c2 = window(n(2),w{:});
catch
  error(lasterr);
end

if isequal(n,[1 1])
   data = dat;
   return
end

%---------------------------------------------------

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

i1 = ( 1 : s1 );
i2 = ( 1 : s2 );

%---------------------------------------------------
% Limits

m0 = NaN * zeros(s1,s2);
m1 = NaN * zeros(s1,s2);

for ii1 = 1 : n(1)
 for ii2 = 1 : n(2)

   m0 = min(m0,dat(i1+(ii1-1),i2+(ii2-1)));
   m1 = max(m1,dat(i1+(ii1-1),i2+(ii2-1)));

 end
end

m1 = m1 - m0;  % Intervall

if v >= 0
   if v == 0 
      m1 = m1/2;
   end
   m0 = m0 + m1;
end

m1 = m1 + ( m1 == 0 );  % Check for ZERO-Intervalls

%---------------------------------------------------
% Prepare Variables

quo = ~isnan(dat);

dat( find(~quo) ) = 0;

data = zeros(s1,s2);
quot = zeros(s1,s2);

%---------------------------------------------------
% Mean

for ii1 = 1 : n(1)

  for ii2 = 1 : n(2)

    cc = ( dat(i1+(ii1-1),i2+(ii2-1)) - m0 ) ./ m1;

    %---------------------------------------------------
    switch lower(m(1))
    %---------------------------------------------------

      %-------------------------------------------------
      case 'c' % Cosine

        cc = ( 1 + cos( pi * cc ) ) / 2;

      %-------------------------------------------------
      case 'g' % Gauss
 
        cc = exp( (-1) * (pi*cc).^2 / 2 );

        cc = ( cc - 1 ) / ( 1 - exp( (-1) * (pi).^2 / 2 ) ) + 1;

      %-------------------------------------------------
      otherwise % Triangle | Potenz

        if isinf(e)
           cc = ( abs(cc) < 1 );
        else
           if e == 0
              cc = double( cc == 0 );
           elseif e == 1
              cc = 1 - abs(cc);
           else    
              cc = 1 - abs(cc) .^ e;
           end
        end

    %---------------------------------------------------
    end
    %---------------------------------------------------

    cc = ( 1 - off ) * ( cc - 1 ) + 1;

    if ~( p == 1 )
        cc = cc .^ p;
    end

    cc = c1(ii1) * c2(ii2) * cc;

    data = data + cc .* dat(i1+(ii1-1),i2+(ii2-1));
    quot = quot + cc .* quo(i1+(ii1-1),i2+(ii2-1));

  end

end

quo = find(~quo(i1+n2(1),i2+n2(2)));

data = data ./ ( quot + ( quot == 0 ) );

data(quo) = NaN;

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [c,x] = window(n,m,p)

% WINDOW  returns an symetric, normalized window
%
% [C,X] = WINDOW( N , Mode , Potenz )
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
% Potenz  = Potenz for pre-normalized ( max == 1 ) window
%
%
% C       = normalized window, [ 1 by N ], sum(C) == 1 
%
% X       = X-Coordinates, [ 1 by N ], [ -1 .. 1 ]
%
%
%  defaults:  N = 1
%             Mode = 'linear'
%             Potenz = 1
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

