function [data,dev] = meanval1(dat,n,w,m,v,p);

% MEANVAL1  Running Value-Mean along 1. Dimension
%
% [Y,STDEV] = MEANVAL1( X , N , W , Mode , V , P );
%
% N = Length of window, oddinary  Integer
%
% W = Inputs for N-based window: { Mode  Potenz }
%
%    Mode = Mode of window: {'linear'} | 'triangle' | 'cosine' | ...
%                            'gauss'   | 'binomial'
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
% A POSitive IMAGinary part of N will use a  LEFTside window only. 
%  (looking backward, aligned at begin to origin)
%
% A NEGative IMAGinary part of N will use a RIGHTside window only.
%  (looking foreward, aligned at end to origin)
%
% If the imaginary part of N is ZERO, the values at Begin and End 
%  will aligned to the origin, a negative Value of N prevents this.
%
% NaN's in X make no Problems.
%
%------------------------------------------------------
% Example
%
%  p = 50;   % Length of Period
%
%  x = ( 0 : 5.5*p );
%  r = randn(size(x)) / (2*pi);      % Noise
%  y = (1+r) .* cos( 2*pi * x./(p*(1+r/5)) );
%
%  n = 2*p+1; % SmoothIntervall
%  m = 'c'; % Mode
%
%  figure, hold on
%
%  plot( x , y , 'k' )
%
%  plot( x , meanind1(y, n  ,m) , 'b' );
%  plot( x , meanind1(y,-n  ,m) , 'c' );
%  plot( x , meanind1(y, n+i,m) , 'g' );
%  plot( x , meanind1(y, n-i,m) , 'r' );
%
%


Nin  = nargin;
Nout = nargout;

%***************************************************************
% Check Inputs

msg = cell(0,1);

%---------------------------------------------------
% Check N

if Nin < 2
   n = 1;
else
  ok = ( isnumeric(n)  &  ( prod(size(n)) == 1 ) );
  if ok
     ok = ( mod(real(n),2) == 1 );
  end
  if ~ok
      msg = cat( 1 , msg , {'N must be an odd Integer.'} );
  end
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

if isempty(dat)
   data = dat;
   dev  = dat;
   return
end

%---------------------------------------------------
% Check N

in = sign(imag(n));
 n = real(n);

is_half = ~( in == 0 );

is_neg = ( ( n < 0 ) & ~is_half );

 n = abs(n);

if abs(n) == 1
   data = dat;
    dev = zeros(size(dat));
   return
end

%---------------------------------------------------
% Check for Deviation and Offset

off = abs(imag(v));

v   = sign(real(v));

off = min(max(off,0),1);

if off == 0
   off = 1/n;
end

%--------------------------------------
% Check for Window-Inputs

if isempty(w) 
   w = {};
elseif ~iscell(w)
   w = {w};
end

%--------------------------------------
% Get Window, expand if imag(N)  

n = n + ( n - 1 ) * is_half;

try
  c = window( (1-2*is_half)*n , w{:} );
catch
  error(lasterr);
end

if in == -1
   c = c( n : -1 : 1 );   % RightSide Window
end

%---------------------------------------------------
% Permute and Reshape to 2D

si = size(dat);
ps = prod(si);
ns = size(si,2);

if ( ps == 1 ) | ( ( ps == 2 ) & ~is_neg & ~is_half )
   data = dat;
    dev = zeros(si);
   return
end

dim  = sum(cumprod(double(si==1))) + 1;
prm = ( dim > 1 );
if prm
   perm = cat( 2 , ( dim : ns ) , ( 1 : dim-1) );
    dat = permute( dat , perm );
     si = si(perm);
end
     
rsh = ~( ps == si(1)*si(2) );
if rsh
   dat = reshape(dat,si(1),ps/si(1));
end

%---------------------------------------------------
% Expand dat with -n2  ..  +n2 - Elements

n2  = ( n - 1 ) / 2; 

s1 = size(dat,1);
s2 = size(dat,2);

if is_neg | is_half

   dat = cat( 1 , NaN*zeros(n2,s2) , dat , NaN*zeros(n2,s2) );

else

    s0 = s1;
    n0 = n2;

    while s0-s1 < 2*n2

       n1  = min(n0,s1-1);
       o1  = ones(n1,1);

        ds = ( s0 - s1 ) / 2;

       i01 = [ 1  s1 ] + ds;

       k1  = i01(1) + ds + ( n1 : -1 :  1 );
       k2  = i01(2) - ds - (  1 :  1 : n1 );

       dat = cat( 1 , 2*dat(i01(1)*o1,:)-dat(k1,:) , dat , ...
                      2*dat(i01(2)*o1,:)-dat(k2,:)             );

       s0 = s0 + 2*n1;
       n0 = n0 -   n1;

    end

end 

nind = ( 1+n2*(in==-1) : n-n2*(in==1) );

 ind  = ( 1 : s1 )';

%---------------------------------------------------
% Limits

m0 = NaN*zeros(s1,s2);
m1 = NaN*zeros(s1,s2);

for ii = nind
    m0 = min(m0,dat(ind+(ii-1),:));
    m1 = max(m1,dat(ind+(ii-1),:));
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
 
dat(find(~quo)) = 0;

data = zeros(s1,s2);
quot = zeros(s1,s2);

if Nout == 2
  nn = zeros(s1,s2);
end

%---------------------------------------------------
% Mean

for ii = nind;

    jj = ind + ( ii - 1 );

    cc = ( dat(jj,:) - m0 ) ./ m1;

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

    cc = c(ii) * cc;

    data = data + cc .* dat(jj,:);
    quot = quot + cc .* quo(jj,:);

    if Nout == 2
       nn = nn + quo(jj,:);
    end

end

data = data ./ ( quot + ( quot == 0 ) );

quo = find(~quo(ind+n2,:));

%------------------------------------------
% Set BorderValues

if ~is_neg

    ind = [ 1  s1 ];
  
    if is_half
       ind = ind( 1 + ( in == -1 ) );
    end

    data(ind,:) = dat(ind+n2,:);

end

%------------------------------------------
if Nout == 2

 dev = zeros(s1,s2);

 for ii = 1 : n
     dev = dev + (dat(ind+(ii-1),:)-data).^2;
 end

 qadd  = 2*(nn==0)+(nn==1);

 dev = sqrt( dev ./ (nn+qadd-1) );

 dev = dev .* (~(qadd==1));

 dev(quo) = NaN;

   if rsh
      dev = reshape(dev,si);
   end
   if prm
      dev = ipermute(dev,perm);
   end

end
%------------------------------------------

data(quo) = NaN;

if rsh
   data = reshape(data,si);
end
if prm
   data = ipermute(data,perm);
end

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

