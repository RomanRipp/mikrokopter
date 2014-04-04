function [data,nn,dv,c] = meanind1(dat,n,m,p);

% MEANIND1  1-Dimensional Running Mean
%
% Y = MEANIND1( X , N , Mode , P );
%
% X    = Vector or Matrice to filter along 1. Dimension
%
% N    = Length of window, oddinary Integer if ZERO Imaginary Part
%
%        A POSitive imaginary part of N uses a LEFTside window. 
%        (looking backward, aligned at begin to origin)
%
%        A NEGative imaginary part of N uses a RIGHTside window.
%        (looking foreward, aligned at end to origin)
%
%        A ZERO imaginary part of N uses a symetric window,
%        the values at Begin and End are aligned to the origin, 
%        a negative Value of N prevents this.
%        (doesn't affect median-Mode)
%
% Mode = Mode of window:  'cosine'   | 'gauss'    | ...
%                         'triangle' | 'binomial' | ...
%                         'linear'   | 'median'   | ...
%                         'min'      | 'max'
%
% Mode = positive Numeric for Power of Window in Range ( -1 .. 1 )
%          Mode == 0    ==>  Impuls
%          Mode == 1    ==>  Triangle
%          Mode == Inf  ==>  Linear (constant)
%
% Mode = Vector of Weights for window with Length N
%
% P    = Power for pre-normalized window
%
%  defaults:  N    = 1         (No mean)
%             Mode = 'linear'  (constant)
%             P    = 1
%
% For Matrices X MEANIND1 works along the 1. Dimension,
%  NaN's in X make no Problems.
%
%------------------------------------------------------
%
% [ Y , NN , STdv , Coeff ] = MEANIND1( X , N , ... )
%
% returns the Number NN of averaged Elements, the Standarddeviation 
%  per Element and the Coefficients of the selected Window.
%
%------------------------------------------------------
%
% See also: MEANIND2, MFILTER, MEDMEAN, WINDOW, NOISE
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

is_nr = ( Nout > 1 );
is_dv = ( Nout > 2 );

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
    dv =  zeros(si);
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

%---------------------------------------------------
% Check N

msg = {};

sn = size(n); pn = prod(sn);

ok = ( isnumeric(n)  &  ( pn == 1 ) );
if ok
   % Single Number
   sg = sign(imag(n));
   rn = real(n);
   ok = ( ( mod(rn,1) == 0 ) & ( ~( sg == 0 ) | ( mod(rn,2) == 1 ) ) );
   if ~ok
       msg = cat(1,msg,{'N must be an Integer and oddinary if symetric Window.'});
   end
elseif isnumeric(n) & ( Nin == 2 )
   % Check for Vector with odd. Number of Elements if NO further Inputs
   ok = ( ( pn == max(sn) ) & ( mod(pn,2) == 1 ) );
   if ok
      m = n;
      n = pn;
   else
       msg = cat(1,msg,{'A Vector for N must have an oddinary Number of Elements.'});
   end
else
    msg = cat(1,msg,{'N must be an Integer or a Vector.'});
end

%---------------------------------------------------
% Check Mode

sm = size(m); pm = prod(sm);

ok = ( ( ischar(m) & ( pm == sm(2) ) ) | isnumeric(m) );
ok = ( ok & ~isempty(m) );

if ~ok
    msg = cat(1,msg,{'Mode must be a String or a Vector.'});
end

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Invalid Inputs.\n%s',msg));
end

%---------------------------------------------------
% Get Window, expand if imag(N)  

is_med = ischar(m);
if is_med
   is_med = strcmp(lower(m(1)),'m');
end

sg = sign(imag(n));  % Left- or RightSide Window
 n = real(n);

is_half = ~( sg == 0 );

is_neg = ( ( n < 0 ) & ~is_half );

 n = abs(n);

n0 = n;

n = n + ( n - 1 ) * is_half;

if isnumeric(m) & ( pm > 1 )
   % Check for Vector
   if ~( ( pm == n0 ) | ( is_half & ( pm == n ) ) )
       error('Length of Vector for Mode must correspond with N.')
   end
   try
      c = window( n , 'l' , p );
   catch
      error(sprintf('Error create Window.\n%s',lasterr));
   end
   c = m(:)';
   m = max(m);
   c = ( c / (m+(m==0)) ) .^ p;
   m = sum(c);
   c = c / (m+(m==0));
   if is_half & ( pm == n0 )
      c = cat( 2 , c , zeros(1,n0-1) );
   end
elseif is_med
   c = ones(1,n);
else
   % Get Window
   try
      c = window( (1-2*is_half)*n , m , p );
   catch
      error(sprintf('Error create Window.\n%s',lasterr));
   end
end

if is_half
   c((n0+1):n) = 0;
   if sg == -1
      c = c( n : -1 : 1 );   % RightSide Window
   end
end

%--------------------------------------

if n == 1
   data = dat;
   if Nout > 1
      dv =  zeros(si);
       nn = double(~isnan(dat));
   end
   return
end

%---------------------------------------------------
% Permute and Reshape to 2D

ns = size(si,2);

if ( ( ps == 2 ) & ~is_neg & ~is_half )
   data = dat;
   if Nout > 1
    dv =  zeros(si);
     nn = ~isnan(dat);
   end
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

s1 = size(dat,1);
s2 = size(dat,2);

%***************************************************
if is_med  % MedianFilter
%***************************************************

   n = n0;

   md = lower([ m 'e' ]);

   md = ( md(2) == 'a' )  - ( md(2) == 'i' );

   %--------------------------------------------------------
   % Sort first 

   [dat,si] = sort(dat,1,'ascend');

   s0 = ( sg == 0 );  % Symetric
   fw = ( sg <= 0 );  % Looking ForeWard
   bw = ( sg >= 0 );  % Looking BackWard

   n2 = ( n - s0 ) / ( 1 + s0 );

   data = zeros(s1,s2);

   if is_nr, nn = data; end
   if is_dv, dv = data; end

   for ii = 1 : min(n,s1)

       x = ceil( ( ( 1 : s1 )' + bw*n2 - (ii-fw) ) / n );

       jj = ( ii : n : s1 );

% fprintf(1,'%2.0f/%.0f: %s  %2.0f..%.0f\n',ii,n,sprintf('%2.0f ',x(1:10)),jj([1 end])); 

       [data(jj,:),nr,dy] = grpmed(x(si),dat,size(jj,2),s0*n2,md,is_dv);

       if is_nr,  nn(jj,:) = nr; end
       if is_dv,  dv(jj,:) = dy; end

   end

%***************************************************
else
%***************************************************

   %---------------------------------------------------
   % Expand dat with -n2  ..  +n2 - Elements

   n2  = ( n - 1 ) / 2; 

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

   %---------------------------------------------------
   % Prepare Variables

   quo = ~isnan(dat);
 
   dat(find(~quo)) = 0;

   ind  = ( 1 : s1 )';
   data = zeros(s1,s2);
   quot = zeros(s1,s2);

   if is_nr
      nn = zeros(s1,s2);
   end

   %---------------------------------------------------
   % Mean

   for ii = 1+n2*(sg==-1) : n-n2*(sg==1)
       data = data + c(ii)*dat(ind+(ii-1),:);
       quot = quot + c(ii)*quo(ind+(ii-1),:);
       if is_nr
          nn = nn + quo(ind+(ii-1),:);
       end
   end

   data = data ./ ( quot + ( quot == 0 ) );

   quo = find(~quo(ind+n2,:));

   %------------------------------------------
   % Set BorderValues

   if ~is_neg

       ii = [ 1  s1 ];
  
       if is_half
          ii = ii( 1 + ( sg == -1 ) );
       end

       data(ii,:) = dat(ii+n2,:);

   end

   %------------------------------------------
   if is_dv

      dv = zeros(s1,s2);

      for ii = 1 : n
          dv = dv + (dat(ind+(ii-1),:)-data).^2;
      end

      qadd  = 2*(nn==0)+(nn==1);

      dv = sqrt( dv ./ (nn+qadd-1) );

      dv = dv .* (~(qadd==1));

      dv(quo) = NaN;

   end
   %------------------------------------------

   data(quo) = set_nan;

   if is_nr
      nn(quo) = 0;
   end

%***************************************************
end
%***************************************************

if rsh
   data = reshape(data,si);
   if is_nr, nn = reshape(nn,si); end
   if is_dv, dv = reshape(dv,si); end
end

if prm
   data = ipermute(data,perm);
   if is_nr, nn = ipermute(nn,perm); end
   if is_dv, dv = ipermute(dv,perm); end
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


%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [y,n,dy] = grpmed(x,y,nx,nl,md,is_std)

dy = [];

s1 = size(y,1);
s2 = size(y,2);

%--------------------------------------------------------
% Sort via X

[n,jj] = sort(x,1,'ascend');
 x     = n(:,1);
 n     = ones(s1,1) * ( 0 : (s2-1) ) * s1;
 y     = y(jj+n);


%--------------------------------------------------------
% Find Groups of X

i0 = ~( diff(x,1,1) == 0 );
ng = sum(i0) + 1;

i0 = cat( 1 , 1 , find(i0)+1 , size(y,1)+1 );
lg = diff(i0,1,1);
i0 = i0(1:ng);

x = x(i0);

%--------------------------------------------------------
% Check for good Groups

ii = ( ( 1 <= x ) & ( x <= nx ) & ( lg >= nl ) );

if ~all(ii)

    ng = sum(ii);

    ii = find(ii);

    i0 = i0(ii);
    lg = lg(ii);
    x  =  x(ii);

    %-------------------------------------------------
    %%% ii = grp2ind(i0,lg)

    s1 = sum(lg);

    ii = ones(s1,1);
    jj = cumsum( cat(1,1,lg) , 1 );

    ii(jj(1:ng)) = i0;


    if ng > 1
       kk = ( 1 : (ng-1) );
       ii(jj(kk+1)) = ii(jj(kk+1)) - (i0(kk,1)+(lg(kk)-1));
    end

    ii = cumsum(ii,1);

    %-------------------------------------------------

     y = y(ii,:);

    i0 = jj(1:ng);  % jj = cumsum( cat(1,1,lg) , 1 );

end

i1 = i0 + lg - 1;

j1 = i1( 1 : (ng-1) );

z2 = zeros(1,s2);

%--------------------------------------------------------
% Check for NaN's in Y

n = isnan(y);

isn = any(n(:));
nn  = [];
if isn
   nn = find(n);
end

n = ~n;  % True for NOT-NaN !!!

%--------------------------------------------------------
% Number of Elements per Group

n = cumsum(n,1);

n = n(i1,:) - cat( 1 , z2 , n(j1,:) );

n0 = ( n == 0 );

%--------------------------------------------------------
% Prepare for MEAN & STD

if is_std

   dy = y; 
   if isn
      dy(nn) = 0;
   end

   dy = cumsum(dy,1);

   dy = dy(i1,:) - cat( 1 , z2 , dy(j1,:) );

   dy = dy ./ (n+n0);  % MeanValue

      ii = zeros(s1,1);
      ii(i0) = 1;
      ii = cumsum(ii,1);

      ii(find(ii==0)) = 1;

      dy = ( y - dy(ii,:) );   % Remove Mean or Median

      if isn
         dy(nn) = 0;
      end

      dy = cumsum( conj(dy) .* dy , 1 );
      dy = dy(i1,:) - cat( 1 , z2 , dy(j1,:) );

      n1 = ( n == 1 );

      dy = sqrt( dy ./ ( n+2*n0+n1 - 1 ) );

      if any(n1(:))
             n1  = find(n1);
          dy(n1) = 0;
      end

end

%--------------------------------------------------------
if ( md == 0 )
%--------------------------------------------------------
% Median of Groups

   k2 = mod(n,2);          % Even or Odd Number of Elements

   k1 = ( n + k2 ) / 2 + n0;

   k1 = k1 + ( i0(:,ones(1,s2)) - 1 );
   k1 = k1 + ones(ng,1) * ( 0 : (s2-1) ) * s1;

   k2 = k1 + 1 - k2 - n0;  % k2 = k1 + 1  if  Even; k2 == k1 if Odd

   y = ( y(k1) + y(k2) ) / 2;

%--------------------------------------------------------
else
%--------------------------------------------------------

   ii = i0(:,ones(1,s2));
   ii = ii + ones(ng,1) * ( 0 : (s2-1) ) * s1;

   if md > 0
      ii = ii + n - 1 + n0;
   end

   y = y(ii);

%--------------------------------------------------------
end
%--------------------------------------------------------

if any(n0(:))
     n0  = find(n0);
   y(n0) = NaN;
   if is_std,  dy(n0) = NaN; end
end


