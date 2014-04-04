function  [y,nn,dy,ym,y0,y1] = medmean(y,n,d,f)

% MEDMEAN  1-Dimensional Median Filter
%
% [ Y , NN ] = MEDMEAN( X , N , DIM )
%
%   X    Matrice to filter
%
%   N    Length of Window
%
%        A NEGative imaginary part of N uses a LEFTside window. 
%        (looking backward, aligned at begin to origin)
%
%        A POSitive imaginary part of N uses a RIGHTside window.
%        (looking foreward, aligned at end to origin)
%
%   DIM  Median along Dimension 
%        (optional, default: longest Dimension)
%
%   Y    Median-filtred X
%
%   NN   Number of non-NaN-Elements used for Median
%
%------------------------------------------------------
%
% [ Y , NN , DEV , XMean , XMin , XMax ] = MEDMEAN( ... )
%
%   DEV     Standard Deviation in Window
%   XMean   Mean in Window
%   XMin    Min in Window
%   XMax    Max in Window
%
% DEV is normalized by  N-1 (see STD), 
%  use MEDMEAN(X,N,DIM,1) to normalize with N.
%
% See also: MEANIND1, GRPMEAN, MEDFILT1 (Signal Processing Toolbox)
%

Nin  = nargin;
Nout = nargout;

if Nin < 2
   error('Not enough InputArguments.');
end

if Nin < 3
   d = [];
end

if Nin < 4
   f = 0;
end

f = 1 - isequal(f,1); % !!!

is_std  = ( Nout >= 3 );
is_mean = ( Nout >= 4 );
is_min  = ( Nout >= 5 );
is_max  = ( Nout >= 6 );

%--------------------------------------------------------

ok = ( isnumeric(n) & ( prod(size(n)) == 1 ) );
if ok
   sg = sign(imag(n));
    n = real(n);
   ok = ( ( n > 0 ) & ( mod(n,1) == 0 ) );
   ok = ( ok & ( ~( sg == 0 ) | ( mod(n,2) == 1 ) ) );
end

if ~ok
    error('N must be an Integer and oddinary if symetric Window')
end

nn = [];
dy = [];
ym = [];
y0 = [];
y1 = [];

if isempty(y)
   return
end

%--------------------------------------------------------
% Get DIM

sz = size(y);

if isempty(d)
   d = 0;
end

if d == 0
   [m,d] = max(sz);
   if ~( m == prod( sz + ( sz == 0 ) ) )  % Not a Vector
       d = 1;
   end
else
   sz = cat( 2 , sz , ones(1,d-size(sz,2)) );
end

%--------------------------------------------------------
% Permute Y to 1. Dimension at DIM

m = size(sz,2);

p = cat( 2 , d , ( 1 : d-1 ) , ( d+1 : m ) );

s1 = sz(d);
s2 = prod(sz(p(2:m)));

y = reshape(permute(y,p),s1,s2);

%--------------------------------------------------------
% Sort first via Y

[z,si] = sort(y,1,'ascend');

s0 = ( sg == 0 );  % Symetric
fw = ( sg >= 0 );  % Looking ForeWard
bw = ( sg <= 0 );  % Looking BackWard

n2 = ( n - s0 ) / ( 1 + s0 );

nn = zeros(s1,s2);

if is_std,  dy = nn; end
if is_mean, ym = nn; end
if is_min,  y0 = nn; end
if is_max,  y1 = nn; end

for ii = 1 : min(n,s1)

   x = ceil( ( ( 1 : s1 )' + bw*n2 - (ii-fw) ) / n );

   jj = ( ii : n : s1 );

% fprintf(1,'%2.0f/%.0f: %s  %2.0f..%.0f\n',ii,n,sprintf('%2.0f ',x(1:10)),jj([1 end])); 

   [y(jj,:),nn(jj,:),dx,xm,x0,x1] = grpmed(x(si),z,size(jj,2),s0*n2,f,...
                                           is_std,is_mean,is_min,is_max);

   if is_std,  dy(jj,:) = dx; end
   if is_mean, ym(jj,:) = xm; end
   if is_min,  y0(jj,:) = x0; end
   if is_max,  y1(jj,:) = x1; end

end

%--------------------------------------------------------
% Permute back

%%% sz(p(1)) = ng;

sz = sz(p);

y = reshape(y,sz);
nn = reshape(nn,sz);

if is_std,  dy = reshape(dy,sz); end
if is_mean, ym = reshape(ym,sz); end
if is_min,  y0 = reshape(y0,sz); end
if is_max,  y1 = reshape(y1,sz); end

p = cat( 2 , ( 1 : d-1 ) + 1 , 1 , ( d+1 : m ) );

y = permute( y , p );
n = permute( n , p );

if is_std,  dy = permute(dy,p); end
if is_mean, ym = permute(ym,p); end
if is_min,  y0 = permute(y0,p); end
if is_max,  y1 = permute(y1,p); end

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [y,n,dy,ym,y0,y1] = grpmed(x,y,nx,nl,flg,...
                                    is_std,is_mean,is_min,is_max);

dy = [];
ym = [];
y0 = [];
y1 = [];

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
% Min and Max

if is_min | is_max
   ii = i0(:,ones(1,s2));
   ii = ii + ones(ng,1) * ( 0 : (s2-1) ) * s1;
   if is_min
      y0 = y(ii);
   end
   if is_max
      ii = ii + n - 1 + n0;
      y1 = y(ii);
   end
end

%--------------------------------------------------------
% Prepare for MEAN & STD

if is_mean | is_std

   ym = y; 
   if isn
      ym(nn) = 0;
   end

   ym = cumsum(ym,1);

   ym = ym(i1,:) - cat( 1 , z2 , ym(j1,:) );

   ym = ym ./ (n+n0);

   if is_std  

      ii = zeros(s1,1);
      ii(i0) = 1;
      ii = cumsum(ii,1);

      ii(find(ii==0)) = 1;

      dy = ( y - ym(ii,:) );   % Remove Mean or Median

      if isn
         dy(nn) = 0;
      end

      dy = cumsum( conj(dy) .* dy , 1 );
      dy = dy(i1,:) - cat( 1 , z2 , dy(j1,:) );

      n1 = ( n == 1 );

      dy = sqrt( dy ./ ( n+2*n0+n1 - flg ) );

      if any(n1(:))
             n1  = find(n1);
          dy(n1) = 0;
      end

   end

end

%--------------------------------------------------------
% Median of Groups

k2 = mod(n,2);          % Even or Odd Number of Elements

k1 = ( n + k2 ) / 2 + n0;

k1 = k1 + ( i0(:,ones(1,s2)) - 1 );
k1 = k1 + ones(ng,1) * ( 0 : (s2-1) ) * s1;

k2 = k1 + 1 - k2 - n0;  % k2 = k1 + 1  if  Even; k2 == k1 if Odd

y = ( y(k1) + y(k2) ) / 2;

%--------------------------------------------------------

if any(n0(:))
     n0  = find(n0);
   y(n0) = NaN;
   if is_std,  dy(n0) = NaN; end
   if is_mean, ym(n0) = NaN; end
   if is_min,  y0(n0) = NaN; end
   if is_max,  y1(n0) = NaN; end
end


