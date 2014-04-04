function [x,y,n,dy,y0,y1] = grpmean(x,y,d,f)

% GRPMEAN  Mean Value at specified Groups
%
% [ XI , YMean , N , DEV , YMin , YMax ] = GRPMEAN( X , Y , DIM )
%
%   X    GroupMarker, Number of Elements correspond
%          with Length of Y at DIM
%   Y    Matrice to find mean of Groups
%
%   DIM  Mean along Dimension 
%        (optional, default: longest Dimension)
%
% A nonzero imaginary part of DIM calculates the Median Value
%   of each Group.
%
%   XI      GroupValue
%   YMean   Mean of Group
%   N       Number of Averaged Elements
%   DEV     Standard Deviation of Group
%   YMin    Min of Group
%   YMax    Max of Group
%
% DEV is normalized by  N-1 (see STD), 
%  use GRPMEAN(X,Y,DIM,1) to normalize with N.
%
%---------------------------------------------------------------
%
% Find Groups of same Elements in X:
%
% [ XI , XI , N , ZEROS , XI , XI ] = GRPMEAN( X , DIM )
%
% [ XI , XI , N , ZEROS , XI , XI ] = GRPMEAN( X , X , DIM )
%
%---------------------------------------------------------------
% 
% See also: MEDMEAN, IND2GRP, GRP2IND, CUMSUM, CUMPROD, UNIQUED
%

Nin  = nargin;
Nout = nargout;

if Nin < 2
   error('Not enough InputArguments.');
end

if Nin < 3
   if ~( prod(size(x)) == 1 ) & ( prod(size(y)) == 1 )
       d = y
       y = x;
   else 
       d = [];
   end
end

if Nin < 4
   f = 0;
end

f = 1 - isequal(f,1);

is_std = ( Nout >= 4 );
is_min = ( Nout >= 5 );
is_max = ( Nout >= 6 );

%--------------------------------------------------------

n  = [];
dy = [];
y0 = [];
y1 = [];

if isempty(x) & isempty(y)
   return
end

%--------------------------------------------------------
% Get DIM

si = size(y);

if isempty(d)
   d = 0;
end

is_med = ~(imag(d)==0); % True for Median

d = real(d);

if d == 0
   [m,d] = max(si);
   if ~( m == prod( si + ( si == 0 ) ) )  % Not a Vector
       d = 1;
   end
else
   si = cat( 2 , si , ones(1,d-size(si,2)) );
end

%--------------------------------------------------------
% Permute X and Y to 1. Dimension at DIM

m = size(si,2);

p = cat( 2 , d , ( 1 : d-1 ) , ( d+1 : m ) );

s1 = si(d);
s2 = prod(si(p(2:m)));

y = reshape(permute(y,p),s1,s2);

x = x(:);

if ~( size(x,1) == size(y,1) )
    error(sprintf('Number of Elements of X must match Size of Y in Dimension %.0f.',d));
end

%--------------------------------------------------------
% Sort via X (and first via Y if Median, Min or Max)

if is_med | is_min | is_max
   [y,ii] = sort(y,1,'ascend');
   [n,jj] = sort(x(ii),1,'ascend');
    x     = n(:,1);
    n     = ones(s1,1) * ( 0 : (s2-1) ) * s1;
    y     = y(jj+n);
else
   [x,ii] = sort(x,1,'ascend');
   y      = y(ii,:);
end

%--------------------------------------------------------
% Find Groups of X

i0 = ~( diff(x,1,1) == 0 );
ng = sum(i0) + 1;

i0 = cat( 1 , 1 , find(i0)+1 , size(y,1)+1 );
lg = diff(i0,1,1);
i0 = i0(1:ng);

i1 = i0 + lg - 1;

j1 = i1( 1 : (ng-1) );

z2 = zeros(1,s2);

x = x(i0);

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
% Prepare for STD

if is_std
   dy = y;
end

%********************************************************
if is_med
%********************************************************
% Median of Groups

k2 = mod(n,2);          % Even or Odd Number of Elements

k1 = ( n + k2 ) / 2 + n0;

k1 = k1 + ( i0(:,ones(1,s2)) - 1 );
k1 = k1 + ones(ng,1) * ( 0 : (s2-1) ) * s1;

k2 = k1 + 1 - k2 - n0;  % k2 = k1 + 1  if  Even; k2 == k1 if Odd

y = ( y(k1) + y(k2) ) / 2;

%********************************************************
else
%********************************************************
% Mean of Groups in Y

if isn
   y(nn) = 0;
end

y = cumsum(y,1);

y = y(i1,:) - cat( 1 , z2 , y(j1,:) );

y = y ./ (n+n0);

%********************************************************
end
%********************************************************

%--------------------------------------------------------
if is_std

   ii = zeros(s1,1);
   ii(i0) = 1;
   ii = cumsum(ii,1);

   dy = ( dy - y(ii,:) );   % Remove Mean or Median

   if isn
      dy(nn) = 0;
   end

   dy = cumsum( conj(dy) .* dy , 1 );
   dy = dy(i1,:) - cat( 1 , z2 , dy(j1,:) );

   n1 = ( n == 1 );

   dy = sqrt( dy ./ ( n+2*n0+n1 - f ) );

   if any(n1(:))
          n1  = find(n1);
       dy(n1) = 0;
   end

end

%--------------------------------------------------------

if any(n0(:))
     n0  = find(n0);
   y(n0) = NaN;
   if is_min, y0(n0) = NaN; end
   if is_max, y1(n0) = NaN; end
   if is_std, dy(n0) = NaN; end
end

%--------------------------------------------------------
% Permute back

sx       = ones(1,m);
sx(p(1)) = ng;
si(p(1)) = ng;

sx = sx(p);
si = si(p);

x = reshape(x,sx);
y = reshape(y,si);
n = reshape(n,si);

if is_min,  y0 = reshape(y0,si); end
if is_max,  y1 = reshape(y1,si); end
if is_std,  dy = reshape(dy,si); end

p = cat( 2 , ( 1 : d-1 ) + 1 , 1 , ( d+1 : m ) );

x = permute( x , p );
y = permute( y , p );
n = permute( n , p );

if is_min,  y0 = permute(y0,p); end
if is_max,  y1 = permute(y1,p); end
if is_std,  dy = permute(dy,p); end
