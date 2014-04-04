function [y,d] = chk_lon(x,l);

% CHK_LON  Check if Longitude in correct Period for Limits
%
% [ LonC , CHK ] = CHK_LON( Lon , Limit )
%
%  Lon   LongitudeData to check
%  Lim   Limit: [ LonMin  LonMax ]
%        AxesHandle or Children to get XLimits
%
%  LonC  Transformed Lon which match Limits
%        if neccessary multiple Periods are appended
%         after the last nonsingelton dimension
%
%  CHK   Appended Dimension: Number + i*Size
%

y = [];
d = 0;

if prod(size(l)) == 1
   ok = ishandle(l);
   if ok
      if ~strcmp( get(l,'type'),'axes' ); 
          l = get(l,'parent');
      end
      ok = strcmp( get(l,'type'),'axes' );
      if ok
         l = get(l,'xlim');
      end
   end
   if ~ok 
       error('Single Input must be an Handle of Axes or Children of it.')
   end
end


while l(2) < l(1)
      l(2) = l(2) + 360;
end

%--------------------------------------------
% Check X

if ~all(isfinite(x(:)))
   y = x;
   return
end

% Extension of X

e = [ max(x(:))  min(x(:)) ];

%--------------------------------------------
% Check Extension !!!

dx = ( e(1) - e(2) );

if dx > 180
   dx = 360 * ceil(dx/360);
   x = x + dx * ( ( e(1) - x ) > 180 );
   e = [ max(x(:))  min(x(:)) ];
end

%------------------------------------------------------------------
% Check from Max of X touch Min of XLim from left
%         to Min of X touch Max of XLim from right
%------------------------------------------------------------------
% max(x) + k1 * 360 <= min(l)
%   e(1) + k1 * 360 <= l(1)
%     k1 <= ( l(1) - e(1) ) / 360 
% 
%  kk(1) =  ceil( ( l(1) - e(1) ) / 360 )
%  kk(1) =  kk(1) + ( l(1) == e(1)+kk(1)*360 )
%------------------------------------------------------------------
% min(x) + k2 * 360 >= max(l)
%   e(2) + k2 * 360 >= l(2)
%                k2 >= ( l(2) - e(2) ) / 360 
% 
%  kk(2) = floor( ( l(2) - e(2) ) / 360 )
%  kk(2) = kk(2) - ( l(2) == e(2)+kk(2)*360 )
%------------------------------------------------------------------

op = [ -1  1 ];

k = floor( op .* ( l - e ) / 360 );
k = k - ( l == e+k*360 );
k = k .* op;

% Complete Outside
if ( k(1) > k(2) ) 
   y = x;
end


m =   k(2) - k(1) + 1;
k = ( k(1) : k(2) );

s = size(x);
n = size(s,2);   % NDims

while ( s(n) == 1 ) & ( n > 1 )
        n = n - 1;
end

d = n+1 + m*i;

s = cat( 2 , s(1:n) , m );  %  Append Dimension
   
y = zeros(s);

s = prod(s(1:n));


for ii = 1 : m

    jj  = ( 1 : s ) + (ii-1) * s;

  y(jj) = x + k(ii)*360;

end

