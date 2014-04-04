function ok = online(xy,x,y,d)

% ONLINE  True for Points near a Curve
%
%  ok = ONLINE( XY , XP , YP , Limit )
%
%  XY: 2-Column Matrice, defines the Curve in X and Y
%
%  XP: Vector or Matrice with X-Coordinates of Points
%  YP: Vector or Matrice with Y-Coordinates of Points
%
%  Limit:  Maximum allowed Distance between the Points and the Curve
%
%  ok: Matrice with same Size of XP and YP, Points inside the Limits
%       get the Value ONE, the other ZERO
%
% see also: INPOLYGON
%
% example:
% 
%    xy = linspace(0,2*pi,20)';    % SIN-Curve
%    xy(:,2) = 2 * sin(xy(:,1));
%
%    xp = ( -0.5 : 0.1 : 7   );
%    yp = ( -2.5 : 0.1 : 2.5 );
%    
%    [xp,yp] = meshgrid( xp , yp );
%
%    lim = 0.4;   % Look for Points with Distance from Curve max. 0.4
%
%    ok = online(xy,xp,yp,0.4);
%
%    ii = find(ok);
%
%    figure, axes('dataaspectratio',[1 1 1],'nextplot','add')
%
%     plot(xp,yp,'k.','markersize',4);   % all Points black
%     plot(xp(ii),yp(ii),'g*');          % good Points green "on" SIN-Curve
%     plot(xy(:,1),xy(:,2),'r-+');       % SIN-Curve red with cross
%
%    % Make Circles around the CurvePoints
%
%     p = linspace(0,2*pi,100);
%     x = lim * sin(p);  y = lim * cos(p);
%     for ii = 1 : size(xy,1)
%         plot(x+xy(ii,1),y+xy(ii,2),'r-');
%     end
%

%-------------------------------------------------
% Check XY, Make ColumnVector

si = size(xy);
ii = find( si == 2 );
ns = prod(size(si));

if isempty(ii) | ( prod(si) < 2 )
   error('XY must define an 2-dimensional Line.');
end

ii = ii(1) + ( 2 - ii(1) ) * any(ii==2);

perm = cat( 2 , ( 1 : ii-1 ) , ( ii+1 : ns ) , ii );

si = si(perm);
xy = permute( xy , perm );

xy = reshape( xy , prod(si(1:ns-1)) , si(ns) );

%-------------------------------------------------
% Check X and Y

if isempty(x) & isempty(y)
   ok = zeros(0,0);
   return
end

sx = size(x);
sy = size(y);

if ~isequal(sx,sy)
   error('X and Y must have the same Size')
end

%-------------------------------------------------
% Make ROW-Vector

x = permute( x(:) , [2 1] );
y = permute( y(:) , [2 1] );

n2 = prod(sx);

%-------------------------------------------------

n1 = size(xy,1) - 1 ;
n2 = size(x,2);

jj = ( 1 : n1 );

o1 = ones(n1,1);
o2 = ones(1,n2);

%-------------------------------------------------

if n1 == 0

   d0 = 0;

else

  dx = xy(jj+0,1) - xy(jj+1,1);
  dy = xy(jj+0,2) - xy(jj+1,2);

  % Normalize: [ ca , sa ]

  d0 = sqrt( dx.^2 + dy.^2 );

end

%-------------------------------------------------
% Check for ZERO-Lines (same Points)

%---------------------------------------------
% Single Point

if     all( d0 == 0 )

  d1 = sqrt( ( x(1,:) - xy(1,1) ).^2 + ...
             ( y(1,:) - xy(1,2) ).^2       );

  ok = ( d1 <= d );

  ok = reshape(ok,sx);

  return

%---------------------------------------------
elseif any( d0 == 0 )

  d0(find(d0==0)) = NaN;
 
end

%-------------------------------------------------

ca = -dy ./ d0;
sa = +dx ./ d0;

% Hesse-Form:  ca*x + sa*y - p = 0

p = xy(jj,1).*ca + xy(jj,2).*sa;

%-------------------------------------------------------
% Check Distance of points to line-segments

ok = ( abs( ca * x + sa * y - p(:,o2) ) <= d );

%-------------------------------------------------------
% Check of points inside line-segments

% Kosinussatz: a^2 = b^2 + c^2 -2*b*c*cos(alpha)
%
%    cos(alpha) = ( b^2 + c^2 - a^2 ) / ( 2 * b * c )
%
% Punkt liegt NICHT zw. P1 und P2 wenn einer der Winkel > 90°
%  also:   ( cos(alpha) < 0 )  ==  ( a^2 > ( b^2 + c^2 ) ) 
%
%  also:   bad = ( ( ( a^2 > ( b^2 + c^2 ) ) & ( b > d ) ) | ...
%                  ( ( b^2 > ( a^2 + c^2 ) ) & ( a > d ) )        );
%
%           ok = ( ~( ( a^2 > ( b^2 + c^2 ) ) & ( b > d ) ) & ...
%                  ~( ( b^2 > ( a^2 + c^2 ) ) & ( a > d ) )        );
%           ok = ( ( ~( a^2 > ( b^2 + c^2 ) ) | ~( b > d ) ) & ...
%                  ( ~( b^2 > ( a^2 + c^2 ) ) | ~( a > d ) )        );
%
%           ok = ( ( ( a^2 <= ( b^2 + c^2 ) ) | ( b <= d ) ) & ...
%                  ( ( b^2 <= ( a^2 + c^2 ) ) | ( a <= d ) )        );
%
%  mit  c == d0
%

% Distance to Points of LineSegments (a,b)

d1 = sqrt( ( x(o1,:) - xy(jj+0,1*o2) ).^2 + ...
           ( y(o1,:) - xy(jj+0,2*o2) ).^2       );

d2 = sqrt( ( x(o1,:) - xy(jj+1,1*o2) ).^2 + ...
           ( y(o1,:) - xy(jj+1,2*o2) ).^2       );


ok = ( ok & ...
       ( ( ( d2.^2 <= ( d1.^2 + d0(:,o2).^2 ) ) | ( d1 <= d ) ) & ...
         ( ( d1.^2 <= ( d2.^2 + d0(:,o2).^2 ) ) | ( d2 <= d ) )        ) );
               

ok = ( sum(ok,1) > 0 );  % Ok if any LineSegment is ok.


ok = reshape(ok,sx);
