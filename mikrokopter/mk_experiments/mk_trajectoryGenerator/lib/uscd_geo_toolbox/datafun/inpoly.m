function in = inpoly(x,y,xv,yv)

% INPOLY True for points inside polygonal region.
%
% IN = INPOLY( X , Y , XV , YV ) 
%
% Returns a matrix IN the size of X and Y, which is true
%  for points are inside the polygonal region whose
%  vertices are specified by the vectors XV and YV.
%
% IN(p,q) = 1   if the point ( X(p,q) , Y(p,q) ) is
%                  strictly inside the polygon.
%
% IN(p,q) = 0.5 if the point is on the polygon.
%
% IN(p,q) = 0   if the point is outside the polygon.
%
% Example:
%
%   xv = rand(6,1);    yv = rand(6,1);
%   xv = [xv ; xv(1)]; yv = [yv ; yv(1)];
%   xx = rand(1000,1); yy = rand(1000,1);
%   in = inpoly(xx,yy,xv,yv);
%   figure, plot(xv,yv,xx(in),yy(in),'.r',xx(~in),yy(~in),'.b')
%

%   Steven L. Eddins, July 1994
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 1.5 $  $Date: 1997/11/21 23:40:39 $
%
%   $Revision: 1.6 $  $Date: 2004/03/21 $  Ch.Begler
%     110% Time   in compare with INPOLYGON
%      50% Memory in compare with INPOLYGON


% PolygonVertices --> ColumnVectors

xv = xv(:);
yv = yv(:);
n  = size(xv,1);

% Close the Polygon

if ~( ( xv(1) == xv(n) ) & ( yv(1) == yv(n) ) )
  xv = [xv ; xv(1)];
  yv = [yv ; yv(1)];
  n  = n + 1;
end

% Make matrices to vectorize the computation.
% Translate the vertices so that the test points are
% at the origin.

si = size(x);

x = x(:).';
y = y(:).';

m = size(x,2);

xv = xv(:,ones(1,m)) - x(ones(n,1),:);
yv = yv(:,ones(1,m)) - y(ones(n,1),:);

%--------------------------------------------------------------
% Compute the Signum of the cross product and the dot product
% of adjacent vertices.

ii = ( 1 : n-1 );

% SignCrossProduct
scp = sign(  xv(ii  ,:) .* yv(ii+1,:) ...
           - xv(ii+1,:) .* yv(ii  ,:)     );

% DotProduct
dp = xv(ii,:) .* xv(ii+1,:) + yv(ii,:) .* yv(ii+1,:);

%--------------------------------------------------------------
% Compute the quadrant number for the vertices relative
% to the test points.

xv = ( xv > 0 );
yv = ( yv > 0 );

% Compute the vertex quadrant changes for each test point.

dq = diff( ( ~xv & yv ) + 2 * ( ~xv & ~yv ) + 3 * ( xv & ~yv ) );

%--------------------------------------------------------------
% Fix up the quadrant differences.  Replace 3 by -1 and -3 by 1.
% Any quadrant difference with an absolute value of 2 should have
% the same sign as the cross product.

% dq = dq - 4 * sign(dq) .* ( abs(dq) == 3 );
% dq = dq - 4 * ( dq == 3 ) + 4 * ( dq == -3 );

 dq(find(dq== 3)) = -1;
 dq(find(dq==-3)) =  1;

dq = dq + (  2 * scp - dq ) .* ( abs(dq) == 2 );

dq(find(isnan(dq))) = 0;

%--------------------------------------------------------------
% Find the points on the polygon.  If the cross product is 0 and
% the dot product is nonpositive anywhere, then the corresponding
% point must be on the contour.

in = ( sum(dq,1) ~= 0 );

% Indicate points on the polygon with a value of 0.5.

in = in  - 0.5 * any( ( scp == 0 ) & ( dp <= 0 ) );

% Reshape output matrix.

in = reshape( in , si );

