function p = arcang(x,n);

% ARCANG   Calculate Angles of an CircleSegment
%
%  ANG = ARCANG( XY , N );
%
%  XY defines an Angle with 3 Points: P1 --> P0 --> P2
%
%  XY = [ X1 Y1 ; X0 Y0 ; X2 Y2 ];
%
%  ANG returns an Vector of Angles for the Segment,
%   coursing anticlockwise from P1-->P2
%
%  N gives the Number of Elements of ANG
%
% example:
%
%    xy = randn(3,2);
%
%    p1 = arcang(        xy  , 10 );
%    p2 = arcang( flipud(xy) , 10 );
%
%    r = 0.5 * min(sqrt(sum(diff(xy,1,1).^2,2)));
%
%    c = xy(2,:);
%
%    figure
%    axes('nextplot','add');
%
%    plot( xy(:,1) , xy(:,2) , 'b-*' );
%
%    plot( c(1)+r*cos(p1) , c(2)+r*sin(p1) , 'g-' );
%    plot( c(1)+r*cos(p2) , c(2)+r*sin(p2) , 'r-' );
%     
%    plot( c(1)+r*cos(p1(end)) , c(2)+r*sin(p1(end)) , 'g.' );
%    plot( c(1)+r*cos(p2(end)) , c(2)+r*sin(p2(end)) , 'ro' );
%     

if nargin < 1
   p = [];
   return
end

if nargin < 2
   n = 2;
end

if isequal(size(x),[2 3])
   x = permute(x,[2 1]);
end

if ~isequal(size(x),[3 2])
   error('XY must define 3 Points of an Angle in [X Y].');
end

if ~( ( prod(size(n)) == 1 ) & ( mod(n,1) == 0 ) & ( n > 1 ) )
   error('N must be an Integer larger ONE.')
end


x = x([1 3],:) - x([2 2],:);

x = atan2(x(:,2),x(:,1));  % [ pi .. -pi )

x(2) = x(2) + 2*pi * ( x(2) < x(1) );

n = n-1;

if n == 1
   p = x;
   return
end

p = x(1) + ( x(2) - x(1) ) * ( 0 : n ) / n;

p(n+1) = x(2);
