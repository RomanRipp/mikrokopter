function [x,y,z,d] = gcline(p1,p2,n)

% GCLINE  Great Circle between two points on Sphere 
%
% [Lon,Lat,W] = GCLINE( [ Lon1 Lat1 ] ,  [ Lon2 Lat2 ] , N );
%
% Creates a N-Point Great Circle between spherical Positions 
%  Pos1 = [ Lon1 Lat1 ] and Pos2 = [ Lon2 Lat2 ]
% 
% Sperical Coordinates of Line: Lon = [ N by 1 ] ; Lat = [ N by 1 ]
%
% W   Angle Increment between Points of Line, [deg]
%
% Use 3-element Vectors for cartesian Coordinates:
% 
% [X,Y,Z,W] = GCLINE( [ X1 Y1 Z1 ] ,  [ X2 Y2 Z2 ] , N );
%
% Cartesian Coordinates of Line: [ X  Y  Z ] = [ N by 3 ];
%
%----------------------------------------------------------------------
%
% see also: ARCANG, CYCLIC, GCDEMO
%
%----------------------------------------------------------------------
% Example: GCDEMO
%
% >> type gcdemo
%



%**************************************************************
% Check Inputs

s = size(p1); p = prod(s);

if ~( isnumeric(p1) & isnumeric(p2) & ...
      isequal(s,size(p2)) & any( p == [ 2  3 ] ) )
    error('First two Inputs must be numeric Vectors of equal Size with [Lon Lat] or [X Y Z].');
end

if nargin < 3
   n = [];
end

if isempty(n)
   n = 2;
end

ok = ( isnumeric(n) & ( prod(size(n)) == 1 ) );
if ok
   ok = ( ( n >= 0 ) & ( mod(n,1) == 0 ) );
end

if ~ok
    error('N must be a positive Integer.');
end

if n == 0
   x = []; y = []; d = [];
   return
end

n = n + ( n == 1 );

%**************************************************************
% Check for geographical Lon/Lat

p180 = pi/180;

geo = ( p == 2 );

if geo

   x = cat(1,p1(1),p2(1));
   y = cat(1,p1(2),p2(2));

   % Cartesic Coordinates of Grid in Basis System

   xyz = ~cat( 2 , ( round((x+00)/180) == (x+00)/180 ) , ...
                   ( round((x+90)/180) == (x+90)/180 ) , ...
                   ( round((y+90)/180) == (y+90)/180 )        );

   xyz(:,1) = ( xyz(:,1) & xyz(:,3) );
   xyz(:,2) = ( xyz(:,2) & xyz(:,3) );
   xyz(:,3) = ~( round((y+00)/180) == (y+00)/180 );

   xyz = double(xyz);

   xyz(:,1) =  sin(x*p180) .* cos(y*p180) .* xyz(:,1);
   xyz(:,2) = -cos(x*p180) .* cos(y*p180) .* xyz(:,2);
   xyz(:,3) =  sin(y*p180) .* xyz(:,3);

   p1 = xyz(1,:);
   p2 = xyz(2,:);

end

%**************************************************************
% Transformation into Plain by P1,P2

p1 = p1(:); p10 = p1';
p2 = p2(:); p20 = p2';

l1 = sqrt(sum(p1.^2));
l2 = sqrt(sum(p2.^2));

n1 = p1 / l1;
n2 = p2 / l2;

ex = n1;
ez = crs(ex,n2);  % Cross-Product, see subfunction below

if all( ez == 0 )
   warning('P1 and P2 have same orientations.');
   x = NaN * zeros(n,1); y = x; 
   if geo
      z = NaN;
   else
      z = x;
   end
   d = NaN;
   return
end

ez = ez / sqrt(sum(ez.^2));

ey = crs(ez,ex);

T = [ ex  ey  ez ];

p1 = p1' * T;
p2 = p2' * T;

%**************************************************************
% Circle in Plain

v = ( 0 : (n-1) )';

w = atan2(p2(2),p2(1));
w = w + 2*pi*( w < 0 );

d = w / v(n);  % Intervall

p    = v * d;
p(n) = w;

l = ( 1 + cos( 2 * pi * v / (2*v(n)) ) ) / 2;
l = l2 + ( l1 - l2 ) * l;

x = l .* cos(p);
y = l .* sin(p);
z = 0*l;

%**************************************************************
% Reverse Transformation

xyz = [x y z] * T'; 

xyz(1,:) = p10;
cxyz(n,:) = p20;

%**************************************************************

d = d / p180;

if geo

   y = asin( xyz(:,3) );

   xyz(:,3) = cos(y);                 % cos(lat)

   xyz( find( abs(xyz(:,3)) < 1e8*eps ) , 3 ) = NaN;

   xyz(:,1) =  xyz(:,1) ./ xyz(:,3);  % sin(lon)
   xyz(:,2) = -xyz(:,2) ./ xyz(:,3);  % cos(lon)

   x = atan2( xyz(:,1) , xyz(:,2) );

   x = x / p180;
   y = y / p180;

   z = d;

else

   x = xyz(:,1);
   y = xyz(:,2);
   z = xyz(:,3);

end

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function c = crs(a,b)

% Vector cross product

c = [a(2)*b(3)-a(3)*b(2)
     a(3)*b(1)-a(1)*b(3)
     a(1)*b(2)-a(2)*b(1)];
