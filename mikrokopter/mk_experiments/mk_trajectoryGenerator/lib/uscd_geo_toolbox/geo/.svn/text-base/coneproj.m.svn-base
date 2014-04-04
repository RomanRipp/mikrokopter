function [x,y,c,w] = coneproj(x,y,c,a,rev);

% CONEPROJ  Conic Projection of geodesic Coordinates
%
% [ X , Y ] = CONEPROJ( LON , LAT , Pole , Parallels )
%
%  Pole = [ PoleLon  <PoleLat>  <Rotation> ]
%
%     Rotation: Angle of the Y-Axis clockwise to North
%
%     defaults: PoleLon  = mean(Lon)
%               PoleLat  = 90
%               Rotation = 0
%
%  Parallels = Lat | [ Lat1  Lat2 ]  Standard Parallels
%
%     default: Parallels = mean(Lat)
%
%---------------------------------------------------------
%
% [ X , Y , Pol , Par ] = CONEPROJ( ... )
% 
%   returns the used Pole and Parallels,
%    to use for inverse Projection:
%
% [ LON , LAT ] = CONEPROJ( X , Y , Pol , Par , 1 )
%
%
%---------------------------------------------------------
%
% see also: SPH_PROJ (required), GLB_PROJ, VDG_PROJ
%


Nin = nargin;

if Nin < 3
   c = [];
end

if Nin < 4
   a = [];
end

if Nin < 5
   rev = [];
end

%********************************************
% Check Inputs

p180 = pi/180;

%--------------------------------------------
% Reverse

if isempty(rev)
   rev = 0;
else
   rev = ~isequal(rev,0);
   if rev & ( isempty(c) | isempty(a) )
      error('Missing Pole and Parallels for inverse Projection.')
   end
end

%--------------------------------------------
% Pol: [ Center Rotation ]

p = [ 0  90  0 ];   % Basic Pole

if isempty(c)
   c = meannan(x(:));
end

c = c(:)';

nc = size(c,2);
np = size(p,2);

c = cat( 2 , c(1:min(nc,np)) , p((nc+1):np) );

%--------------------------------------------
% Parallels

if isempty(a)
   a = meannan(y(:));
end

%********************************************
% CylinderParameter:
%
%  a = Half Angle
%  c = Height above Centre of Sphere: [ 0  0  0 ]
%

if ~( prod(size(a)) == 1 )

   a = ( a(1) + a(2) ) / 2;

%   t1 = tan( p180 * a(1) / 2 );
%   t2 = tan( p180 * a(2) / 2 );
%   b  = ( 1 + t1*t2 ) / ( t1 + t2 );

end

w = a;

a = a * p180;

b = cos(a);
b = b + ( abs(b) < 1e3*eps );


%********************************************
if rev
%********************************************

   %--------------------------------------------
   % Reverse Projection

   x = atan2(x,-y);

   y = -y ./ cos(x);

   x = x / sin(a);

   y = atan( ( b./y - cos(a) ) / sin(a) );

   x = x / p180;
   y = y / p180;

   x = x - 360 * floor( (x+180) / 360 );

   %--------------------------------------------
   % Transform Coordinates with Pole

   if all( c([2 3]) == p([2 3]) )

      x = ( x + c(1) );

   else

      [x,y,z] = sph_proj( p , x , y );

      [x,y]   = sph_proj( c , x , y , z );

   end

%********************************************
else
%********************************************

   %--------------------------------------------
   % Transform Coordinates with Pole

   if all( c([2 3]) == p([2 3]) )

      x = ( x - c(1) );

   else

      [x,y,z] = sph_proj( c , x , y );
%figure,plot3(x,y,z)
      [x,y]   = sph_proj( p , x , y , z );
%figure,plot(x,y)
   end

   %--------------------------------------------
   % Cylindric Projection

   x = x - 360 * floor( (x+180) / 360 );

   x = p180 * x;
   y = p180 * y;

   % Distance from Top of Cylinder

   d = b ./ ( cos(a) + sin(a) * tan(y) );

   % Flat Projection

   x = x * sin(a);

   y = -d .* cos(x);
   x =  d .* sin(x);

%********************************************
end
%********************************************

