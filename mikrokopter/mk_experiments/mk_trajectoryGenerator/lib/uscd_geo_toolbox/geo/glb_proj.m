function [x,y] = glb_proj(x,y,f);

% GLB_PROJ  Elliptic Projection (Hammer-Aithoff), using WGS84
%
% The  Hammer-Aithoff-Projection is equal-area.  The world is
% displayed as an ellipse, with curved parallels and meridians.
% It is neither conformal nor equal area.  
% The only point free of distortion is the center point.
% Distortion of shape and area are moderate throughout.
%
% [ X , Y ] = GLB_PROJ( LON , LAT , Ratio )
%
% The Value of Ratio is valid for the transformed X / Y
%
%
%  default: Ratio = 1    (Hammer-Aithoff-Projection)
%
%------------------------------------------------------------------
%
% Use a NonZero imaginary Part of Ratio for inverse Projection:
%
% [ LON , LAT ] = GLB_PROJ( X , Y , Ratio+i )
%
%------------------------------------------------------------------
%
% see also: VDG_PROJ, CONEPROJ, SPH_PROJ, GEOLAT (required)
%

if nargin < 3
   f = [];
end

if isempty(f)
   f = 0;
end

rev = ~( imag(f) == 0 );

f   = 2 * real(f);

%******************************************************************
% WGS84

a = 6378137.;
b = 6356752.314;        %       A * sqrt(1-E*E);
e = sqrt(a*a - b*b)/a;  %   E = 0.081819191;

% Equal Surface Area Sphere
% f1 = (a/1e3)^2 / 2;
% f2 = (1 - e^2) / (2*e);
% f3 = log((1+e) / (1-e));
% r  = sqrt( f1 * ( 1 + f2 * f3 ) );

acc = 1e8 * eps;

%******************************************************************

p180 = pi/180;

b = 2;  % Basic Scale

if isequal(f,0)
   f = b;
end

f = b/f;

scl = 180 / ( b * sqrt(b) );

%********************************************
if rev
%********************************************

   x = x / scl / b;
   y = y / scl / f;

   %----------------------------------------
   % Check for Border-crossing X at 180 !!!

   xb = sqrt( 2 - y.^2 );  % Border at 180
                           % use: b == 2

   off = sign(x) .* ( abs(x) >  xb );

   %----------------------------------------

   d = b ./ sqrt( 2*b - ( x.^2 + y.^2 ) );

   w = 1 - ( y ./ d ).^2;
   w(find(w<=0)) = NaN;

   w = d .* sqrt(w);
   w(find(w==0)) = NaN;

   jj = find( isnan(x) | isnan(y) );

   x = min( max( x./w , -1 ) , 1 );
   y = min( max( y./d , -1 ) , 1 );

   if ~isempty(jj)
       x(jj) = NaN;
       y(jj) = NaN;
   end

   x = b * asin( x );
   y =     asin( y );

   %----------------------------------------
   % Check BorderCrossing Values

   if ~all( off(:) == 0 )

       jj = find( ~( off == 0 ) );

       x(jj) = 2*pi*off(jj) - x(jj);

   end
   %----------------------------------------

   x = x / p180;
   y = y / p180;

   y = geolat(y,'a',1);

%********************************************
else
%********************************************

   y = geolat(y,'a');

   x = p180 * x;
   y = p180 * y;

     ii  = find( abs( abs(y) - pi/2 ) < acc );
   y(ii) = sign(y(ii)) * ( pi/2 - acc );

   d = b ./ ( 1 + cos(y) .* cos(x/b) );
	
   d = sqrt(d);

   x = d .* cos(y) .* sin(x/b);
   y = d .* sin(y);

   x = b * x * scl;
   y = f * y * scl;

%********************************************
end
%********************************************
