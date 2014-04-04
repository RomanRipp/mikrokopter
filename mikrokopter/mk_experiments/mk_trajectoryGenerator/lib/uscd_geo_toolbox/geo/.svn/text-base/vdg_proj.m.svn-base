function [xx,yy] = vdg_proj(x,y,rev);

% VDG_PROJ  Van Der Grinten Polyconic Projection
%
%  In this projection, the world is enclosed in a circle.  Scale is
%  true along the Equator and increases rapidly away from the Equator.
%  Area distortion is extreme near the poles.  This projection is
%  neither conformal nor equal area.
%
%  This projection was presented by Alphons J. van der Grinten 1898.
%  He obtained a U.S. patent for it in 1904.  It is also known simply
%  as the Van der Grinten projection.
%
% [ X , Y ] = VDG_PROJ( LON , LAT )
%
% [ X , Y ] = VDG_PROJ( LON , LAT , XScale*i )
%
%  Multiplies LON with the Factor XScale before Projection
%
%------------------------------------------------------------------
%
% Use ONE as third Input for inverse Projection:
%
% [ LON , LAT ] = VDG_PROJ( X , Y , 1 )
% [ LON , LAT ] = VDG_PROJ( X , Y , 1+XScale*i )
%
%------------------------------------------------------------------
%
% see also: GLB_PROJ, CONEPROJ, SPH_PROJ
%

if nargin < 3
   rev = 0;
end

scl = 1;
if ~( imag(rev) == 0 )
    scl = imag(rev);
end

rev = isequal(real(rev),1);

acc = 1e4*eps;

p180 = pi/180;

sx = size(x); px = prod(sx);
sy = size(y); py = prod(sy);

if     ( px == 1 ) & ~( py == 1 )
        x = x(ones(sy));
elseif ( py == 1 ) & ~( px == 1 )
        y = y(ones(sx));
elseif ~isequal(sx,sy)
        if ( px == max(sx) ) & ( py == max(sy) )
           x = ones(py,1)*x(:)';
           y = y(:) * ones(1,px);
        else
           error('Matrix Dimensions must be agree.');
        end
end

xx = zeros(size(x));
yy = zeros(size(y));

if isempty(x)
   return
end

%********************************************
if rev
%********************************************

   x = x/180;
   y = y/180;

   ii = find( ~( ( x == 0 ) & ( y == 0 ) ) );

   if ~isempty(ii)

       f1 = x(ii).^2 + y(ii).^2;
       c1 = -abs(y(ii)) .* (1 + f1);
       c2 = c1 - 2* y(ii).^2 + x(ii).^2;
       c3 = -2*c1 + 1 + 2*y(ii).^2 + f1.^2;

       d = y(ii).^2 ./ c3 + (2*c2.^3./c3.^3 - 9*c1.*c2./c3.^2)/27;
       a = (c1 - c2.^2./(3*c3)) ./ c3;
       b = 2*sqrt(-a/3);

       t = acos( 3*d ./ (a.*b))/3;
       yy(ii) = 180 * sign(y(ii)) .* (-b.*cos(t+pi/3)-c2./(3*c3));

       jj = find( ~( x == 0 ) );

       if ~isempty(jj)
           c1 = x(jj).^2 + y(jj).^2;
           c2 = x(jj).^2 - y(jj).^2;
           c3 = (c1 - 1 + sqrt(1 + 2*c2 + c1.^2));
           xx(jj) = 180 * c3 ./ (2*x(jj));
       end

   end

   if ~( scl == 1 )
       x = x / scl;
   end

   return

%********************************************
end
%********************************************

if ~( scl == 1 )
    x = x * scl;
end

x = x * p180;
y = y * p180;

ii = find(abs(abs(y)-pi/2)<=acc);
if ~isempty(ii)
    y(ii) = ( pi/2 - acc ) * sign(y(ii));
end

t = asin(2*abs(y)/pi);

ii = find(abs(x)<=acc);
if ~isempty(ii)
    xx(ii) = 0;
    yy(ii) = tan(t(ii)/2);
end

ii = find(abs(y)<=acc);
if ~isempty(ii)
    xx(ii) = abs(x(ii))/pi;
    yy(ii) = 0;
end

ii = find( ( abs(x) > acc) & ( abs(y) > acc ) );

if ~isempty(ii)

    a = sign(x(ii)) .* ( pi./x(ii) - x(ii)/pi) / 2;
    b = cos(t(ii)) ./ (sin(t(ii)) + cos(t(ii)) - 1);
    c = b .* (2./sin(t(ii)) - 1);
    d = a.^2 + b;

    f1 = a .* (b - c.^2);
    f2 = c.^2 + a.^2;
    f3 = f1 + sqrt( f1.^2  - f2.*(b.^2 - c.^2));
    f4 = c.*d - a.*sqrt((a.^2 + 1).*f2 - d.^2);

    xx(ii) = f3 ./ f2;
    yy(ii) = f4 ./ f2;

end

xx = 180 * xx .* sign(x);
yy = 180 * yy .* sign(y);


