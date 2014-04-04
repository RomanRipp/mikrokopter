function g = gravity(lat,z)

% GRAVITY Gravitational acceleration from Latitude and Z
%
%  G = GRAVITY(LAT) returns acceleration due to gravity as function of
%  latitude. The latitude must be given in decimal degrees north.
%
%  G = GRAVITY(LAT,Z) takes height z (m) above the sea level into account.
%
%  References:
%  Gill, A. E., Atmosphere-Ocean Dynamics, Volume 30 of International
%    Geophysics Series, Academic Press Inc., New York, 1982.
%  Fofonoff, N. P., and R. C. Millard Jr.,  Algorithms for computation of
%    fundamental properties of seawater. Unesco Tech. Pap. in Mar. Sci.,
%    No. 44, 54 pp.

%  Christian Mertens, IfM Kiel
%  $Revision: 1.0 $ $Date: 1997/07/28 20:55:18 $

if nargin == 1
  z = 0;
end

% mean radius of earth
a = 6371000;

lat = sin(pi/180*lat).^2;
g = 9.780318*(1 + (5.2788e-3 + 2.36e-5*lat).*lat)./((1 + z/a).^2);

