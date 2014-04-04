function latout = cen2geod(latin,geoid,units)

%CEN2GEOD  Converts from geocentric latitude to geodetic latitude
%
%  lat = CEN2GEOD(lat0) converts from the geocentric latitude to the
%  geodetic latitude.  The geocentric latitude is the angle made by a
%  line from a point on the surface of the ellipsoid to the center of
%  the ellipsoid and the equatorial plane.  The geodetic latitude is the
%  angle that a line perpendicular to the surface of the ellipsoid at
%  the given point makes with the equatorial plane.
%
%  lat = CEN2GEOD(lat0,geoid) uses the ellipsoid definition given in
%  the input vector geoid.  Geoid can be determined from the ALMANAC
%  function.  If omitted, the default Earth geoid is assumed.
%
%  lat = CEN2GEOD(lat0,'units') uses the units defined by the input string
%  'units'.  If omitted, default units of degrees are assumed.
%
%  lat = CEN2GEOD(lat0,geoid,'units') uses the geoid and 'units'
%  definitions provided by the corresponding inputs.
%
%  See also:  GEOD2CEN, AUT2GEOD, CNF2GEOD, ISO2GEOD,
%             PAR2GEOD, REC2GEOD, ALMANAC

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.8 $    $Date: 1998/08/10 17:47:32 $

%   The formulae employed were taken from:  J. P. Snyder,
%   "Map Projections - A Working Manual,"  US Geological
%   Survey Professional Paper 1395, US Government Printing
%   Office, Washington, DC, 1987, pp. 13-18.


if nargin == 0
    error('Incorrect number of arguments')
elseif nargin == 1
	units = [];     geoid = [];
elseif nargin == 2
    if isstr(geoid);    units = geoid;    geoid = [];
         else;          units = 'degrees';
    end
end

%  Empty argument tests

if isempty(units);   units = 'degrees';       end
if isempty(geoid);   geoid = almanac('earth','geoid');  end

%  Test the geoid input

[geoid,msg] = geoidtst(geoid);
if ~isempty(msg);   error(msg);   end
e = geoid(2);

%  Perform the transformation

latin  = angledim(latin,units,'radians');     %  Convert to radians

halfpi = pi/2;
latin(find(latin==halfpi)) = halfpi-eps;
latin(find(latin==-halfpi)) = -halfpi+eps;

latout = atan (tan(latin) / (1-e^2) );        %  Conversion formula
latout = angledim(latout,'radians',units);    %  Convert to output units