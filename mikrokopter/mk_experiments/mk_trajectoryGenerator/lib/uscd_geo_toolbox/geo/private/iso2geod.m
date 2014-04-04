function latout = iso2geod(latin,geoid,units)

%ISO2GEOD  Computes the geodetic latitude given the isometric latitude
%
%  lat = ISO2GEOD(lat0) computes the geodetic latitude given the
%  isometric latitude.  The isometric latitude is a nonlinear function of
%  the geodetic latitude.  It is directly proportional to the spacing of
%  parallels of geodetic latitude from the Equator on the ellipsoidal Mercator
%  projection.  The geodetic latitude is the angle that a line perpendicular
%  to the surface of the ellipsoid at the given point makes with the equatorial
%  plane.
%
%  lat = ISO2GEOD(lat0,geoid) uses the ellipsoid definition given in
%  the input vector geoid.  Geoid can be determined from the ALMANAC
%  function.  If omitted, the default Earth geoid is assumed.
%
%  lat = ISO2GEOD(lat0,'units') uses the units defined by the input string
%  'units'.  If omitted, default units of degrees are assumed.
%
%  lat = ISO2GEOD(lat0,geoid,'units') uses the geoid and 'units'
%  definitions provided by the corresponding inputs.
%
%  See also:  GEOD2ISO, AUT2GEOD, CEN2GEOD, CNF2GEOD,
%             PAR2GEOD, REC2GEOD, ALMANAC

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.8 $    $Date: 1998/08/10 17:47:50 $

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


latin = angledim(latin,units,'radians');   %  Convert to radians

latcnf = 2 * atan(exp(latin)) - pi/2;       %  Compute the conformal lat
latout = cnf2geod(latcnf,geoid,'radians');  %  Transform to geodetic

latout = angledim(latout,'radians',units);  %  Convert to output units