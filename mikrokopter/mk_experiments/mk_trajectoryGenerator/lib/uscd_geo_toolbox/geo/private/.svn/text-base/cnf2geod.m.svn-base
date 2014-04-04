function latout = cnf2geod(latin,geoid,units)

%CNF2GEOD  Converts from conformal latitude to geodetic latitude
%
%  lat = CNF2GEOD(lat0) converts from the conformal latitude to the
%  geodetic latitude.  The conformal latitude is used to map an
%  ellipsoid conformally onto a sphere. Conformal latitudes are used in
%  place of the geodetic latitudes when projecting the ellipsoid using
%  a conformal projection.  The geodetic latitude is the angle which a
%  line perpendicular to the surface of the ellipsoid at the given point
%  makes with the equatorial plane.
%
%  lat = CNF2GEOD(lat0,geoid) uses the ellipsoid definition given in
%  the input vector geoid.  Geoid can be determined from the ALMANAC
%  function.  If omitted, the default Earth geoid is assumed.
%
%  lat = CNF2GEOD(lat0,'units') uses the units defined by the input string
%  'units'.  If omitted, default units of degrees are assumed.
%
%  lat = CNF2GEOD(lat0,geoid,'units') uses the geoid and 'units'
%  definitions provided by the corresponding inputs.
%
%  See also:  GEOD2CNF, AUT2GEOD, CEN2GEOD, ISO2GEOD,
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
if e > 0.5
   warning('Auxilary sphere approximation weakens with eccentricity > 0.5')
end

%  Compute the series expansion terms

fact1 = e^2 /2 + 5*e^4 /24 + e^6 /12 + 13*e^8 /360;
fact2 = 7*e^4 /48 + 29*e^6 /240 + 811*e^8 /11520;
fact3 = 7*e^6 /120 + 81*e^8 /1120;
fact4 = 4279*e^8 /161280;


latin = angledim(latin,units,'radians');   %  Convert to radians

latout = latin + ...                %  Series expansion for the
         fact1*sin(2*latin) + ...   %  transformation.  This is
		 fact2*sin(4*latin) + ...   %  an approximation of an
         fact3*sin(6*latin) + ...   %  infinite series
		 fact4*sin(8*latin);

latout = angledim(latout,'radians',units);  %  Convert to output units