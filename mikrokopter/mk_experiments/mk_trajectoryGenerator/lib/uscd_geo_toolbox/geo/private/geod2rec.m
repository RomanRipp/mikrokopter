function latout = geod2rec(latin,geoid,units)

%GEOD2REC  Converts from rectifying latitude to geodetic latitude
%
%  lat = GEOD2REC(lat0) converts from the geodetic latitude to the
%  rectifying latitude.  The rectifying latitude is used to map an
%  ellipsoid to a sphere in such a way that correct distances along
%  meridians are preserved.  Rectifying latitudes are used in place
%  of the geodetic latitudes when projecting the ellispoid using an
%  equal distant projection.  The geodetic latitude is the angle which
%  a line perpendicular to the surface of the ellipsoid at the given
%  point makes with the equatorial plane.
%
%  lat = GEOD2REC(lat0,geoid) uses the ellipsoid definition given in
%  the input vector geoid.  Geoid can be determined from the ALMANAC
%  function.  If omitted, the default Earth geoid is assumed.
%
%  lat = GEOD2REC(lat0,'units') uses the units defined by the input string
%  'units'.  If omitted, default units of degrees are assumed.
%
%  lat = GEOD2REC(lat0,geoid,'units') uses the geoid and 'units'
%  definitions provided by the corresponding inputs.
%
%  To properly project rectified latitudes, the radius must also
%  be scaled to assure the equal meridional distance property.
%  This is accomplished by RSPHERE.
%
%  See also:  REC2GEOD, GEOD2AUT, GEOD2CEN, GEOD2CNF,
%             GEOD2ISO, GEOD2PAR, ALMANAC

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.8 $    $Date: 1998/08/10 17:47:45 $

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

n     = ecc2n(e);
fact1 = 3*n /2 - 9*n^3 /16;
fact2 = 15*n^2 /16 - 15*n^4 /32;
fact3 = 35*n^3 /48;
fact4 = 315*n^4 /512;

latin = angledim(latin,units,'radians');   %  Convert to radians

latout = latin - ...                %  Series expansion for the
         fact1*sin(2*latin) + ...   %  transformation.  This is
		 fact2*sin(4*latin) - ...   %  an approximation of an
         fact3*sin(6*latin) + ...   %  infinite series
		 fact4*sin(8*latin);

latout = angledim(latout,'radians',units);  %  Convert to output units