function latout = rec2geod(latin,geoid,units)

%REC2GEOD  Converts from rectifying latitude to geodetic latitude
%
%  lat = REC2GEOD(lat0) converts from the rectifying latitude to the
%  geodetic latitude.  The rectifying latitude is used to map an
%  ellipsoid to a sphere in such a way that correct distances along
%  meridians are preserved.  Rectifying latitudes are used in place
%  of the geodetic latitudes when projecting the ellipsoid using an
%  equal distant projection.  The geodetic latitude is the angle that
%  a line perpendicular to the surface of the ellipsoid at the given
%  point makes with the equatorial plane.
%
%  lat = REC2GEOD(lat0,geoid) uses the ellipsoid definition given in
%  the input vector geoid.  Geoid can be determined from an appropriate
%  data almanac function.  If omitted, the default Earth geoid is assumed.
%
%  lat = REC2GEOD(lat0,'units') uses the units defined by the input string
%  'units'.  If omitted, default units of degrees are assumed.
%
%  lat = REC2GEOD(lat0,geoid,'units') uses the geoid and 'units'
%  definitions provided by the corresponding inputs.
%
%  See also:  GEOD2REC, AUT2GEOD, CEN2GEOD, CNF2GEOD,
%             ISO2GEOD, PAR2GEOD, ALMANAC

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.8 $    $Date: 1998/08/10 17:48:03 $

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
fact1 = 3*n /2 - 27*n^3 /32;
fact2 = 21*n^2 /16 - 55*n^4 /32;
fact3 = 151*n^3 /96;
fact4 = 1097*n^4 /512;

latin = angledim(latin,units,'radians');   %  Convert to radians

latout = latin + ...                %  Series expansion for the
         fact1*sin(2*latin) + ...   %  transformation.  This is
		 fact2*sin(4*latin) + ...   %  an approximation of an
         fact3*sin(6*latin) + ...   %  infinite series
		 fact4*sin(8*latin);

latout = angledim(latout,'radians',units);  %  Convert to output units