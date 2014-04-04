function latout = aut2geod(latin,geoid,units)

%AUT2GEOD  Converts from authalic latitude to geodetic latitude
%
%  lat = AUT2GEOD(lat0) converts from the authlic latitude to the
%  geodetic latitude.  The authalic latitude is used to map an
%  ellipsoid to a sphere in such a way that the sphere has equal surface
%  area as the ellipsoid.  Authalic latitudes are used in place of the
%  geodetic latitudes when projecting the ellipsoid using an equal
%  area projection.  The geodetic latitude is the angle that a line
%  perpendicular to the surface of the ellipsoid at the given point makes
%  with the equatorial plane.
%
%  lat = AUT2GEOD(lat0,geoid) uses the ellipsoid definition given in
%  the input vector geoid.  Geoid can be determined from the ALMANAC
%  function.  If omitted, the default Earth geoid is assumed.
%
%  lat = AUT2GEOD(lat0,'units') uses the units defined by the input string
%  'units'.  If omitted, default units of degrees are assumed.
%
%  lat = AUT2GEOD(lat0,geoid,'units') uses the geoid and 'units'
%  definitions provided by the corresponding inputs.
%
%  See also:  GEOD2AUT, CEN2GEOD, CNF2GEOD, ISO2GEOD,
%             PAR2GEOD, REC2GEOD, ALMANAC

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.8 $    $Date: 1998/08/10 17:47:30 $

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

fact1 = e^2 /3 + 31*e^4 /180 + 517*e^6 /5040;
fact2 = 23*e^4 /360 + 251*e^6 /3780;
fact3 = 761*e^6 /45360;


latin = angledim(latin,units,'radians');  %  Convert to radians

latout = latin + ...                %  Series expansion for the
         fact1*sin(2*latin) + ...   %  transformation.  This is
		 fact2*sin(4*latin) + ...   %  an approximation of an
         fact3*sin(6*latin);        %  infinite series

latout = angledim(latout,'radians',units);  %  Convert to output units