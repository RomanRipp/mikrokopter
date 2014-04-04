function latout = geod2par(latin,geoid,units)

%GEOD2PAR  Converts from geodetic latitude to parametric latitude
%
%  lat = GEOD2PAR(lat0) converts from the geodetic latitude to the
%  parametric latitude.  The parametric latitude of a point on the
%  ellipsoid is the latitude on a sphere (of radius r = semimajor axis
%  of ellipsoid) for which the parallel has the same radius as the
%  parallel of geodetic latitude.  The geodetic latitude is the angle
%  which a line perpendicular to the surface of the ellipsoid at the
%  given point makes with the equatorial plane.
%
%  lat = GEOD2PAR(lat0,geoid) uses the ellipsoid definition given in
%  the input vector geoid.  Geoid can be determined from the ALMANAC
%  function.  If omitted, the default Earth geoid is assumed.
%
%  lat = GEOD2PAR(lat0,'units') uses the units defined by the input string
%  'units'.  If omitted, default units of degrees are assumed.
%
%  lat = GEOD2PAR(lat0,geoid,'units') uses the geoid and 'units'
%  definitions provided by the corresponding inputs.
%
%  See also:  PAR2GEOD, GEOD2AUT, GEOD2CEN, GEOD2CNF,
%             GEOD2ISO, GEOD2REC, ALMANAC

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

%  Perform the transformation

latin  = angledim(latin,units,'radians');     %  Convert to radians

halfpi = pi/2;
latin(find(latin==halfpi)) = halfpi-eps;
latin(find(latin==-halfpi)) = -halfpi+eps;

latout = atan (sqrt(1-e^2) * tan(latin));     %  Conversion formula
latout = angledim(latout,'radians',units);    %  Convert to output units