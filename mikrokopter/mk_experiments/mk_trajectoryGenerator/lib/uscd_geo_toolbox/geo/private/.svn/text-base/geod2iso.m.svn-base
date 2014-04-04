function latout = geod2iso(latin,geoid,units)

%GEOD2ISO  Comutes the isometric latitude given the geodetic latitude
%
%  lat = GEOD2ISO(lat0) computes the isometric latitude given the
%  geodetic latitude.  The isometric latitude is a nonlinear function
%  of the geodetic latitude.  It is directly proportional to the
%  spacing of parallels of geodetic latitude from the Equator on the
%  ellipsoidal Mercator projection.  The geodetic latitude is the angle
%  which a line perpendicular to the surface of the ellipsoid at the
%  given point makes with the equatorial plane.
%
%  lat = GEOD2ISO(lat0,geoid) uses the ellipsoid definition given in
%  the input vector geoid.  Geoid can be determined from the ALMANAC
%  function.  If omitted, the default Earth geoid is assumed.
%
%  lat = GEOD2ISO(lat0,'units') uses the units defined by the input string
%  'units'.  If omitted, default units of degrees are assumed.
%
%  lat = GEOD2ISO(lat0,geoid,'units') uses the geoid and 'units'
%  definitions provided by the corresponding inputs.
%
%  See also:  GEOD2ISO, AUT2GEOD, CEN2GEOD, CNF2GEOD,
%             PAR2GEOD, REC2GEOD, ALMANAC

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


latin = angledim(latin,units,'radians');   %  Convert to radians

%  Eliminate any log of zeros.  Points exactly at -90.  Replace
%  with +90 and then negate these points

indx = find(abs(latin + pi/2) <= eps );
if ~isempty(indx);   latin(indx) = pi/2;   end


latcnf = geod2cnf(latin,geoid,'radians');   %  Compute the conformal lat
latout = log ( tan(pi/4 + latcnf/2) );      %  Transform to isometric

%  Correct for points at -90;

if ~isempty(indx);    latout(indx) = -latout(indx);    end

indx = find(isnan(latout));               % operations on NaNs produce
if ~isempty(indx); latout(indx)=NaN; end  % unwanted complex results

latout = angledim(latout,'radians',units);  %  Convert to output units