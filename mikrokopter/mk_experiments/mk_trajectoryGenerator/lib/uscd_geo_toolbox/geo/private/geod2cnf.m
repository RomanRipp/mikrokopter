function latout = geod2cnf(latin,geoid,units)

%GEOD2CNF  Converts from geodetic latitude to conformal latitude
%
%  lat = GEOD2CNF(lat0) converts from the geodetic latitude to the
%  conformal latitude.  The conformal latitude is used to map an
%  ellipsoid conformally onto a sphere. Conformal latitudes are used
%  in place of the geodetic latitudes when projecting the ellipsoid
%  using a conformal projection.  The geodetic latitude is the angle
%  which a line perpendicular to the surface of the ellipsoid at the
%  given point makes with the equatorial plane.
%
%  lat = GEOD2CNF(lat0,geoid) uses the ellipsoid definition given in
%  the input vector geoid.  Geoid can be determined from the ALMANAC
%  function.  If omitted, the default Earth geoid is assumed.
%
%  lat = GEOD2CNF(lat0,'units') uses the units defined by the input string
%  'units'.  If omitted, default units of degrees are assumed.
%
%  lat = GEOD2CNF(lat0,geoid,'units') uses the geoid and 'units'
%  definitions provided by the corresponding inputs.
%
%  See also:  CNF2GEOD, GEOD2AUT, GEOD2CEN, GEOD2ISO,
%             GEOD2PAR, GEOD2REC, ALMANAC

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.8 $    $Date: 1998/08/10 17:47:44 $

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

%  Compute the transformation

fact1  = 1 - e*sin(latin);     fact2 = 1 + e*sin(latin);
fact3  = 1 - sin(latin);       fact4 = 1 + sin(latin);

%  Eliminate any divide by zeros.  Points exactly at +90

if ~isempty(fact3)            %  Input matrix may be empty
    indx = find(fact3 == 0);  %  when called from places like reckon
    if ~isempty(indx);   fact3(indx) = eps;   end
else
    indx = [];
end

%  Conformal latitudes

latout = 2 * atan(sqrt((fact4./fact3) .* ((fact1./fact2).^e)) ) - pi/2;

%  Correct for points at +90;

if ~isempty(indx);    latout(indx) = pi/2;    end

indx = find(isnan(latout));               % operations on NaNs produce
if ~isempty(indx); latout(indx)=NaN; end  % unwanted complex results

latout = angledim(latout,'radians',units);  %  Convert to output units