function [latout,lonout] = reckon(varargin)

%RECKON  Computes points at a specified azimuth and range
%
%  [lat,lon] = RECKON(lat0,lon0,rng,az) computes the latitude and
%  longitude positions for selected ranges and azimuths from a
%  starting point along a great circle path on a globe.  The range
%  is input as degrees of arc length on a sphere.  The input
%  azimuth is measured clockwise from due north.
%
%  [lat,lon] = RECKON(lat0,lon0,rng,az,geoid) computes the new point
%  assuming that the great circle points lie on the ellipsoid defined by
%  the input geoid.  The geoid vector is of the form
%  [semimajor axes, eccentricity].
%
%  [lat,lon] = RECKON(lat0,lon0,rng,az,'units') uses the input string 'units'
%  to define the angle units of the input and output data.  In this
%  form, the input range is measured as an arc length in the units
%  specified by 'units'.  If 'units' is omitted, 'degrees' are assumed.
%
%  [lat,lon] = RECKON(lat0,lon0,rng,az,geoid,'units') is a valid calling
%  form.  In this case, the input range is in the same units as
%  the semimajor axes of the geoid vector.
%
%  [lat,lon] = RECKON('track',...) uses the input string 'track' to define
%  either a great circle or rhumb line reckon calculation.  If
%  'track' = 'gc', then the new point is computed to along a great
%  circle track.  If 'track' = 'rh', then the new point is computed
%  along a rhumb line track.  If omitted, 'gc' is assumed.
%
%  See also AZIMUTH, DISTANCE

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.9 $    $Date: 1998/08/10 17:48:03 $

if nargin < 1
     error('Incorrect number of arguments')
else
    if isstr(varargin{1})
	     str = varargin{1};    varargin(1) = [];
	else
	     str = [];
	end
end


%  Test the track string and call the appropriate function

if isempty(str)
    [newlat,newlon,msg] = reckongc(varargin{:});
else
    validstr = strvcat('gc','rh');
	indx     = strmatch(lower(str),validstr);
	if length(indx) ~= 1
	      error('Unrecognized track string')
    elseif indx == 1
          [newlat,newlon,msg] = reckongc(varargin{:});
    elseif indx == 2
          [newlat,newlon,msg] = reckonrh(varargin{:});
    end
end

%  Error out if necessary

if ~isempty(msg);   error(msg);   end

%  Set the output arguments

if nargout <= 1
     latout = [newlat newlon];
elseif nargout == 2
     latout = newlat;  lonout = newlon;
end


%***********************************************************************
%***********************************************************************
%***********************************************************************


function [latout,lonout,msg] = reckongc(lat,lon,rng,az,in5,in6)

%RECKONGC:  Computes points along a specified azimuth and range on a sphere
%
%  Purpose
%
%  Computes Greenwich frame latitude/longitude positions
%  for selected ranges and azimuths along great circle paths
%  from defined starting latitude/longitude positions.
%  The computations are performed assuming a spherical geoid.
%  The default input angle units is degrees.  The default
%  range input is in the angular distance in the same units
%  as the input angles.
%
%  Synopsis
%
%       [newlat,newlon] = reckongc(lat,lon,rng,az)
%       [newlat,newlon] = reckongc(lat,lon,rng,az,geoid)
%       [newlat,newlon] = reckongc(lat,lon,rng,az,'units')
%       [newlat,newlon] = reckongc(lat,lon,rng,az,geoid,'units')
%
%       [newlat,newlon,errmsg] = reckongc(....
%            If three output arguments are supplied, then error condition
%            messages are returned to the calling function for processing.

%   The formulae employed were taken from:  J. P. Snyder,
%   "Map Projections - A Working Manual,"  US Geological
%   Survey Professional Paper 1395, US Government Printing
%   Office, Washington, DC, 1987, pp. 29-32.

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown




%  Initialize outputs

if nargout ~= 0;  latout = [];   lonout = [];   msg = [];  end

%  Test inputs

if nargin == 4
	units  = [];         geoid = [];

elseif nargin == 5

	if isstr(in5)
	    units = in5;        geoid = [];

    else
	    units = [];         geoid = in5;
    end


elseif nargin == 6

	geoid = in5;        units  = in6;

else
    msg = 'Incorrect number of arguments';
	if nargout < 3;  error(msg);  end
	return
end

%  Empty argument tests

if isempty(geoid);   geoid = [0 0];       end
if isempty(units);   units = 'degrees';   end

%  Dimension tests

if ~isequal(size(lat),size(lon),size(az),size(rng))
      msg = 'Lat, long, azimuth and range inputs must have same dimension';
	  if nargout < 3;  error(msg);  end
	  return
end


%  Test the geoid parameter

[geoid,msg] = geoidtst(geoid);
if ~isempty(msg)
	  if nargout < 3;  error(msg);  end
	  return
end

%  Angle unit conversion

lat = angledim(lat,units,'radians');
lon = angledim(lon,units,'radians');
az  = angledim(az,units,'radians');

%  Compute the range in radians

 if ~isreal(rng)
%     warning('Imaginary parts of complex RANGE argument ignored')
     rng = real(rng);
 end

%  If this routine ever becomes applicable for ellipsoids, be
%  sure to make the appropriate normalizing changes in trackui and scirclui.

if geoid(1) == 0;    rng = angledim(rng,units,'radians');
    else;            rng = rng / geoid(1);
end

%  Ensure correct azimuths at either pole.

epsilon = 10*epsm('radians');     % Set tolerance
indx = find(lat >= pi/2-epsilon);   az(indx) = pi;    % starting at north pole
indx = find(lat <= epsilon-pi/2);   az(indx) = 0;     % starting at south pole

%  Perform calculation only for spherical geoids.
%  So far, we have been unable to find the closed form
%  solution for reckoning on an elliptical geoid.  We tend
%  to think that this is an extremely hard problem in geodesy.  EVB 3/4/96

if geoid(2) ~= 0
	geoid(2) = 0;    warning('Only spherical reckoning allowed')
end

%  Compute the new points

if geoid(2) == 0                   %  Spherical geoid

    temp1  = sin(lat).*cos(rng);          %  Compute latitude
    temp2  = cos(lat).*sin(rng).*cos(az);
    newlat = asin(temp1+temp2);

    temp1  = sin(rng).*sin(az);            %  Compute longitude
    temp2  = cos(lat).*cos(rng);
    temp3  = sin(lat).*sin(rng).*cos(az);
    newlon = lon + atan2(temp1,temp2-temp3);

end


%  Convert to desired units

newlat = angledim(newlat,'radians',units);
newlon = npi2pi(newlon,'radians','exact');
newlon = angledim(newlon,'radians',units);


%  Set the output arguments

if nargout <= 2
     latout = [newlat newlon];
else
     latout = newlat;  lonout = newlon;
end


%***********************************************************************
%***********************************************************************
%***********************************************************************


function [latout,lonout,msg] = reckonrh(lat,lon,rng,az,in5,in6)

%RECKONRH:  Computes points at a rhumb line azimuth and range
%
%  Purpose
%
%  Computes Greenwich frame latitude/longitude positions
%  for selected ranges and azimuths along rhumb line paths
%  from defined starting latitude/longitude positions.
%  The default range input is the angular
%  distance in the same units as the input angles.
%  The default geoid is a sphere, but this can be
%  redefined to an ellipsoid using the geoid input.
%
%  Synopsis
%
%       [newlat,newlon] = reckonrh(lat,lon,rng,az)
%       [newlat,newlon] = reckonrh(lat,lon,rng,az,geoid)
%       [newlat,newlon] = reckonrh(lat,lon,rng,az,'units')
%       [newlat,newlon] = reckonrh(lat,lon,rng,az,geoid,'units')
%
%       [newlat,newlon,errmsg] = reckonrh(....
%            If three output arguments are supplied, then error condition
%            messages are returned to the calling function for processing.


%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Brown, E. Byrns



%  Initialize outputs

if nargout ~= 0;  latout = [];   lonout = [];   msg = [];  end

%  Test inputs

if nargin == 4
	units  = [];       geoid = [];

elseif nargin == 5

	if isstr(in5)
	    units = in5;       geoid = [];
    else
	    units = [];        geoid = in5;
    end

elseif nargin == 6
	geoid = in5;      units  = in6;

else
	error('Incorrect number of arguments.')
end


%  Empty argument tests

if isempty(geoid);   geoid = [0 0];       end
if isempty(units);   units = 'degrees';   end

%  Dimension tests

if ~isequal(size(lat),size(lon),size(az),size(rng))
      msg = 'Lat, long, azimuth and range inputs must have same dimension';
	  if nargout < 3;  error(msg);  end
	  return
end

%  Test the geoid parameter

[geoid,msg] = geoidtst(geoid);
if ~isempty(msg)
	  if nargout < 3;  error(msg);  end
	  return
end

%  Angle unit conversion

lat = angledim(lat,units,'radians');
lon = angledim(lon,units,'radians');
az  = angledim(az,units,'radians');

%  Compute the rectifying sphere latitudes and radius

rec    = geod2rec(lat,geoid,'radians');
radius = rsphere('rectifying',geoid);

%  Compute the range in radians

% if ~isreal(rng)
%     warning('Imaginary parts of complex RANGE argument ignored')
%     rng = real(rng);
% end

if geoid(1) == 0;    rng = angledim(rng,units,'radians');
    else;            rng = rng / radius;
end

newlat  = zeros(size(lat));     % Preallocate memory for output
newlon  = zeros(size(lon));
epsilon = 10*epsm('radians');     % Set tolerance

%  Ensure correct azimuths at either pole.

indx = find(lat >= pi/2-epsilon);   az(indx) = pi;    % starting at north pole
indx = find(lat <= epsilon-pi/2);   az(indx) = 0;     % starting at south pole

%  Compute the new points

cosaz  = cos(az);
newlat = rec + rng.*cosaz;                    %  Compute rectifying latitude
newlat = rec2geod(newlat,geoid,'radians');    %  Transform to geodetic latitude

%  Latitudes cannot lie at or beyond a pole.  A latitude at the
%  pole causes the general longitude equation (below) to become
%  indeterminate.  Easiest to back away from the pole by a little bit.

indx=find(newlat >=  pi/2);     newlat(indx) = pi/2-epsilon;
indx=find(newlat <= -pi/2);     newlat(indx) = epsilon-pi/2;

indx1=find(abs(cosaz)<=epsilon); % find degenerate cases
								 % #1: cos(course)=0
indx=1:prod(size(az));  %  identify non-degenerate cases
indx([indx1])=[];

% handle the degenerate case (Go east or west only)

if geoid(1) ~= 0; rng(indx1) = rng(indx1) * radius / geoid(1);  end
newlon(indx1) = lon(indx1) + ...
                sign(sin(az(indx1))) .* rng(indx1) .* ...
				sqrt( 1 - (geoid(2)*sin(lat(indx1))).^2) ./ ...
                cos(lat(indx1));

% handle the general case:

newiso = geod2iso(newlat(indx),geoid,'radians');
iso    = geod2iso(lat(indx),geoid,'radians');

newlon(indx) = lon(indx) + tan(az(indx)) .* (newiso - iso);

%  Convert to the geodetic latitudes and then to desired units

newlat = angledim(newlat,'radians',units);
newlon = npi2pi(newlon,'radians','exact');
newlon = angledim(newlon,'radians',units);

%  Set the output arguments

if nargout <= 2
     latout = [newlat newlon];
else
     latout = newlat;  lonout = newlon;
end
