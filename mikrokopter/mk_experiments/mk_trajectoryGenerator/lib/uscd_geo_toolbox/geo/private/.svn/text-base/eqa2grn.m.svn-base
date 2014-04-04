function [lat,lon] = eqa2grn(x,y,origin,geoid,units)

%EQA2GRN  Equal area coordinates to Greenwich coordinates
%
%  [lat,lon] = EQA2GRN(x,y) transforms from equal area cylindrical
%  coordinates to Greenwich latitude and longitude coordinates.
%  The transformation assumes that the origin of the cylindrical
%  coordinate plane is at (lat,lon) = (0,0).  The outputs
%  are returned in degrees.
%
%  [lat,lon] = EQA2GRN(x,y,origin) assumes that the origin of the
%  cylindrical coordinate plate is given by the input origin.
%  This input is of the form [lat0,lon0] where the 2 elements are
%  in degrees.
%
%  [lat,lon] = EQA2GRN(x,y,origin,geoid) assumes that the data
%  is distributed on the ellipsoid defined by the input geoid.
%  The geoid vector is of the form [semimajor axes, eccentricity].
%  If omitted, the unit sphere, geoid = [1 0], is assumed.
%
%  [lat,lon] = EQA2GRN(x,y,origin,'units') and
%  [lat,lon] = EQA2GRN(x,y,origin,geoid,'units') assumes that the
%  origin input and output vectors are in the angle units specified
%  by 'units'.  If omitted, 'degrees' are assumed.
%
%  mat = EQA2GRN(...) returns a single output matrix,
%  where mat = [lat lon].
%
%  See also GRN2EQA

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.7 $    $Date: 1998/08/10 17:47:41 $

if nargin < 2
    error('Incorrect number of arguments')
elseif nargin == 2
    origin = [];   units = [];    geoid = [];
elseif nargin == 3
    units = [];    geoid = [];
elseif nargin == 4 & isstr(geoid)
    units = geoid;  geoid = [];
elseif nargin == 4 & ~isstr(geoid)
    units = [];
end

%  Empty argument tests

if isempty(units);   units  = 'degrees';   end
if isempty(geoid);   geoid  = [1 0];       end
if isempty(origin);  origin = [0 0 0];     end

%  Dimension tests

if ~isequal(size(x),size(y))
    error('Inconsistent dimensions on latitude and longitude')
end

%  Test the geoid input

[geoid,msg] = geoidtst(geoid);
if ~isempty(msg);   error(msg);   end

%  Transform the location to an equal area coordinate frame

[lat0,lon0] = eqacalc(x,y,origin,'inverse',units,geoid);

%  Set the output arguments

if nargout <= 1;           lat = [lat0 lon0];
    else;                  lat = lat0;  lon = lon0;
end