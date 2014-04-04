function nm = deg2nm(deg,radius)

%DEG2NM Converts distances from degrees to nautical miles
%
%  nm = DEG2NM(deg) converts distances from degrees to nautical miles.
%  A degree of distance is measured along a great circle of a sphere.
%
%  nm = DEG2NM(deg,radius) uses the second input to determine the
%  radius of the sphere.  If radius is a string, then it is evaluated
%  as an ALMANAC body to determine the spherical radius.  If numerical,
%  it is the radius of the desired sphere in nautical miles.  If omitted,
%  the default radius of the Earth is used.
%
%  See also NM2DEG, DEG2RAD, DEG2KM, DEG2SM, DISTDIM

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Brown, E. Byrns
%   $Revision: 1.8 $    $Date: 1998/08/10 17:47:34 $

if nargin==0
	error('Incorrect number of arguments')
elseif nargin == 1
    radius = almanac('earth','radius','km');
elseif nargin == 2 & isstr(radius)
	radius = almanac(radius,'radius','km');
elseif nargin == 2 & ~isstr(radius)
    radius = nm2km(radius);
end


nm = km2nm(deg2km(deg,radius) );