function rad=sm2rad(sm,radius)

%SM2RAD Converts distances from statute miles to radians
%
%  rad = SM2RAD(sm) converts distances from statute miles to radians.
%  A radian of distance is measured along a great circle of a sphere.
%
%  rad = SM2RAD(sm,radius) uses the second input to determine the
%  radius of the sphere.  If radius is a string, then it is evaluated
%  as an ALMANAC body to determine the spherical radius.  If numerical,
%  it is the radius of the desired sphere in statute miles.  If omitted,
%  the default radius of the Earth is used.
%
%  See also RAD2SM, SM2DEG, SM2NM, SM2SM, DISTDIM

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Brown, E. Byrns
%   $Revision: 1.8 $    $Date: 1998/08/10 17:48:08 $

if nargin==0
	error('Incorrect number of arguments')
elseif nargin == 1
    radius = almanac('earth','radius','km');
elseif nargin == 2 & isstr(radius)
	radius = almanac(radius,'radius','km');
elseif nargin == 2 & ~isstr(radius)
    radius = sm2km(radius);
end


rad = km2rad(sm2km(sm),radius);