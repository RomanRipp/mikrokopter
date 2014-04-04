function nm = km2nm(km)

%KM2NM Converts distances from kilometers to nautical miles
%
%  nm = KM2NM(km) converts distances from kilometers to nautical miles.
%
%  See also NM2KM, KM2DEG, KM2RAD, KM2SM, DISTDIM

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Brown, E. Byrns
%   $Revision: 1.7 $    $Date: 1998/08/10 17:47:51 $


if nargin==0
	error('Incorrect number of arguments')
elseif ~isreal(km)
     warning('Imaginary parts of complex DISTANCE argument ignored')
     km = real(km);
end


nm=km/1.8520;