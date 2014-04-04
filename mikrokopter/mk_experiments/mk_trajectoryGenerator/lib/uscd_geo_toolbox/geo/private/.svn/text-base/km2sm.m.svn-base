function sm = km2sm(km)

%KM2SM Converts distances from kilometers to statute miles
%
%  sm = KM2SM(km) converts distances from kilometers to statute miles.
%
%  See also SM2KM, KM2DEG, KM2RAD, KM2NM, DISTDIM

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Brown, E. Byrns
%   $Revision: 1.7 $    $Date: 1998/08/10 17:47:51 $

if nargin==0
	error('Incorrect number of arguments')
elseif ~isreal(km)
     warning('Imaginary parts of complex DISTANCE argument ignored')
     km = real(km);
end


sm=km/1.6093;