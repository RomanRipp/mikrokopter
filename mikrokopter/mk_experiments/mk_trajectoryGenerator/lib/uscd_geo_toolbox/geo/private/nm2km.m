function km = nm2km(nm)

%NM2KM Converts distances from nautical miles to kilometers
%
%  km = NM2KM(nm) converts distances from nautical miles to kilometers.
%
%  See also KM2NM, NM2DEG, NM2RAD, NM2SM, DISTDIM

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Brown, E. Byrns
%   $Revision: 1.7 $    $Date: 1998/08/10 17:47:59 $

if nargin==0
	error('Incorrect number of arguments')
elseif ~isreal(nm)
     warning('Imaginary parts of complex DISTANCE argument ignored')
     nm = real(nm);
end


km=1.8520*nm;