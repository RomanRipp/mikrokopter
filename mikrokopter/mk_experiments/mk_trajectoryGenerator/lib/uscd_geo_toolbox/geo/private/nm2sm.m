function sm = nm2sm(nm)

%NM2SM Converts distances from nautical miles to statute miles
%
%  sm = NM2SM(nm) converts distances from nautical miles to statute miles.
%
%  See also SM2NM, NM2DEG, NM2RAD, NM2KM, DISTDIM

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Brown, E. Byrns
%   $Revision: 1.7 $    $Date: 1998/08/10 17:47:59 $

if nargin==0
	error('Incorrect number of arguments')
elseif ~isreal(nm)
     warning('Imaginary parts of complex distance argument ignored')
     nm = real(nm);
end


sm=nm*1.150779;