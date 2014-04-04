function nm = sm2nm(sm)

%SM2NM Converts distances from statute miles to nautical miles
%
%  nm = SM2NM(sm) converts distances from statute miles to nautical miles.
%
%  See also NM2SM, SM2DEG, SM2RAD, SM2KM, DISTDIM

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Brown, E. Byrns
%   $Revision: 1.7 $    $Date: 1998/08/10 17:48:08 $

if nargin==0
	error('Incorrect number of arguments')
elseif ~isreal(sm)
     warning('Imaginary parts of complex DISTANCE argument ignored')
     sm = real(sm);
end


nm=sm/1.150779;