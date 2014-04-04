function km = sm2km(sm)

%SM2KM Converts distances from statute miles to kilometers
%
%  km = SM2KM(sm) converts distances from statute miles to kilometers.
%
%  See also KM2SM, SM2DEG, SM2RAD, SM2NM, DISTDIM

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Brown, E. Byrns
%   $Revision: 1.7 $    $Date: 1998/08/10 17:48:08 $

if nargin==0
	error('Incorrect number of arguments')
elseif ~isreal(sm)
     warning('Imaginary parts of complex DISTANCE argument ignored')
     sm = real(sm);
end


km = 1.6093*sm;