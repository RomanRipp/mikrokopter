function eccent = axes2ecc(in1,in2)

%AXES2ECC  Computes eccentricity given semimajor and semiminor axes
%
%  e = AXES2ECC(semimajor,semiminor) computes the eccentricity
%  of an ellipse (or ellipsoid of revolution) given the semimajor
%  and semiminor axes.  The input data can be scalar or matrices
%  of equal dimensions.
%
%  e = AXES2ECC(vec) assumes a 2 element vector (vec) is supplied,
%  where vec = [semimajor semiminor].
%
%  See also  MAJAXIS, MINAXIS, ECC2FLAT, ECC2N

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.7 $    $Date: 1998/08/10 17:47:31 $


if nargin == 0
    error('Incorrect number of arguments')

elseif nargin == 1
    parm = in1;

    if ~isequal(sort(size(parm)),[1 2])
         error('Input must be a 2 element vector')
    else
         semimajor = max(parm);    semiminor = min(parm);
    end

elseif nargin == 2
    semimajor = in1;     semiminor = in2;

    if ~isequal(size(semimajor),size(semiminor))
         error('Inconsistent input dimensions')
	end
end

%  Ensure real inputs

if any([~isreal(semimajor) ~isreal(semiminor)])
   warning('Imaginary parts of complex arguments ignored')
   semimajor = real(semimajor);   semiminor = real(semiminor);
end

%  Compute the eccentricity

eccent = sqrt(semimajor.^2 - semiminor.^2) ./ semimajor;