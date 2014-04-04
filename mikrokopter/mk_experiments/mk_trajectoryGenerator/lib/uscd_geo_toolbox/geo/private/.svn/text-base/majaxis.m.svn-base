function semimajor = majaxis(in1,in2)

%MAJAXIS  Computes semimajor axis given a semiminor axis and an eccentricity
%
%  a = MAJAXIS(semiminor,e) computes the semimajor axis of an ellipse
%  (or ellipsoid of revolution) given the semiminor axis and eccentricity.
%  The input data can be scalar or matrices of equal dimensions.
%
%  a = MAJAXIS(vec) assumes a 2 element vector (vec) is supplied,
%  where vec = [semiminor, e].
%
%  See also:  AXES2ECC, FLAT2ECC, MINAXIS, N2ECC

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.7 $    $Date: 1998/08/10 17:47:53 $


if nargin == 0
    error('Incorrect number of arguments')

elseif nargin == 1
    if ~isequal(sort(size(in1)),[1 2])
         error('Input must be a 2 element vector')
    else
         semiminor = in1(1);    eccent    = in1(2);
    end

elseif nargin == 2
    if ~isequal(size(in1),size(in2))
         error('Inconsistent input dimensions')
	else
          semiminor = in1;    eccent  = in2;
	end
end

%  Ensure real inputs

if any([~isreal(eccent) ~isreal(semiminor)])
   warning('Imaginary parts of complex arguments ignored')
   eccent = real(eccent);   semiminor = real(semiminor);
end

%  Compute the semimajor axis

semimajor = semiminor ./ sqrt(1 - eccent.^2);