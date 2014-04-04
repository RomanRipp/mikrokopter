function val = sphcalc(rad,calculation)

%SPHCALC  Computes volume and surface area for a sphere
%
%  SPHCALC(r,'volume') computes the volume of a sphere
%  defined by the input radius.  The units are defined
%  by the input radius.
%
%  SPHCALC(r,'surfarea') computes the surface area of a sphere
%  defined by the input radius.
%
%  See also ELPCALC

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.7 $    $Date: 1998/08/10 17:48:09 $


if nargin ~= 2
    error('Incorrect number of arguments')
elseif max(size(rad)) ~= 1
    error('Radius input must be a scalar')
end


if strcmp(calculation,'volume')
    val = (4*pi/3) * rad^3;

elseif strcmp(calculation,'surfarea')
    val = 4*pi*rad^2;

else
    error('Unrecognized calculation string')
end