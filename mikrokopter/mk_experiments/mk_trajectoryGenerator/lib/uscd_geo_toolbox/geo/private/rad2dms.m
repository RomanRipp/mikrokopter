function dms=rad2dms(rad)

%RAD2DMS Converts angles from radians to deg:min:sec vector format
%
%  dms = RAD2DMS(rad) converts angles from radians to deg:min:sec vector
%  format.
%
%  See also DMS2RAD, RAD2DEG, MAT2DMS, DMS2MAT, ANGLEDIM, ANGL2STR

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.7 $    $Date: 1998/08/10 17:48:01 $

if nargin==0;  error('Incorrect number of arguments');   end

%  Compute the angle in dms by first transforming from rad to deg.

dms = deg2dms(rad2deg(rad));
