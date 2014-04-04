function dm=rad2dm(rad)

%RAD2DM Converts angles from radians to deg:min vector format
%
%  dm = RAD2DM(rad) converts angles from radians to deg:min vector
%  format.  The angle is rounded to the nearest minute.
%
%  See also DM2RAD, DMS2RAD, RAD2DMS, MAT2DMS, DMS2MAT, ANGLEDIM, ANGL2STR

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.7 $    $Date: 1998/08/10 17:48:01 $

if nargin==0;   error('Incorrect number of arguments');   end

%  Compute the angle in dm.
%  0.2 is used to round seconds.  0.3+0.2 = 0.5 which will round up
%  to an additional minute.  0.29+0.2 = 0.49 which will stay at
%  the curren minute.

dms = round(rad2dms(rad)+0.2);
[d,m,s] = dms2mat(dms);
dm = mat2dms(d,m);         %  Round 60 minutes to 1 degree here