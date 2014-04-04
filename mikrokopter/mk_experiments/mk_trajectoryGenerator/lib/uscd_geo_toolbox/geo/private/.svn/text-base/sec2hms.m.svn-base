function hms=sec2hms(sec)

%SEC2HMS Converts time from seconds to hrs:min:sec vector format
%
%  hms = SEC2HMS(sec) converts time from seconds to hrs:min:sec
%  vector format.
%
%  See also HMS2SEC, SEC2HR, MAT2HMS, HMS2MAT, TIMEDIM, TIME2STR

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Brown, E. Byrns
%  $Revision: 1.7 $    $Date: 1998/08/10 17:48:06 $


if nargin==0;   error('Incorrect number of arguments');   end

%  Compute the time in hms by first transforming from sec to hrs.

hms = hr2hms(sec2hr(sec));