function hm=sec2hm(sec)

%SEC2HM Converts time from seconds to hrs:min vector format
%
%  hm = SEC2HMS(sec) converts time from seconds to hrs:min vector format.
%
%  See also HR2HM, SEC2HMS, TIMEDIM, TIME2STR

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Brown, E. Byrns
%  $Revision: 1.7 $    $Date: 1998/08/10 17:48:06 $

if nargin==0;   error('Incorrect number of arguments');   end

%  Compute the time in hm.
%  0.2 is used to round seconds.  0.3+0.2 = 0.5 which will round up
%  to an additional minute.  0.29+0.2 = 0.49 which will stay at
%  the curren minute.

hms = round(sec2hms(sec)+0.2);
[h,m,s] = hms2mat(hms);
hm = mat2hms(h,m);         %  Round 60 minutes to 1 hour here