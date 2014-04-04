function sec=hms2sec(hms, in2, in3)

%HMS2SEC Converts time from hrs:min:sec to seconds
%
%  sec = HMS2SEC(hms) converts time from the hrs:min:sec vector format
%  to seconds.
%
%  sec = HMS2SEC(h,m,s) converts time from hours (h), minute (m) and
%  second (s) format to degrees.  The input matrices h, m and s must
%  be of equal size.  Minutes and seconds must be between 0 and 60.
%
%  See also SEC2HMS, HMS2HR, MAT2HMS, HMS2MAT, TIMEDIM, TIME2STR

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Brown, E. Byrns
%  $Revision: 1.7 $    $Date: 1998/08/10 17:47:48 $


if nargin == 0 | nargin == 2
	error('Incorrect number of arguments')
elseif nargin == 1
    [h,m,s] = hms2mat(hms);
elseif nargin == 3
    h = hms;   m = in2;    s = in3;
end

%  Compute the time in seconds by first transforming from hms to hrs.

sec = hr2sec(hms2hr(h,m,s));