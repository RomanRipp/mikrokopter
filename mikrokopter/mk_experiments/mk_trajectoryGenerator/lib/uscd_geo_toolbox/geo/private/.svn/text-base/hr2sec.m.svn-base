function seconds = hr2sec(hrs)

%HR2SEC Converts time from hours to seconds
%
%  sec = HR2SEC(hr) converts time from hours to seconds.
%
%  See also SEC2HR, HR2HMS, TIMEDIM, TIME2STR

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Brown, E. Byrns
%  $Revision: 1.7 $    $Date: 1998/08/10 17:47:49 $

if nargin==0
	error('Incorrect number of arguments')
elseif ~isreal(hrs)
     warning('Imaginary parts of complex TIME argument ignored')
     hrs = real(hrs);
end

seconds=hrs*3600;