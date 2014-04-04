function hrs = sec2hr(seconds)

%SEC2HR Converts time from seconds to hours
%
%  hr = SEC2HR(sec) converts time from seconds to hours.
%
%  See also HR2SEC, SEC2HMS, TIMEDIM, TIME2STR

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Brown, E. Byrns
%  $Revision: 1.7 $    $Date: 1998/08/10 17:48:06 $

if nargin==0
	error('Incorrect number of arguments')
elseif ~isreal(seconds)
     warning('Imaginary parts of complex TIME argument ignored')
     seconds = real(seconds);
end

hrs=seconds/3600;