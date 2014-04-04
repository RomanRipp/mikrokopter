function [map,maplegend] = zerom(latlim,lonlim,scale)

%ZEROM  Constructs a regular matrix map of all zeros.
%
%  map = ZEROM(latlim,lonlim,scale) constructs a regular matrix
%  map of all 0's.  The inputs are the latitude and longitude limits
%  as 2 element vectors, and a scale factor.  The latitude and longitude
%  limits are in degrees.  The scale factor represents the number of
%  matrix entries per single unit of latitude and longitude (eg:  10
%  entries per degree, 100 entries per degree).  The scale input must
%  be scalar.
%
%  [map,maplegend] = ZEROM(...) returns the maplegend vector for the
%  constructed regular matrix map.
%
%  See also NANM, SPZEROM, ONEM, ZEROS

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.7 $    $Date: 1998/08/10 17:48:14 $


if nargin ~= 3;  error('Incorrect number of arguments');  end


%  Determine the size of the desired map.
%  Then construct the appropriate zero map matrix


eval('[r,c,legend0] = sizem(latlim,lonlim,scale);','error(lasterr)')
map0 = zeros([r,c]);

%  Set the output arguments

if nargout == 1
    map = map0;

elseif nargout == 2
    map = map0;   maplegend = legend0;
end