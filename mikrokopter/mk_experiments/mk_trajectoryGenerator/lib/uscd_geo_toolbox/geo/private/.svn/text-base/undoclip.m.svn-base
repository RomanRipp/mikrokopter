function [lat,long] = undoclip(lat,long,splitpts,object)

%UNDOCLIP  Removes object clips introduced by CLIPDATA
%
%  [lat,long] = undoclip(lat,long,clippts,'object') will remove
%  the object clips introduced by CLIPDATA.  This function is necessary
%  to properly invert projected data from the cartesian space to the
%  original lat, long data points.  The input variable, clippts, must
%  be constructed by the function CLIPDATA.
%
%  Allowable object string are:  'surface' for undoing clipped graticules;
%  'light' for undoing clipped lights; 'line' for undoing clipped lines;
%  'patch' for undoing clipped patches; and 'text' for undoing clipped
%  text object location points.
%
%  See also CLIPDATA, TRIMDATA, UNDOTRIM

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.6 $    $Date: 1998/08/10 17:51:23 $


%  Argument tests

if nargin ~= 4;  error('Incorrect number of arguments');  end

%  Return if nothing to undo

if isempty(splitpts);    return;   end

%  Switch according to the correct object

switch object
    case 'surface'
	     lat(splitpts(:,1)) = splitpts(:,2);
	     long(splitpts(:,1)) = splitpts(:,3);

    case 'light'
         lat = lat;    long = long;

    case 'line'
	     lat(splitpts(:,1)) = [];
	     long(splitpts(:,1)) = [];

    case 'patch'           %  Simply replace the original patch data
         lat = splitpts(:,1);    long = splitpts(:,2);

    case 'text'
         lat = lat;    long = long;

    otherwise
         error(['Unrecognized object:  ',object])
end