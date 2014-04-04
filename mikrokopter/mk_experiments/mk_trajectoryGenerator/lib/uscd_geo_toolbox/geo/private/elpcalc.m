function val = elpcalc(geoid,calculation)

%ELPCALC  Computes volume and surface area for an oblate spheroid
%
%  ELPCALC(geoid,'volume') computes the volume of an oblate
%  spheriod defined by geoid = [SemimajorAxis  Eccentricity].
%  The units are defined by the input semimajor axis.
%
%  ELPCALC(geoid,'surfarea') computes the surface area of an oblate
%  spheriod defined by geoid = [SemimajorAxis  Eccentricity].
%
%  See also SPHCALC


%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.7 $    $Date: 1998/08/10 17:47:40 $

if nargin ~= 2;    error('Incorrect number of arguments');   end

%  Test the geoid input

[geoid,msg] = geoidtst(geoid);
if ~isempty(msg);   error(msg);   end

semiminor =  minaxis(geoid);
semimajor = geoid(1);
eccent    = geoid(2);

if strcmp(calculation,'volume')
    val = (4*pi/3) * semiminor * semimajor^2;

elseif strcmp(calculation,'surfarea')
    if eccent > 1E-10
	    fact = log((1+eccent)/(1-eccent))/eccent;
	 else
	    fact = 2;
	 end

	 val = 2*pi*semimajor^2 + fact*pi*semiminor^2;

else
    error('Unrecognized calculation string')
end