function mat = almanac(object,parm,units,refbody)

%ALMANAC  Almanac data for the Solar System
%
%  ALMANAC, without any input arguments, displays a list of the
%  heavenly objects recognized by ALMANAC.
%
%  ALMANAC('object') displays recognized parameters, units, and
%  reference body strings for the planet.
%
%  ALMANAC('object','parameter') returns the specified parameter from the
%  almanac.  Available parameters are the spherical radius of the planet,
%  surface area and volume of the sphere, the definition of the ellipsoid
%  (semimajor axis and eccentricity), the volume and surface area of the
%  ellipsoid, and tabulated surface area and volume.
%
%  ALMANAC('object','parameter','units') returns the corresponding parameter
%  in the units defined by 'units'.  If omitted, kilometers are used.
%
%  ALMANAC('object','parameter','units','referencebody') returns the
%  corresponding parameter given the spherical and or elliptical reference
%  body specified by 'referencebody'.  If omitted, a sphere is assumed where
%  appropriate.
%

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown, W. Stumpf
%   $Revision: 1.6 $    $Date: 1998/08/10 17:50:29 $

% List valid objects

if nargin == 0;

	disp(strvcat('The heavenly objects recognized by ALMANAC are:', ...
	'  ', ...
	'  sun', ...
	'  mercury',...
	'  venus',...
	'  earth',...
	'  moon',...
	'  mars',...
	'  jupiter',...
	'  saturn',...
	'  uranus', ...
	'  neptune',...
	'  pluto'...
	)   )

	return
end

% check for valid objects

if ~isstr(object); error('Object must be a string'); end

% handle the most common case for faster execution

if strcmp(object,'earth') & nargin == 2
	mat = earth(parm);
	return
end

% check for valid object string

indx = strmatch(lower(object),{...
			'sun',...
			'mercury',...
			'venus',...
			'earth',...
			'moon',...
			'mars',...
			'jupiter',...
			'saturn',...
			'uranus',...
			'neptune',...
			'pluto'},...
				'exact');

if isempty(indx); error('Unrecognized object string'); end

% call the appropriate local function

if nargin==1;
	feval(object);
elseif nargin==2
	mat = feval(object,parm);
elseif nargin==3
	mat = feval(object,parm,units);
elseif nargin==4
	mat = feval(object,parm,units,refbody);
else
	error('Incorrect number of arguments')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mat = earth(parm,units,refbody)

%EARTH  Calculates the almanac data for the planet Earth
%
%  EARTH, without any input arguments, displays recognized
%  parameters, units and reference body strings.
%
%  EARTH('parameter') returns the specified parameter from the almanac.
%  Available parameters are the spherical radius of the planet, surface
%  area and volume of the sphere, the definition of the ellipsoid
%  (semimajor axis and eccentricity), the volume and surface area of the
%  ellipsoid, and actual surface area and volume.
%
%  EARTH('parameter','units') returns the corresponding parameter in
%  the units defined by 'units'.  If omitted, kilometers are used.
%
%  EARTH('parameter','units','referencebody') returns the corresponding
%  parameter given the spherical and or elliptical reference body
%  specified by 'referencebody'.  If omitted, a sphere is assumed
%  where appropriate.
%
%  See also:  JUPITER, MARS, MERCURY, MOON, NEPTUNE, PLUTO, SATURN
%             SUN, URANUS, VENUS

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown, W. Stumpf
%   $Revision: 1.6 $    $Date: 1998/08/10 17:50:29 $

%  Define the default geoid as a string

defaultgeoid = 'grs80';

%  Define valid parameter strings.  Pre-padding is
%  much faster than strvcat

validparm = [
             'list    '
			 'radius  '
			 'geoid   '
			 'sphere  '
			 'surfarea'
			 'volume  '
		    ];

%  Define valid reference body strings

validref = [
            'sphere'
			'geoid '
			'actual'
		   ];

%  Define valid ellipsoid strings

validelps = [
             'everest      '
			 'bessel       '
			 'airy         '
			 'clarke66     '
			 'clarke80     '
			 'international'
			 'krasovsky    '
			 'wgs60        '
			 'iau65        '
			 'wgs66        '
			 'iau68        '
			 'wgs72        '
			 'grs80        '
			];

%  Test the input arguments

if nargin == 0;           parm    = [];     units   = [];     refbody = [];
    elseif nargin == 1;   units   = [];     refbody = [];
    elseif nargin == 2;   refbody = [];
end

%  Initialize output argument
%  Initialization needed for to avoid warning with isempty test at end

mat = [];

%  Empty units tests

if isempty(units);     units   = 'km';       end

%  Test if parameter input, then search for a match

if isempty(parm)
     parm = 'list';

elseif ~isstr(parm)
     error('Input argument PARAMETER must be a string')
else
     strmat = str2mat(validparm,validelps);
     strindx = strmatch(lower(parm),strmat);  %  String match

     if length(strindx) == 1
	       parm = deblank(strmat(strindx,:));
		   if strcmp(parm,'geoid') & isempty(refbody);  parm = defaultgeoid; end
	 else
	       error(['Unrecognized parameter:  ',parm])
     end
end

%  Test if parameter input, then search for a match

if isempty(refbody)
     strindx = strmatch(parm,validelps);  %  String match on ellipsoids
	 if length(strindx) == 1;   refbody = parm;
	    else;                   refbody = 'sphere';
	 end

elseif ~isstr(refbody)
     error('Input argument REFBODY must be a string')
else
     strmat = str2mat(validref,validelps);
     strindx = strmatch(lower(refbody),strmat);  %  String match

     if length(strindx) == 1
	       refbody = deblank(strmat(strindx,:));
		   if strcmp(refbody,'geoid');  refbody = defaultgeoid;  end
	 else
	       error(['Unrecognized reference body:  ',refbody])
     end
end

%  Parameters
%    radius, geoid definitions
%    actual volume,  actual surface area
%    spherical volume,  spherical surface area
%    geoid volume, geoid surface area

%  Basic surface area data taken from Encyclopaedia Britannica, 1995.
%  The earth radius is that of a sphere with a volume equivalent to that
%  of an ellipsoid with the equatorial and polar radii tabulated in the
%  Encyclopaedia Britannica. The volume is based on the same ellipsoid.

volume   = 1.0832e+12;     %  Earth volume in kilometers^3
surfarea = 510100000;      %  Earth surface area in kilometers^2
radius   = 6.3710e+03;     %  Earth radius in kilometers


%  Geoid definitions.
%  From:  D. H. Maling, Coordinate Systems and Map Projections, 2nd Edition
%         Pergamon Press, 1992, pp. 10-11, Table 1.01.  Semimajor axes
%         reported in kilometers.  Error in the Clarke1866 entry in this
%         reference.  Corrected using J.P. Synder, Map Projections A Working
%         Manual, U.S. Geological Survey Paper 1395, US Government Printing
%         Office, Washington, DC, 1987, Table 1, p. 12.  All flattening
%         entries in these references are transformed to eccentricities.

everest1830   = [6377.276345   flat2ecc(1/300.8017)];
bessel1841    = [6377.397155   flat2ecc(1/299.1528)];
airy1849      = [6377.563396   flat2ecc(1/299.3249)];
clarke1866    = [6378.206400   flat2ecc(1/294.97866982)];
clarke1880    = [6378.249170   flat2ecc(1/293.465)];
internatl1924 = [6378.388      flat2ecc(1/297.0)];
krasovsky1940 = [6378.245      flat2ecc(1/298.3)];
wgs60         = [6378.165      flat2ecc(1/298.3)];
iau65         = [6378.160      flat2ecc(1/298.25)];
wgs66         = [6378.145      flat2ecc(1/298.25)];
iau68         = [6378.16000    flat2ecc(1/298.2472)];
wgs72         = [6378.135      flat2ecc(1/298.26)];
grs80         = [6378.13700    flat2ecc(1/298.257222101)];



if strcmp(parm,'list')
    disp('Function Call:   mat = almanac(''earth'',''parameter'',''units'',''referencebody'')')
    disp(' ')
    disp('Valid parameter strings are:')
    disp('   ''radius''       for the planet radius')
    disp('   ''geoid''        for the planet geoid vector')
    disp('   ''volume''       for the planet volume')
    disp('   ''surfarea''     for the planet surface area')
	disp('   Or any valid ellipsoid definition string defined below')
    disp(' ')
    disp('Valid units strings are:')
    disp('   ''degrees''       or ''deg'' for degrees')
    disp('   ''kilometers''    or ''km''  for kilometers   (default)')
    disp('   ''nauticalmiles'' or ''nm''  for nautical miles')
    disp('   ''radians''       or ''rad'' for radians')
    disp('   ''statutemiles''  or ''sm''  for statute miles')
    disp(' ')
    disp('Valid reference body strings are:')
    disp('   ''sphere''       for a sphere (default)')
    disp('   ''geoid''        for a ellipsoid')
    disp('   ''actual''       for the tabulated volume and surface area of the planet')
	disp('   Or any valid ellipsoid definition string defined below')
    disp(' ')
    disp('Valid ellipsoid definition strings are:')
    disp('   ''everest''        for the 1830 Everest ellipsoid')
    disp('   ''bessel''         for the 1841 Bessel ellipsoid')
    disp('   ''airy''           for the 1849 Airy ellipsoid')
    disp('   ''clarke66''       for the 1866 Clarke ellipsoid')
    disp('   ''clarke80''       for the 1880 Clarke ellipsoid')
    disp('   ''international''  for the 1924 International ellipsoid')
    disp('   ''krasovsky''      for the 1940 Krasovsky ellipsoid')
    disp('   ''wgs60''          for the 1960 World Geodetic System ellipsoid')
    disp('   ''iau65''          for the 1965 International Astronomical Union ellipsoid')
    disp('   ''wgs66''          for the 1966 World Geodetic System ellipsoid')
    disp('   ''iau68''          for the 1968 International Astronomical Union ellipsoid')
    disp('   ''wgs72''          for the 1972 World Geodetic System ellipsoid')
    disp('   ''grs80''          for the 1980 Geodetic Reference System ellipsoid')
    disp('   An input of ''geoid'' is equivalent to ''grs80'' ')

    return
end

switch refbody
case 'actual'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'earth');
	  case 'geoid'
        eval(['mat = ',defaultgeoid,';']);
		mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'sphere'
	    mat = distdim(radius,'km',units,'earth');  mat = [mat 0];
	  case 'volume'
        fact = distdim(1.0,'km',units,'earth');
        mat = volume * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'earth');
        mat = surfarea * fact^2;
	end

case 'sphere'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'earth');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'earth');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'earth');  mat = [mat 0];
	  case 'volume'
	    rad = distdim(radius,'km',units,'earth');
        mat = sphcalc(rad,parm);
	  case 'surfarea'
	    rad = distdim(radius,'km',units,'earth');
        mat = sphcalc(rad,parm);
	end

case 'everest'
    switch parm
	  case 'radius'
        mat = everest1830;
		mat(1) = distdim(mat(1),'km',units,'earth');  mat(2) = [];
		warning('Semimajor axis returned for radius parameter')
	  case 'everest'
        mat = everest1830;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'geoid'
        mat = everest1830;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'volume'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(everest1830,parm);      mat = mat * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(everest1830,parm);      mat = mat * fact^2;
	end

case 'bessel'
    switch parm
	  case 'radius'
        mat = bessel1841;
		mat(1) = distdim(mat(1),'km',units,'earth');  mat(2) = [];
		warning('Semimajor axis returned for radius parameter')
	  case 'bessel'
        mat = bessel1841;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'geoid'
        mat = bessel1841;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'volume'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(bessel1841,parm);      mat = mat * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(bessel1841,parm);      mat = mat * fact^2;
	end

case 'airy'
    switch parm
	  case 'radius'
        mat = airy1849;
		mat(1) = distdim(mat(1),'km',units,'earth');  mat(2) = [];
		warning('Semimajor axis returned for radius parameter')
	  case 'airy'
        mat = airy1849;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'geoid'
        mat = airy1849;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'volume'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(airy1849,parm);      mat = mat * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(airy1849,parm);      mat = mat * fact^2;
	end

case 'clarke66'
    switch parm
	  case 'radius'
        mat = clarke1866;
		mat(1) = distdim(mat(1),'km',units,'earth');  mat(2) = [];
		warning('Semimajor axis returned for radius parameter')
	  case 'clarke66'
        mat = clarke1866;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'geoid'
        mat = clarke1866;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'volume'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(clarke1866,parm);      mat = mat * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(clarke1866,parm);      mat = mat * fact^2;
	end

case 'clarke80'
    switch parm
	  case 'radius'
        mat = clarke1880;
		mat(1) = distdim(mat(1),'km',units,'earth');  mat(2) = [];
		warning('Semimajor axis returned for radius parameter')
	  case 'clarke80'
        mat = clarke1880;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'geoid'
        mat = clarke1880;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'volume'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(clarke1880,parm);      mat = mat * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(clarke1880,parm);      mat = mat * fact^2;
	end

case 'international'
    switch parm
	  case 'radius'
        mat = internatl1924;
		mat(1) = distdim(mat(1),'km',units,'earth');  mat(2) = [];
		warning('Semimajor axis returned for radius parameter')
	  case 'international'
        mat = internatl1924;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'geoid'
        mat = internatl1924;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'volume'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(internatl1924,parm);      mat = mat * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(internatl1924,parm);      mat = mat * fact^2;
	end

case 'krasovsky'
    switch parm
	  case 'radius'
        mat = krasovsky1940;
		mat(1) = distdim(mat(1),'km',units,'earth');  mat(2) = [];
		warning('Semimajor axis returned for radius parameter')
	  case 'krasovsky'
        mat = krasovsky1940;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'geoid'
        mat = krasovsky1940;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'volume'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(krasovsky1940,parm);      mat = mat * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(krasovsky1940,parm);      mat = mat * fact^2;
	end

case 'wgs60'
    switch parm
	  case 'radius'
        mat = wgs60;
		mat(1) = distdim(mat(1),'km',units,'earth');  mat(2) = [];
		warning('Semimajor axis returned for radius parameter')
	  case 'wgs60'
        mat = wgs60;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'geoid'
        mat = wgs60;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'volume'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(wgs60,parm);      mat = mat * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(wgs60,parm);      mat = mat * fact^2;
	end

case 'iau65'
    switch parm
	  case 'radius'
        mat = iau65;
		mat(1) = distdim(mat(1),'km',units,'earth');  mat(2) = [];
		warning('Semimajor axis returned for radius parameter')
	  case 'iau65'
        mat = iau65;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'geoid'
        mat = iau65;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'volume'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(iau65,parm);      mat = mat * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(iau65,parm);      mat = mat * fact^2;
	end

case 'wgs66'
    switch parm
	  case 'radius'
        mat = wgs66;
		mat(1) = distdim(mat(1),'km',units,'earth');  mat(2) = [];
		warning('Semimajor axis returned for radius parameter')
	  case 'wgs66'
        mat = wgs66;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'geoid'
        mat = wgs66;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'volume'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(wgs66,parm);      mat = mat * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(wgs66,parm);      mat = mat * fact^2;
	end

case 'iau68'
    switch parm
	  case 'radius'
        mat = iau68;
		mat(1) = distdim(mat(1),'km',units,'earth');  mat(2) = [];
		warning('Semimajor axis returned for radius parameter')
	  case 'iau68'
        mat = iau68;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'geoid'
        mat = iau68;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'volume'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(iau68,parm);      mat = mat * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(iau68,parm);      mat = mat * fact^2;
	end

case 'wgs72'
    switch parm
	  case 'radius'
        mat = wgs72;
		mat(1) = distdim(mat(1),'km',units,'earth');  mat(2) = [];
		warning('Semimajor axis returned for radius parameter')
	  case 'wgs72'
        mat = wgs72;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'geoid'
        mat = wgs72;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'volume'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(wgs72,parm);      mat = mat * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(wgs72,parm);      mat = mat * fact^2;
	end

case 'grs80'
    switch parm
	  case 'radius'
        mat = grs80;
		mat(1) = distdim(mat(1),'km',units,'earth');  mat(2) = [];
		warning('Semimajor axis returned for radius parameter')
	  case 'grs80'
        mat = grs80;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'geoid'
        mat = grs80;    mat(1) = distdim(mat(1),'km',units,'earth');
	  case 'volume'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(grs80,parm);      mat = mat * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'earth');
        mat = elpcalc(grs80,parm);      mat = mat * fact^2;
	end

otherwise
    error('Unrecognized reference body string')
end

if isempty(mat)
    error('Unrecognized parameter and reference body combination')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mat = jupiter(parm,units,refbody)

%JUPITER  Calculates the almanac data for the planet Jupiter
%
%  JUPITER, without any input arguments, displays recognized
%  parameters, units and reference body strings.
%
%  JUPITER('parameter') returns the specified parameter from the almanac.
%  Available parameters are the spherical radius of the planet, surface
%  area and volume of the sphere, the definition of the ellipsoid
%  (semimajor axis and eccentricity), the volume and surface area of the
%  ellipsoid, and actual surface area and volume.
%
%  JUPITER('parameter','units') returns the corresponding parameter in
%  the units defined by 'units'.  If omitted, kilometers are used.
%
%  JUPITER('parameter','units','referencebody') returns the corresponding
%  parameter given the spherical and or elliptical reference body
%  specified by 'referencebody'.  If omitted, a sphere is assumed
%  where appropriate.
%
%  See also:  EARTH, MARS, MERCURY, MOON, NEPTUNE, PLUTO, SATURN
%             SUN, URANUS, VENUS

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown, W. Stumpf
%   $Revision: 1.6 $    $Date: 1998/08/10 17:50:29 $


%  Define valid parameter strings.  Pre-padding is
%  much faster than strvcat

validparm = [
             'list    '
			 'radius  '
			 'geoid   '
			 'sphere  '
			 'surfarea'
			 'volume  '
		    ];

%  Define valid reference body strings

validref = [
            'sphere'
			'geoid '
			'actual'
		   ];

%  Test the input arguments

if nargin == 0;           parm    = [];     units   = [];     refbody = [];
    elseif nargin == 1;   units   = [];     refbody = [];
    elseif nargin == 2;   refbody = [];
end

%  Initialize output argument
%  Initialization needed for to avoid warning with isempty test at end

mat = [];

%  Empty units tests

if isempty(units);     units   = 'km';       end

%  Test if parameter input, then search for a match

if isempty(parm)
     parm = 'list';
elseif ~isstr(parm)
     error('Input argument PARAMETER must be a string')
else
     strindx = strmatch(lower(parm),validparm);  %  String match
     if length(strindx) == 1
	       parm = deblank(validparm(strindx,:));
	 else
	       error(['Unrecognized parameter:  ',parm])
     end
end

%  Test if parameter input, then search for a match

if isempty(refbody)
	 if strcmp(parm,'geoid');   refbody = 'geoid';
	    else;                   refbody = 'sphere';
	 end
elseif ~isstr(refbody)
     error('Input argument REFBODY must be a string')
else
     strindx = strmatch(lower(refbody),validref);  %  String match
     if length(strindx) == 1
	       refbody = deblank(validref(strindx,:));
	 else
	       error(['Unrecognized reference body:  ',refbody])
     end
end

%  Parameters
%    radius, geoid definitions
%    actual volume,  actual surface area
%    spherical volume,  spherical surface area
%    geoid volume, geoid surface area

%  Geoid data derived from polar and equatorial radii tabulated in the
%  Encyclopaedia Britannica, 1995.  Radius is the radius of the sphere
%  with the same volume as the ellipsoid.  Volume and surface area
%  are derived from the geoid.

%  Geoid definitions.

geoid   = [71492    0.3574];

radius   = 6.9882e+04;           %  Jupiter radius in kilometers
volume   = 1.4295e+15;           %  Jupiter volume in kilometers^3
surfarea = 6.1419e+10;           %  Jupiter surface area in kilometers^2



if strcmp(parm,'list')
    disp('Function Call:   mat = almanac(''jupiter'',''parameter'',''units'',''referencebody'')')
    disp(' ')
    disp('Valid parameter strings are:')
    disp('   ''radius''       for the planet radius')
    disp('   ''geoid''        for the planet geoid vector')
    disp('   ''volume''       for the planet volume')
    disp('   ''surfarea''     for the planet surface area')
    disp(' ')
    disp('Valid units strings are:')
    disp('   ''degrees''       or ''deg'' for degrees')
    disp('   ''kilometers''    or ''km''  for kilometers   (default)')
    disp('   ''nauticalmiles'' or ''nm''  for nautical miles')
    disp('   ''radians''       or ''rad'' for radians')
    disp('   ''statutemiles''  or ''sm''  for statute miles')
    disp(' ')
    disp('Valid reference body strings are:')
    disp('   ''sphere''       for a sphere (default)')
    disp('   ''geoid''        for a ellipsoid')
    disp('   ''actual''       for the tabulated volume and surface area of the planet')


    return
end

switch refbody
case 'actual'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'jupiter');
	  case 'geoid'
        mat = geoid;    mat(1) = distdim(mat(1),'km',units,'jupiter');
	  case 'sphere'
	    mat = distdim(radius,'km',units,'jupiter');  mat = [mat 0];
	  case 'volume'
        fact = distdim(1.0,'km',units,'jupiter');
        mat = volume * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'jupiter');
        mat = surfarea * fact^2;
	end

case 'sphere'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'jupiter');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'jupiter');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'jupiter');  mat = [mat 0];
	  case 'volume'
	    rad = distdim(radius,'km',units,'jupiter');
        mat = sphcalc(rad,parm);
	  case 'surfarea'
	    rad = distdim(radius,'km',units,'jupiter');
        mat = sphcalc(rad,parm);
	end

case 'geoid'
    switch parm
	  case 'radius'
        mat = geoid;
		mat(1) = distdim(mat(1),'km',units,'jupiter');  mat(2) = [];
		warning('Semimajor axis returned for radius parameter')
	  case 'geoid'
        mat = geoid;    mat(1) = distdim(mat(1),'km',units,'jupiter');
	  case 'volume'
        fact = distdim(1.0,'km',units,'jupiter');
        mat = elpcalc(geoid,parm);      mat = mat * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'jupiter');
        mat = elpcalc(geoid,parm);      mat = mat * fact^2;
	end

otherwise
    error('Unrecognized reference body string')
end

if isempty(mat)
    error('Unrecognized parameter and reference body combination')
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mat = mars(parm,units,refbody)

%MARS  Calculates the almanac data for the planet mars
%
%  MARS, without any input arguments, displays recognized
%  parameters, units and reference body strings.
%
%  MARS('parameter') returns the specified parameter from the almanac.
%  Available parameters are the spherical radius of the planet, surface
%  area and volume of the sphere, the definition of the ellipsoid
%  (semimajor axis and eccentricity), the volume and surface area of the
%  ellipsoid, and actual surface area and volume.
%
%  MARS('parameter','units') returns the corresponding parameter in
%  the units defined by 'units'.  If omitted, kilometers are used.
%
%  MARS('parameter','units','referencebody') returns the corresponding
%  parameter given the spherical and or elliptical reference body
%  specified by 'referencebody'.  If omitted, a sphere is assumed
%  where appropriate.
%
%  See also:  EARTH, JUPITER, MERCURY, MOON, NEPTUNE, PLUTO, SATURN
%             SUN, URANUS, VENUS

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown, W. Stumpf
%   $Revision: 1.6 $    $Date: 1998/08/10 17:50:29 $


%  Define valid parameter strings.  Pre-padding is
%  much faster than strvcat

validparm = [
             'list    '
			 'radius  '
			 'geoid   '
			 'sphere  '
			 'surfarea'
			 'volume  '
		    ];

%  Define valid reference body strings

validref = [
            'sphere'
			'geoid '
			'actual'
		   ];

%  Test the input arguments

if nargin == 0;           parm    = [];     units   = [];     refbody = [];
    elseif nargin == 1;   units   = [];     refbody = [];
    elseif nargin == 2;   refbody = [];
end

%  Initialize output argument
%  Initialization needed for to avoid warning with isempty test at end

mat = [];

%  Empty units tests

if isempty(units);     units   = 'km';       end

%  Test if parameter input, then search for a match

if isempty(parm)
     parm = 'list';
elseif ~isstr(parm)
     error('Input argument PARAMETER must be a string')
else
     strindx = strmatch(lower(parm),validparm);  %  String match
     if length(strindx) == 1
	       parm = deblank(validparm(strindx,:));
	 else
	       error(['Unrecognized parameter:  ',parm])
     end
end

%  Test if parameter input, then search for a match

if isempty(refbody)
	 if strcmp(parm,'geoid');   refbody = 'geoid';
	    else;                   refbody = 'sphere';
	 end
elseif ~isstr(refbody)
     error('Input argument REFBODY must be a string')
else
     strindx = strmatch(lower(refbody),validref);  %  String match
     if length(strindx) == 1
	       refbody = deblank(validref(strindx,:));
	 else
	       error(['Unrecognized reference body:  ',refbody])
     end
end

%  Parameters
%    radius, geoid definitions
%    actual volume,  actual surface area
%    spherical volume,  spherical surface area
%    geoid volume, geoid surface area

%  Geoid data derived from polar and equatorial radii tabulated in the
%  Encyclopaedia Britannica, 1995.  Radius is the radius of the sphere
%  with the same volume as the ellipsoid.  Volume and surface area
%  are from the Encyclopaedia Britannica.

%  Geoid definitions.

geoid   = [3396.9    0.1105];

radius   = 3.39e+03;           %  Mars radius in kilometers
volume   = 1.63e+11;           %  Mars volume in kilometers^3
surfarea = 1.44e+08;           %  Mars surface area in kilometers^2



if strcmp(parm,'list')
    disp('Function Call:   mat = almanac(''mars'',''parameter'',''units'',''referencebody'')')
    disp(' ')
    disp('Valid parameter strings are:')
    disp('   ''radius''       for the planet radius')
    disp('   ''geoid''        for the planet geoid vector')
    disp('   ''volume''       for the planet volume')
    disp('   ''surfarea''     for the planet surface area')
    disp(' ')
    disp('Valid units strings are:')
    disp('   ''degrees''       or ''deg'' for degrees')
    disp('   ''kilometers''    or ''km''  for kilometers   (default)')
    disp('   ''nauticalmiles'' or ''nm''  for nautical miles')
    disp('   ''radians''       or ''rad'' for radians')
    disp('   ''statutemiles''  or ''sm''  for statute miles')
    disp(' ')
    disp('Valid reference body strings are:')
    disp('   ''sphere''       for a sphere (default)')
    disp('   ''geoid''        for a ellipsoid')
    disp('   ''actual''       for the tabulated volume and surface area of the planet')


    return
end

switch refbody
case 'actual'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'mars');
	  case 'geoid'
        mat = geoid;    mat(1) = distdim(mat(1),'km',units,'mars');
	  case 'sphere'
	    mat = distdim(radius,'km',units,'mars');  mat = [mat 0];
	  case 'volume'
        fact = distdim(1.0,'km',units,'mars');
        mat = volume * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'mars');
        mat = surfarea * fact^2;
	end

case 'sphere'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'mars');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'mars');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'mars');  mat = [mat 0];
	  case 'volume'
	    rad = distdim(radius,'km',units,'mars');
        mat = sphcalc(rad,parm);
	  case 'surfarea'
	    rad = distdim(radius,'km',units,'mars');
        mat = sphcalc(rad,parm);
	end

case 'geoid'
    switch parm
	  case 'radius'
        mat = geoid;
		mat(1) = distdim(mat(1),'km',units,'mars');  mat(2) = [];
		warning('Semimajor axis returned for radius parameter')
	  case 'geoid'
        mat = geoid;    mat(1) = distdim(mat(1),'km',units,'mars');
	  case 'volume'
        fact = distdim(1.0,'km',units,'mars');
        mat = elpcalc(geoid,parm);      mat = mat * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'mars');
        mat = elpcalc(geoid,parm);      mat = mat * fact^2;
	end

otherwise
    error('Unrecognized reference body string')
end

if isempty(mat)
    error('Unrecognized parameter and reference body combination')
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mat = mercury(parm,units,refbody)

%MERCURY  Calculates the almanac data for the planet Mercury
%
%  MERCURY, without any input arguments, displays recognized
%  parameters, units and reference body strings.
%
%  MERCURY('parameter') returns the specified parameter from the almanac.
%  Available parameters are the spherical radius of the planet, surface
%  area and volume of the sphere, the definition of the ellipsoid
%  (semimajor axis and eccentricity), the volume and surface area of the
%  ellipsoid, and actual surface area and volume.
%
%  MERCURY('parameter','units') returns the corresponding parameter in
%  the units defined by 'units'.  If omitted, kilometers are used.
%
%  MERCURY('parameter','units','referencebody') returns the corresponding
%  parameter given the spherical and or elliptical reference body
%  specified by 'referencebody'.  If omitted, a sphere is assumed
%  where appropriate.
%
%  See also:  EARTH, JUPITER, MARS, MOON, NEPTUNE, PLUTO, SATURN
%             SUN, URANUS, VENUS

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown, W. Stumpf
%   $Revision: 1.6 $    $Date: 1998/08/10 17:50:29 $


%  Define valid parameter strings.  Pre-padding is
%  much faster than strvcat

validparm = [
             'list    '
			 'radius  '
			 'geoid   '
			 'sphere  '
			 'surfarea'
			 'volume  '
		    ];

%  Define valid reference body strings

validref = [
            'sphere'
			'geoid '
			'actual'
		   ];

%  Test the input arguments

if nargin == 0;           parm    = [];     units   = [];     refbody = [];
    elseif nargin == 1;   units   = [];     refbody = [];
    elseif nargin == 2;   refbody = [];
end

%  Initialize output argument
%  Initialization needed for to avoid warning with isempty test at end

mat = [];

%  Empty units tests

if isempty(units);     units   = 'km';       end

%  Test if parameter input, then search for a match

if isempty(parm)
     parm = 'list';
elseif ~isstr(parm)
     error('Input argument PARAMETER must be a string')
else
     strindx = strmatch(lower(parm),validparm);  %  String match
     if length(strindx) == 1
	       parm = deblank(validparm(strindx,:));
	 else
	       error(['Unrecognized parameter:  ',parm])
     end
end

%  Test if parameter input, then search for a match

if isempty(refbody)
     refbody = 'sphere';
elseif ~isstr(refbody)
     error('Input argument REFBODY must be a string')
else
     strindx = strmatch(lower(refbody),validref);  %  String match
     if length(strindx) == 1
	       refbody = deblank(validref(strindx,:));
	 else
	       error(['Unrecognized reference body:  ',refbody])
     end
end


%  Parameters
%    radius, geoid definitions
%    actual volume,  actual surface area
%    spherical volume,  spherical surface area
%    geoid volume, geoid surface area

%  Basic data taken from J.P. Synder, Map Projections A Working
%  Manual, U.S. Geological Survey Paper 1395, US Government Printing
%  Office, Washington, DC, 1987, Table 2, p. 14.  Volume and
%  surface area computed assuming a perfect sphere. Radius agrees
%  with the value tabulated in the Encyclopaedia Britannica. 1995.

radius   = 2439.0;         %  Mercury radius in kilometers
volume   = 6.0775e+10;     %  Mercury volume in kilometers^3
surfarea = 7.4754e+07;     %  Mercury surface area in kilometers^2


if strcmp(parm,'list')
    disp('Function Call:   mat = almanac(''mercury'',''parameter'',''units'',''referencebody'')')
    disp(' ')
    disp('Valid parameter strings are:')
    disp('   ''radius''       for the planet radius')
    disp('   ''geoid''        for the planet geoid vector (e = 0 for mercury)')
    disp('   ''volume''       for the planet volume')
    disp('   ''surfarea''     for the planet surface area')
    disp(' ')
    disp('Valid units strings are:')
    disp('   ''degrees''       or ''deg'' for degrees')
    disp('   ''kilometers''    or ''km''  for kilometers   (default)')
    disp('   ''nauticalmiles'' or ''nm''  for nautical miles')
    disp('   ''radians''       or ''rad'' for radians')
    disp('   ''statutemiles''  or ''sm''  for statute miles')
    disp(' ')
    disp('Valid reference body strings are:')
    disp('   ''sphere''       for a sphere (default)')
    disp('   ''geoid''        for a ellipsoid (e = 0 for mercury)')
    disp('   ''actual''       for the tabulated volume and surface area of the planet')

    return
end

switch refbody
case 'actual'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'mercury');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'mercury');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'mercury');  mat = [mat 0];
	  case 'volume'
        fact = distdim(1.0,'km',units,'mercury');
        mat = volume * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'mercury');
        mat = surfarea * fact^2;
	end

case 'sphere'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'mercury');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'mercury');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'mercury');  mat = [mat 0];
	  case 'volume'
	    rad = distdim(radius,'km',units,'mercury');
        mat = sphcalc(rad,parm);
	  case 'surfarea'
	    rad = distdim(radius,'km',units,'mercury');
        mat = sphcalc(rad,parm);
	end

case 'geoid'
    warning('Spherical reference body only')

    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'mercury');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'mercury');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'mercury');  mat = [mat 0];
	  case 'volume'
	    rad = distdim(radius,'km',units,'mercury');
        mat = sphcalc(rad,parm);
	  case 'surfarea'
	    rad = distdim(radius,'km',units,'mercury');
        mat = sphcalc(rad,parm);
	end

otherwise
    error('Unrecognized reference body string')
end

if isempty(mat)
    error('Unrecognized parameter and reference body combination')
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mat = moon(parm,units,refbody)

%MOON  Calculates the almanac data for the Moon
%
%  MOON, without any input arguments, displays recognized
%  parameters, units and reference body strings.
%
%  MOON('parameter') returns the specified parameter from the almanac.
%  Available parameters are the spherical radius of the planet, surface
%  area and volume of the sphere, the definition of the ellipsoid
%  (semimajor axis and eccentricity), the volume and surface area of the
%  ellipsoid, and actual surface area and volume.
%
%  MOON('parameter','units') returns the corresponding parameter in
%  the units defined by 'units'.  If omitted, kilometers are used.
%
%  MOON('parameter','units','referencebody') returns the corresponding
%  parameter given the spherical and or elliptical reference body
%  specified by 'referencebody'.  If omitted, a sphere is assumed
%  where appropriate.
%
%  See also:  EARTH, JUPITER, MARS, MERCURY, NEPTUNE, PLUTO, SATURN
%             SUN, URANUS, VENUS

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown, W. Stumpf
%   $Revision: 1.6 $    $Date: 1998/08/10 17:50:29 $


%  Define valid parameter strings.  Pre-padding is
%  much faster than strvcat

validparm = [
             'list    '
			 'radius  '
			 'geoid   '
			 'sphere  '
			 'surfarea'
			 'volume  '
		    ];

%  Define valid reference body strings

validref = [
            'sphere'
			'geoid '
			'actual'
		   ];

%  Test the input arguments

if nargin == 0;           parm    = [];     units   = [];     refbody = [];
    elseif nargin == 1;   units   = [];     refbody = [];
    elseif nargin == 2;   refbody = [];
end

%  Initialize output argument
%  Initialization needed for to avoid warning with isempty test at end

mat = [];

%  Empty units tests

if isempty(units);     units   = 'km';       end

%  Test if parameter input, then search for a match

if isempty(parm)
     parm = 'list';
elseif ~isstr(parm)
     error('Input argument PARAMETER must be a string')
else
     strindx = strmatch(lower(parm),validparm);  %  String match
     if length(strindx) == 1
	       parm = deblank(validparm(strindx,:));
	 else
	       error(['Unrecognized parameter:  ',parm])
     end
end

%  Test if parameter input, then search for a match

if isempty(refbody)
     refbody = 'sphere';
elseif ~isstr(refbody)
     error('Input argument REFBODY must be a string')
else
     strindx = strmatch(lower(refbody),validref);  %  String match
     if length(strindx) == 1
	       refbody = deblank(validref(strindx,:));
	 else
	       error(['Unrecognized reference body:  ',refbody])
     end
end

%  Parameters
%    radius, geoid definitions
%    actual volume,  actual surface area
%    spherical volume,  spherical surface area
%    geoid volume, geoid surface area

%  Basic data taken from J.P. Synder, Map Projections A Working
%  Manual, U.S. Geological Survey Paper 1395, US Government Printing
%  Office, Washington, DC, 1987, Table 2, p. 14. Radius agrees
%  with the value tabulated in the Encyclopaedia Britannica. 1995.


radius   = 1738.0;         %  Moon radius in kilometers
volume   = 2.1991e+10;     %  Moon volume in kilometers^3
surfarea = 3.7959e+07;     %  Moon surface area in kilometers^2


if strcmp(parm,'list')
    disp('Function Call:   mat = almanac(''moon'',''parameter'',''units'',''referencebody'')')
    disp(' ')
    disp('Valid parameter strings are:')
    disp('   ''radius''       for the planet radius')
    disp('   ''geoid''        for the planet geoid vector (e = 0 for moon)')
    disp('   ''volume''       for the planet volume')
    disp('   ''surfarea''     for the planet surface area')
    disp(' ')
    disp('Valid units strings are:')
    disp('   ''degrees''       or ''deg'' for degrees')
    disp('   ''kilometers''    or ''km''  for kilometers   (default)')
    disp('   ''nauticalmiles'' or ''nm''  for nautical miles')
    disp('   ''radians''       or ''rad'' for radians')
    disp('   ''statutemiles''  or ''sm''  for statute miles')
    disp(' ')
    disp('Valid reference body strings are:')
    disp('   ''sphere''       for a sphere (default)')
    disp('   ''geoid''        for a ellipsoid (e = 0 for moon)')
    disp('   ''actual''       for the tabulated volume and surface area of the planet')

    return
end

switch refbody
case 'actual'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'moon');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'moon');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'moon');  mat = [mat 0];
	  case 'volume'
        fact = distdim(1.0,'km',units,'moon');
        mat = volume * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'moon');
        mat = surfarea * fact^2;
	end

case 'sphere'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'moon');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'moon');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'moon');  mat = [mat 0];
	  case 'volume'
	    rad = distdim(radius,'km',units,'moon');
        mat = sphcalc(rad,parm);
	  case 'surfarea'
	    rad = distdim(radius,'km',units,'moon');
        mat = sphcalc(rad,parm);
	end

case 'geoid'
    warning('Spherical reference body only')

    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'moon');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'moon');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'moon');  mat = [mat 0];
	  case 'volume'
	    rad = distdim(radius,'km',units,'moon');
        mat = sphcalc(rad,parm);
	  case 'surfarea'
	    rad = distdim(radius,'km',units,'moon');
        mat = sphcalc(rad,parm);
	end

otherwise
    error('Unrecognized reference body string')
end

if isempty(mat)
    error('Unrecognized parameter and reference body combination')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mat = neptune(parm,units,refbody)

%NEPTUNE  Calculates the almanac data for the planet Neptune
%
%  NEPTUNE, without any input arguments, displays recognized
%  parameters, units and reference body strings.
%
%  NEPTUNE('parameter') returns the specified parameter from the almanac.
%  Available parameters are the spherical radius of the planet, surface
%  area and volume of the sphere, the definition of the ellipsoid
%  (semimajor axis and eccentricity), the volume and surface area of the
%  ellipsoid, and actual surface area and volume.
%
%  NEPTUNE('parameter','units') returns the corresponding parameter in
%  the units defined by 'units'.  If omitted, kilometers are used.
%
%  NEPTUNE('parameter','units','referencebody') returns the corresponding
%  parameter given the spherical and or elliptical reference body
%  specified by 'referencebody'.  If omitted, a sphere is assumed
%  where appropriate.
%
%  See also:  EARTH, JUPITER, MARS, MERCURY, MOON, NEPTUNE, PLUTO,
%             SUN, URANUS, VENUS

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown, W. Stumpf
%   $Revision: 1.6 $    $Date: 1998/08/10 17:50:29 $


%  Define valid parameter strings.  Pre-padding is
%  much faster than strvcat

validparm = [
             'list    '
			 'radius  '
			 'geoid   '
			 'sphere  '
			 'surfarea'
			 'volume  '
		    ];

%  Define valid reference body strings

validref = [
            'sphere'
			'geoid '
			'actual'
		   ];

%  Test the input arguments

if nargin == 0;           parm    = [];     units   = [];     refbody = [];
    elseif nargin == 1;   units   = [];     refbody = [];
    elseif nargin == 2;   refbody = [];
end

%  Initialize output argument
%  Initialization needed for to avoid warning with isempty test at end

mat = [];

%  Empty units tests

if isempty(units);     units   = 'km';       end

%  Test if parameter input, then search for a match

if isempty(parm)
     parm = 'list';
elseif ~isstr(parm)
     error('Input argument PARAMETER must be a string')
else
     strindx = strmatch(lower(parm),validparm);  %  String match
     if length(strindx) == 1
	       parm = deblank(validparm(strindx,:));
	 else
	       error(['Unrecognized parameter:  ',parm])
     end
end

%  Test if parameter input, then search for a match

if isempty(refbody)
	 if strcmp(parm,'geoid');   refbody = 'geoid';
	    else;                   refbody = 'sphere';
	 end
elseif ~isstr(refbody)
     error('Input argument REFBODY must be a string')
else
     strindx = strmatch(lower(refbody),validref);  %  String match
     if length(strindx) == 1
	       refbody = deblank(validref(strindx,:));
	 else
	       error(['Unrecognized reference body:  ',refbody])
     end
end

%  Parameters
%    radius, geoid definitions
%    actual volume,  actual surface area
%    spherical volume,  spherical surface area
%    geoid volume, geoid surface area

%  Geoid data derived from 1 bar polar and equatorial radii tabulated in the
%  Encyclopaedia Britannica, 1995.  Radius is the radius of the sphere
%  with the same volume as the ellipsoid.  Volume and surface area
%  are derived from the geoid.

%  Geoid definitions.

geoid   = [24764    0.1843];

radius   = 2.4622e+04;           %  Neptune radius in kilometers
volume   = 6.2524e+13;           %  Neptune volume in kilometers^3
surfarea = 7.6185e+09;           %  Neptune surface area in kilometers^2



if strcmp(parm,'list')
    disp('Function Call:   mat = almanac(''neptune'',''parameter'',''units'',''referencebody'')')
    disp(' ')
    disp('Valid parameter strings are:')
    disp('   ''radius''       for the planet radius')
    disp('   ''geoid''        for the planet geoid vector')
    disp('   ''volume''       for the planet volume')
    disp('   ''surfarea''     for the planet surface area')
    disp(' ')
    disp('Valid units strings are:')
    disp('   ''degrees''       or ''deg'' for degrees')
    disp('   ''kilometers''    or ''km''  for kilometers   (default)')
    disp('   ''nauticalmiles'' or ''nm''  for nautical miles')
    disp('   ''radians''       or ''rad'' for radians')
    disp('   ''statutemiles''  or ''sm''  for statute miles')
    disp(' ')
    disp('Valid reference body strings are:')
    disp('   ''sphere''       for a sphere (default)')
    disp('   ''geoid''        for a ellipsoid')
    disp('   ''actual''       for the tabulated volume and surface area of the planet')


    return
end

switch refbody
case 'actual'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'neptune');
	  case 'geoid'
        mat = geoid;    mat(1) = distdim(mat(1),'km',units,'neptune');
	  case 'sphere'
	    mat = distdim(radius,'km',units,'neptune');  mat = [mat 0];
	  case 'volume'
        fact = distdim(1.0,'km',units,'neptune');
        mat = volume * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'neptune');
        mat = surfarea * fact^2;
	end

case 'sphere'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'neptune');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'neptune');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'neptune');  mat = [mat 0];
	  case 'volume'
	    rad = distdim(radius,'km',units,'neptune');
        mat = sphcalc(rad,parm);
	  case 'surfarea'
	    rad = distdim(radius,'km',units,'neptune');
        mat = sphcalc(rad,parm);
	end

case 'geoid'
    switch parm
	  case 'radius'
        mat = geoid;
		mat(1) = distdim(mat(1),'km',units,'neptune');  mat(2) = [];
		warning('Semimajor axis returned for radius parameter')
	  case 'geoid'
        mat = geoid;    mat(1) = distdim(mat(1),'km',units,'neptune');
	  case 'volume'
        fact = distdim(1.0,'km',units,'neptune');
        mat = elpcalc(geoid,parm);      mat = mat * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'neptune');
        mat = elpcalc(geoid,parm);      mat = mat * fact^2;
	end

otherwise
    error('Unrecognized reference body string')
end

if isempty(mat)
    error('Unrecognized parameter and reference body combination')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mat = pluto(parm,units,refbody)

%PLUTO  Calculates the almanac data for the planet Pluto
%
%  PLUTO, without any input arguments, displays recognized
%  parameters, units and reference body strings.
%
%  PLUTO('parameter') returns the specified parameter from the almanac.
%  Available parameters are the spherical radius of the planet, surface
%  area and volume of the sphere, the definition of the ellipsoid
%  (semimajor axis and eccentricity), the volume and surface area of the
%  ellipsoid, and actual surface area and volume.
%
%  PLUTO('parameter','units') returns the corresponding parameter in
%  the units defined by 'units'.  If omitted, kilometers are used.
%
%  PLUTO('parameter','units','referencebody') returns the corresponding
%  parameter given the spherical and or elliptical reference body
%  specified by 'referencebody'.  If omitted, a sphere is assumed
%  where appropriate.
%
%  See also:  EARTH, JUPITER, MARS, MERCURY, MOON, NEPTUNE, SATURN
%             SUN, URANUS, VENUS

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown, W. Stumpf
%   $Revision: 1.6 $    $Date: 1998/08/10 17:50:29 $


%  Define valid parameter strings.  Pre-padding is
%  much faster than strvcat

validparm = [
             'list    '
			 'radius  '
			 'geoid   '
			 'sphere  '
			 'surfarea'
			 'volume  '
		    ];

%  Define valid reference body strings

validref = [
            'sphere'
			'geoid '
			'actual'
		   ];

%  Test the input arguments

if nargin == 0;           parm    = [];     units   = [];     refbody = [];
    elseif nargin == 1;   units   = [];     refbody = [];
    elseif nargin == 2;   refbody = [];
end

%  Initialize output argument
%  Initialization needed for to avoid warning with isempty test at end

mat = [];

%  Empty units tests

if isempty(units);     units   = 'km';       end

%  Test if parameter input, then search for a match

if isempty(parm)
     parm = 'list';
elseif ~isstr(parm)
     error('Input argument PARAMETER must be a string')
else
     strindx = strmatch(lower(parm),validparm);  %  String match
     if length(strindx) == 1
	       parm = deblank(validparm(strindx,:));
	 else
	       error(['Unrecognized parameter:  ',parm])
     end
end

%  Test if parameter input, then search for a match

if isempty(refbody)
     refbody = 'sphere';
elseif ~isstr(refbody)
     error('Input argument REFBODY must be a string')
else
     strindx = strmatch(lower(refbody),validref);  %  String match
     if length(strindx) == 1
	       refbody = deblank(validref(strindx,:));
	 else
	       error(['Unrecognized reference body:  ',refbody])
     end
end

%  Parameters
%    radius, geoid definitions
%    actual volume,  actual surface area
%    spherical volume,  spherical surface area
%    geoid volume, geoid surface area

%  Radius is the value tabulated in the Encyclopaedia Britannica, 1995.
%  Volume and surface area are derived for a sphere.

radius   = 1151;            %  Pluto radius in kilometers
volume   = 6.3873e+09;      %  Pluto volume in kilometers^3
surfarea = 1.6648e+07;      %  Pluto surface area in kilometers^2


if strcmp(parm,'list')
    disp('Function Call:   mat = almanac(''pluto'',''parameter'',''units'',''referencebody'')')
    disp(' ')
    disp('Valid parameter strings are:')
    disp('   ''radius''       for the planet radius')
    disp('   ''geoid''        for the planet geoid vector (e = 0 for pluto)')
    disp('   ''volume''       for the planet volume')
    disp('   ''surfarea''     for the planet surface area')
    disp(' ')
    disp('Valid units strings are:')
    disp('   ''degrees''       or ''deg'' for degrees')
    disp('   ''kilometers''    or ''km''  for kilometers   (default)')
    disp('   ''nauticalmiles'' or ''nm''  for nautical miles')
    disp('   ''radians''       or ''rad'' for radians')
    disp('   ''statutemiles''  or ''sm''  for statute miles')
    disp(' ')
    disp('Valid reference body strings are:')
    disp('   ''sphere''       for a sphere (default)')
    disp('   ''geoid''        for a ellipsoid (e = 0 for pluto)')
    disp('   ''actual''       for the tabulated volume and surface area of the planet')

    return
end

switch refbody
case 'actual'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'pluto');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'pluto');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'pluto');  mat = [mat 0];
	  case 'volume'
        fact = distdim(1.0,'km',units,'pluto');
        mat = volume * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'pluto');
        mat = surfarea * fact^2;
	end

case 'sphere'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'pluto');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'pluto');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'pluto');  mat = [mat 0];
	  case 'volume'
	    rad = distdim(radius,'km',units,'pluto');
        mat = sphcalc(rad,parm);
	  case 'surfarea'
	    rad = distdim(radius,'km',units,'pluto');
        mat = sphcalc(rad,parm);
	end

case 'geoid'
    warning('Spherical reference body only')

    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'pluto');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'pluto');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'pluto');  mat = [mat 0];
	  case 'volume'
	    rad = distdim(radius,'km',units,'pluto');
        mat = sphcalc(rad,parm);
	  case 'surfarea'
	    rad = distdim(radius,'km',units,'pluto');
        mat = sphcalc(rad,parm);
	end

otherwise
    error('Unrecognized reference body string')
end

if isempty(mat)
    error('Unrecognized parameter and reference body combination')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mat = saturn(parm,units,refbody)

%SATURN  Calculates the almanac data for the planet Saturn
%
%  SATURN, without any input arguments, displays recognized
%  parameters, units and reference body strings.
%
%  SATURN('parameter') returns the specified parameter from the almanac.
%  Available parameters are the spherical radius of the planet, surface
%  area and volume of the sphere, the definition of the ellipsoid
%  (semimajor axis and eccentricity), the volume and surface area of the
%  ellipsoid, and actual surface area and volume.
%
%  SATURN('parameter','units') returns the corresponding parameter in
%  the units defined by 'units'.  If omitted, kilometers are used.
%
%  SATURN('parameter','units','referencebody') returns the corresponding
%  parameter given the spherical and or elliptical reference body
%  specified by 'referencebody'.  If omitted, a sphere is assumed
%  where appropriate.
%
%  See also:  EARTH, JUPITER, MARS, MERCURY, MOON, NEPTUNE, PLUTO,
%             SUN, URANUS, VENUS

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown, W. Stumpf
%   $Revision: 1.6 $    $Date: 1998/08/10 17:50:29 $


%  Define valid parameter strings.  Pre-padding is
%  much faster than strvcat

validparm = [
             'list    '
			 'radius  '
			 'geoid   '
			 'sphere  '
			 'surfarea'
			 'volume  '
		    ];

%  Define valid reference body strings

validref = [
            'sphere'
			'geoid '
			'actual'
		   ];

%  Test the input arguments

if nargin == 0;           parm    = [];     units   = [];     refbody = [];
    elseif nargin == 1;   units   = [];     refbody = [];
    elseif nargin == 2;   refbody = [];
end

%  Initialize output argument
%  Initialization needed for to avoid warning with isempty test at end

mat = [];

%  Empty units tests

if isempty(units);     units   = 'km';       end

%  Test if parameter input, then search for a match

if isempty(parm)
     parm = 'list';
elseif ~isstr(parm)
     error('Input argument PARAMETER must be a string')
else
     strindx = strmatch(lower(parm),validparm);  %  String match
     if length(strindx) == 1
	       parm = deblank(validparm(strindx,:));
	 else
	       error(['Unrecognized parameter:  ',parm])
     end
end

%  Test if parameter input, then search for a match

if isempty(refbody)
	 if strcmp(parm,'geoid');   refbody = 'geoid';
	    else;                   refbody = 'sphere';
	 end
elseif ~isstr(refbody)
     error('Input argument REFBODY must be a string')
else
     strindx = strmatch(lower(refbody),validref);  %  String match
     if length(strindx) == 1
	       refbody = deblank(validref(strindx,:));
	 else
	       error(['Unrecognized reference body:  ',refbody])
     end
end

%  Parameters
%    radius, geoid definitions
%    actual volume,  actual surface area
%    spherical volume,  spherical surface area
%    geoid volume, geoid surface area

%  Geoid data derived from 1 bar polar and equatorial radii tabulated in the
%  Encyclopaedia Britannica, 1995.  Radius is the radius of the sphere
%  with the same volume as the ellipsoid.  Volume and surface area
%  are derived from the geoid.

%  Geoid definitions.

geoid   = [60268    0.4317];

radius   = 5.8235e+04;           %  Saturn radius in kilometers
volume   = 8.2711e+14;           %  Saturn volume in kilometers^3
surfarea = 4.2693e+10;           %  Saturn surface area in kilometers^2



if strcmp(parm,'list')
    disp('Function Call:   mat = almanac(''saturn'',''parameter'',''units'',''referencebody'')')
    disp(' ')
    disp('Valid parameter strings are:')
    disp('   ''radius''       for the planet radius')
    disp('   ''geoid''        for the planet geoid vector')
    disp('   ''volume''       for the planet volume')
    disp('   ''surfarea''     for the planet surface area')
    disp(' ')
    disp('Valid units strings are:')
    disp('   ''degrees''       or ''deg'' for degrees')
    disp('   ''kilometers''    or ''km''  for kilometers   (default)')
    disp('   ''nauticalmiles'' or ''nm''  for nautical miles')
    disp('   ''radians''       or ''rad'' for radians')
    disp('   ''statutemiles''  or ''sm''  for statute miles')
    disp(' ')
    disp('Valid reference body strings are:')
    disp('   ''sphere''       for a sphere (default)')
    disp('   ''geoid''        for a ellipsoid')
    disp('   ''actual''       for the tabulated volume and surface area of the planet')


    return
end

switch refbody
case 'actual'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'saturn');
	  case 'geoid'
        mat = geoid;    mat(1) = distdim(mat(1),'km',units,'saturn');
	  case 'sphere'
	    mat = distdim(radius,'km',units,'saturn');  mat = [mat 0];
	  case 'volume'
        fact = distdim(1.0,'km',units,'saturn');
        mat = volume * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'saturn');
        mat = surfarea * fact^2;
	end

case 'sphere'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'saturn');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'saturn');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'saturn');  mat = [mat 0];
	  case 'volume'
	    rad = distdim(radius,'km',units,'saturn');
        mat = sphcalc(rad,parm);
	  case 'surfarea'
	    rad = distdim(radius,'km',units,'saturn');
        mat = sphcalc(rad,parm);
	end

case 'geoid'
    switch parm
	  case 'radius'
        mat = geoid;
		mat(1) = distdim(mat(1),'km',units,'saturn');  mat(2) = [];
		warning('Semimajor axis returned for radius parameter')
	  case 'geoid'
        mat = geoid;    mat(1) = distdim(mat(1),'km',units,'saturn');
	  case 'volume'
        fact = distdim(1.0,'km',units,'saturn');
        mat = elpcalc(geoid,parm);      mat = mat * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'saturn');
        mat = elpcalc(geoid,parm);      mat = mat * fact^2;
	end

otherwise
    error('Unrecognized reference body string')
end

if isempty(mat)
    error('Unrecognized parameter and reference body combination')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mat = sun(parm,units,refbody)

%SUN  Calculates the almanac data for the Sun
%
%  SUN, without any input arguments, displays recognized
%  parameters, units and reference body strings.
%
%  SUN('parameter') returns the specified parameter from the almanac.
%  Available parameters are the spherical radius of the planet, surface
%  area and volume of the sphere, the definition of the ellipsoid
%  (semimajor axis and eccentricity), the volume and surface area of the
%  ellipsoid, and actual surface area and volume.
%
%  SUN('parameter','units') returns the corresponding parameter in
%  the units defined by 'units'.  If omitted, kilometers are used.
%
%  SUN('parameter','units','referencebody') returns the corresponding
%  parameter given the spherical and or elliptical reference body
%  specified by 'referencebody'.  If omitted, a sphere is assumed
%  where appropriate.
%
%  See also:  EARTH, JUPITER, MARS, MERCURY, MOON, NEUPTUNE, PLUTO
%             SATURN, URANUS, VENUS

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown, W. Stumpf
%   $Revision: 1.6 $    $Date: 1998/08/10 17:50:29 $


%  Define valid parameter strings.  Pre-padding is
%  much faster than strvcat

validparm = [
             'list    '
			 'radius  '
			 'geoid   '
			 'sphere  '
			 'surfarea'
			 'volume  '
		    ];

%  Define valid reference body strings

validref = [
            'sphere'
			'geoid '
			'actual'
		   ];

%  Test the input arguments

if nargin == 0;           parm    = [];     units   = [];     refbody = [];
    elseif nargin == 1;   units   = [];     refbody = [];
    elseif nargin == 2;   refbody = [];
end

%  Initialize output argument
%  Initialization needed for to avoid warning with isempty test at end

mat = [];

%  Empty units tests

if isempty(units);     units   = 'km';       end

%  Test if parameter input, then search for a match

if isempty(parm)
     parm = 'list';
elseif ~isstr(parm)
     error('Input argument PARAMETER must be a string')
else
     strindx = strmatch(lower(parm),validparm);  %  String match
     if length(strindx) == 1
	       parm = deblank(validparm(strindx,:));
	 else
	       error(['Unrecognized parameter:  ',parm])
     end
end

%  Test if parameter input, then search for a match

if isempty(refbody)
     refbody = 'sphere';
elseif ~isstr(refbody)
     error('Input argument REFBODY must be a string')
else
     strindx = strmatch(lower(refbody),validref);  %  String match
     if length(strindx) == 1
	       refbody = deblank(validref(strindx,:));
	 else
	       error(['Unrecognized reference body:  ',refbody])
     end
end

%  Parameters
%    radius, geoid definitions
%    actual volume,  actual surface area
%    spherical volume,  spherical surface area
%    geoid volume, geoid surface area

%  Basic radius data taken from Encyclopaedia Britannica, 1995.
%  Volume and surface area computed assuming a perfect sphere.

radius   = 6.9446e+05;         %  Sun radius in kilometers
volume   = 1.4029e+18;    %  Sun volume in kilometers^3
surfarea = 6.0604e+12;    %  Sun surface area in kilometers^2


if strcmp(parm,'list')
    disp('Function Call:   mat = almanac(''sun'',''parameter'',''units'',''referencebody'')')
    disp(' ')
    disp('Valid parameter strings are:')
    disp('   ''radius''       for the planet radius')
    disp('   ''geoid''        for the planet geoid vector (e = 0 for sun)')
    disp('   ''volume''       for the planet volume')
    disp('   ''surfarea''     for the planet surface area')
    disp(' ')
    disp('Valid units strings are:')
    disp('   ''degrees''       or ''deg'' for degrees')
    disp('   ''kilometers''    or ''km''  for kilometers   (default)')
    disp('   ''nauticalmiles'' or ''nm''  for nautical miles')
    disp('   ''radians''       or ''rad'' for radians')
    disp('   ''statutemiles''  or ''sm''  for statute miles')
    disp(' ')
    disp('Valid reference body strings are:')
    disp('   ''sphere''       for a sphere (default)')
    disp('   ''geoid''        for a ellipsoid (e = 0 for sun)')
    disp('   ''actual''       for the tabulated volume and surface area of the planet')

    return
end

switch refbody
case 'actual'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'sun');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'sun');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'sun');  mat = [mat 0];
	  case 'volume'
        fact = distdim(1.0,'km',units,'sun');
        mat = volume * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'sun');
        mat = surfarea * fact^2;
	end

case 'sphere'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'sun');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'sun');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'sun');  mat = [mat 0];
	  case 'volume'
	    rad = distdim(radius,'km',units,'sun');
        mat = sphcalc(rad,parm);
	  case 'surfarea'
	    rad = distdim(radius,'km',units,'sun');
        mat = sphcalc(rad,parm);
	end

case 'geoid'
    warning('Spherical reference body only')

    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'sun');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'sun');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'sun');  mat = [mat 0];
	  case 'volume'
	    rad = distdim(radius,'km',units,'sun');
        mat = sphcalc(rad,parm);
	  case 'surfarea'
	    rad = distdim(radius,'km',units,'sun');
        mat = sphcalc(rad,parm);
	end

otherwise
    error('Unrecognized reference body string')
end

if isempty(mat)
    error('Unrecognized parameter and reference body combination')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mat = uranus(parm,units,refbody)

%URANUS  Calculates the almanac data for the planet Uranus
%
%  URANUS, without any input arguments, displays recognized
%  parameters, units and reference body strings.
%
%  URANUS('parameter') returns the specified parameter from the almanac.
%  Available parameters are the spherical radius of the planet, surface
%  area and volume of the sphere, the definition of the ellipsoid
%  (semimajor axis and eccentricity), the volume and surface area of the
%  ellipsoid, and actual surface area and volume.
%
%  URANUS('parameter','units') returns the corresponding parameter in
%  the units defined by 'units'.  If omitted, kilometers are used.
%
%  URANUS('parameter','units','referencebody') returns the corresponding
%  parameter given the spherical and or elliptical reference body
%  specified by 'referencebody'.  If omitted, a sphere is assumed
%  where appropriate.
%
%  See also:  EARTH, JUPITER, MARS, MERCURY, MOON, NEPTUNE, PLUTO, SATURN
%             SUN, VENUS

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown, W. Stumpf
%   $Revision: 1.6 $    $Date: 1998/08/10 17:50:29 $


%  Define valid parameter strings.  Pre-padding is
%  much faster than strvcat

validparm = [
             'list    '
			 'radius  '
			 'geoid   '
			 'sphere  '
			 'surfarea'
			 'volume  '
		    ];

%  Define valid reference body strings

validref = [
            'sphere'
			'geoid '
			'actual'
		   ];

%  Test the input arguments

if nargin == 0;           parm    = [];     units   = [];     refbody = [];
    elseif nargin == 1;   units   = [];     refbody = [];
    elseif nargin == 2;   refbody = [];
end

%  Initialize output argument
%  Initialization needed for to avoid warning with isempty test at end

mat = [];

%  Empty units tests

if isempty(units);     units   = 'km';       end

%  Test if parameter input, then search for a match

if isempty(parm)
     parm = 'list';
elseif ~isstr(parm)
     error('Input argument PARAMETER must be a string')
else
     strindx = strmatch(lower(parm),validparm);  %  String match
     if length(strindx) == 1
	       parm = deblank(validparm(strindx,:));
	 else
	       error(['Unrecognized parameter:  ',parm])
     end
end

%  Test if parameter input, then search for a match

if isempty(refbody)
	 if strcmp(parm,'geoid');   refbody = 'geoid';
	    else;                   refbody = 'sphere';
	 end
elseif ~isstr(refbody)
     error('Input argument REFBODY must be a string')
else
     strindx = strmatch(lower(refbody),validref);  %  String match
     if length(strindx) == 1
	       refbody = deblank(validref(strindx,:));
	 else
	       error(['Unrecognized reference body:  ',refbody])
     end
end

%  Parameters
%    radius, geoid definitions
%    actual volume,  actual surface area
%    spherical volume,  spherical surface area
%    geoid volume, geoid surface area

%  Geoid data derived from equatorial radius (1 bar) and elipticity (flattening)
%  tabulated in the Encyclopaedia Britannica, 1995.  Radius is the radius
%  of the sphere with the same volume as the ellipsoid.  Volume and surface
%  area are derived from the geoid.

%  Geoid definitions.

geoid   = [25559      flat2ecc(0.0229)];

radius   = 2.5362e+04;           %  Uranus radius in kilometers
volume   = 6.8338e+13;           %  Uranus volume in kilometers^3
surfarea = 8.0841e+09;           %  Uranus surface area in kilometers^2



if strcmp(parm,'list')
    disp('Function Call:   mat = almanac(''uranus'',''parameter'',''units'',''referencebody'')')
    disp(' ')
    disp('Valid parameter strings are:')
    disp('   ''radius''       for the planet radius')
    disp('   ''geoid''        for the planet geoid vector')
    disp('   ''volume''       for the planet volume')
    disp('   ''surfarea''     for the planet surface area')
    disp(' ')
    disp('Valid units strings are:')
    disp('   ''degrees''       or ''deg'' for degrees')
    disp('   ''kilometers''    or ''km''  for kilometers   (default)')
    disp('   ''nauticalmiles'' or ''nm''  for nautical miles')
    disp('   ''radians''       or ''rad'' for radians')
    disp('   ''statutemiles''  or ''sm''  for statute miles')
    disp(' ')
    disp('Valid reference body strings are:')
    disp('   ''sphere''       for a sphere (default)')
    disp('   ''geoid''        for a ellipsoid')
    disp('   ''actual''       for the tabulated volume and surface area of the planet')


    return
end

switch refbody
case 'actual'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'uranus');
	  case 'geoid'
        mat = geoid;    mat(1) = distdim(mat(1),'km',units,'uranus');
	  case 'sphere'
	    mat = distdim(radius,'km',units,'uranus');  mat = [mat 0];
	  case 'volume'
        fact = distdim(1.0,'km',units,'uranus');
        mat = volume * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'uranus');
        mat = surfarea * fact^2;
	end

case 'sphere'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'uranus');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'uranus');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'uranus');  mat = [mat 0];
	  case 'volume'
	    rad = distdim(radius,'km',units,'uranus');
        mat = sphcalc(rad,parm);
	  case 'surfarea'
	    rad = distdim(radius,'km',units,'uranus');
        mat = sphcalc(rad,parm);
	end

case 'geoid'
    switch parm
	  case 'radius'
        mat = geoid;
		mat(1) = distdim(mat(1),'km',units,'uranus');  mat(2) = [];
		warning('Semimajor axis returned for radius parameter')
	  case 'geoid'
        mat = geoid;    mat(1) = distdim(mat(1),'km',units,'uranus');
	  case 'volume'
        fact = distdim(1.0,'km',units,'uranus');
        mat = elpcalc(geoid,parm);      mat = mat * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'uranus');
        mat = elpcalc(geoid,parm);      mat = mat * fact^2;
	end

otherwise
    error('Unrecognized reference body string')
end

if isempty(mat)
    error('Unrecognized parameter and reference body combination')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mat = venus(parm,units,refbody)

%VENUS  Calculates the almanac data for the planet venus
%
%  VENUS, without any input arguments, displays recognized
%  parameters, units and reference body strings.
%
%  VENUS('parameter') returns the specified parameter from the almanac.
%  Available parameters are the spherical radius of the planet, surface
%  area and volume of the sphere, the definition of the ellipsoid
%  (semimajor axis and eccentricity), the volume and surface area of the
%  ellipsoid, and actual surface area and volume.
%
%  VENUS('parameter','units') returns the corresponding parameter in
%  the units defined by 'units'.  If omitted, kilometers are used.
%
%  VENUS('parameter','units','referencebody') returns the corresponding
%  parameter given the spherical and or elliptical reference body
%  specified by 'referencebody'.  If omitted, a sphere is assumed
%  where appropriate.
%
%  See also:  EARTH, JUPITER, MARS, MERCURY, MOON, NEPTUNE, PLUTO,
%             SATURN, SUN, URANUS

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown, W. Stumpf
%   $Revision: 1.6 $    $Date: 1998/08/10 17:50:29 $


%  Define valid parameter strings.  Pre-padding is
%  much faster than strvcat

validparm = [
             'list    '
			 'radius  '
			 'geoid   '
			 'sphere  '
			 'surfarea'
			 'volume  '
		    ];

%  Define valid reference body strings

validref = [
            'sphere'
			'geoid '
			'actual'
		   ];

%  Test the input arguments

if nargin == 0;           parm    = [];     units   = [];     refbody = [];
    elseif nargin == 1;   units   = [];     refbody = [];
    elseif nargin == 2;   refbody = [];
end

%  Initialize output argument
%  Initialization needed for to avoid warning with isempty test at end

mat = [];

%  Empty units tests

if isempty(units);     units   = 'km';       end

%  Test if parameter input, then search for a match

if isempty(parm)
     parm = 'list';
elseif ~isstr(parm)
     error('Input argument PARAMETER must be a string')
else
     strindx = strmatch(lower(parm),validparm);  %  String match
     if length(strindx) == 1
	       parm = deblank(validparm(strindx,:));
	 else
	       error(['Unrecognized parameter:  ',parm])
     end
end

%  Test if parameter input, then search for a match

if isempty(refbody)
     refbody = 'sphere';
elseif ~isstr(refbody)
     error('Input argument REFBODY must be a string')
else
     strindx = strmatch(lower(refbody),validref);  %  String match
     if length(strindx) == 1
	       refbody = deblank(validref(strindx,:));
	 else
	       error(['Unrecognized reference body:  ',refbody])
     end
end

%  Parameters
%    radius, geoid definitions
%    actual volume,  actual surface area
%    spherical volume,  spherical surface area
%    geoid volume, geoid surface area

%  Basic data taken from J.P. Synder, Map Projections A Working
%  Manual, U.S. Geological Survey Paper 1395, US Government Printing
%  Office, Washington, DC, 1987, Table 2, p. 14.  Volume and
%  surface area computed assuming a perfect sphere. Radius agrees
%  with the value tabulated in the Encyclopaedia Britannica. 1995.


radius   = 6051.0;          %  Venus radius in kilometers
volume   = 9.28047e+11;     %  Venus volume in kilometers^3
surfarea = 4.60113e+08;     %  Venus surface area in kilometers^2


if strcmp(parm,'list')
    disp('Function Call:   mat = almanac(''venus'',''parameter'',''units'',''referencebody'')')
    disp(' ')
    disp('Valid parameter strings are:')
    disp('   ''radius''       for the planet radius')
    disp('   ''geoid''        for the planet geoid vector (e = 0 for venus)')
    disp('   ''volume''       for the planet volume')
    disp('   ''surfarea''     for the planet surface area')
    disp(' ')
    disp('Valid units strings are:')
    disp('   ''degrees''       or ''deg'' for degrees')
    disp('   ''kilometers''    or ''km''  for kilometers   (default)')
    disp('   ''nauticalmiles'' or ''nm''  for nautical miles')
    disp('   ''radians''       or ''rad'' for radians')
    disp('   ''statutemiles''  or ''sm''  for statute miles')
    disp(' ')
    disp('Valid reference body strings are:')
    disp('   ''sphere''       for a sphere (default)')
    disp('   ''geoid''        for a ellipsoid (e = 0 for venus)')
    disp('   ''actual''       for the tabulated volume and surface area of the planet')

    return
end

switch refbody
case 'actual'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'venus');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'venus');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'venus');  mat = [mat 0];
	  case 'volume'
        fact = distdim(1.0,'km',units,'venus');
        mat = volume * fact^3;
	  case 'surfarea'
        fact = distdim(1.0,'km',units,'venus');
        mat = surfarea * fact^2;
	end

case 'sphere'
    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'venus');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'venus');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'venus');  mat = [mat 0];
	  case 'volume'
	    rad = distdim(radius,'km',units,'venus');
        mat = sphcalc(rad,parm);
	  case 'surfarea'
	    rad = distdim(radius,'km',units,'venus');
        mat = sphcalc(rad,parm);
	end

case 'geoid'
    warning('Spherical reference body only')

    switch parm
	  case 'radius'
	    mat = distdim(radius,'km',units,'venus');
	  case 'geoid'
	    mat = distdim(radius,'km',units,'venus');  mat = [mat 0];
	  case 'sphere'
	    mat = distdim(radius,'km',units,'venus');  mat = [mat 0];
	  case 'volume'
	    rad = distdim(radius,'km',units,'venus');
        mat = sphcalc(rad,parm);
	  case 'surfarea'
	    rad = distdim(radius,'km',units,'venus');
        mat = sphcalc(rad,parm);
	end

otherwise
    error('Unrecognized reference body string')
end

if isempty(mat)
    error('Unrecognized parameter and reference body combination')
end
