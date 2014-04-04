pi = Math.PI;
deg = pi/180.0;
RAD = 180./pi;


function x = sqr(x); x = x.*x; end
function x = int(x); x = fix(x); end
function x = frac(x); x = x - floor(x); end
function x = mod2pi(x); x = mod(x,2*pi); end
function x = round10(x,p); p = 10^p; x = round(x*p)/p; end

function h = HHMM(hh) 
  if ( hh == 0 ), h = ''; return, end
  m = frac(hh)*60;
  h =  int(hh);
  if ( m >= 59.5 ), h = h + 1; m = m - 60; end
  m = round(m);
  h = sprintf(%2.2d:%2.2d = %.0f',h,m,round10(hh,3))
end


function h = HHMMSS(hh) 
  if ( hh == 0 ), h = ''; return, end
  m = frac(hh)*60;
  h =  int(hh);
  s = frac(m)*60;
  m =  int(m);
  if ( s >= 59.5 ), m = m + 1; s = s - 60; end
  if ( m >= 60   ), h = h + 1; m = m - 60; end
  s = round(s);
  h = sprintf(%2.2d:%2.2d:%2.2d = %.0f',h,m,s,round10(hh,4))
end

function sgn = StarSign(lon)
   sgn = { 'Widder'
          'Stier'
          'Zwillinge'
          'Krebs'
          'Löwe'
          'Jungfrau' 
  	  'Waage'
          'Skorpion'
          'Schütze'
          'Steinbock'
          'Wassermann'
          'Fische'      };
  lon = floor( 180/pi * lon / 30 ) + 1;
  lon = lon - 12 * floor(lon/12);
  lon = lon + ( lon == 0 );
  sgn = sgn{lon};
end

function jd = CalcJD(day,mon,year)
% Calculate Julian date: valid only from 1.3.1901 to 28.2.2100
  jd = 2415020.5-64;  %  1.1.1900 - correction of algorithm
  if ( mon <= 2 ),  { year = year - 1; mon = mon + 12; end
  jd = jd + int( (year-1900)*365.25 );
  jd = jd + int( 30.6001*(1+mon) );
  jd = jd + day;
end

function T = GMST(JD)
% Julian Date to Greenwich Mean Sidereal Time
  UT =  frac( JD - 0.5 ) * 24;    % UT in hours
  JD = floor( JD - 0.5 ) + 0.5;   % JD at 0 hours UT
   T = ( JD - 2451545.0 ) / 36525.0;
   T = 6.697374558 + T * ( 2400.051336 + T * 0.000025862 );
   T = T + UT * 1.002737909;
   T = mod( T , 24 );
end

function UT = GMST2UT(JD, gmst)
% Convert Greenweek mean sidereal time to UT
  JD = floor( JD - 0.5 ) + 0.5;  % JD at 0 hours UT
   T = ( JD - 2451545.0 ) / 36525.0;
   T = 6.697374558 + T * ( 2400.051336 + T * 0.000025862 );
   T = mod( T , 24 );
  UT = 0.9972695663*( gmst - T );
end

function lmst = GMST2LMST(gmst, lon)
% Local Mean Sidereal Time, 
% geographical longitude in radians, 
% East is positive
  lon = 180/pi * lon / 15;
  lmst = mod( gmst + lon , 24.);
end

function [rad,dec] = ecl2equ( lon , lat , T )
% Transform ecliptical coordinates (lon/lat) to 
% equatorial coordinates (RA/dec)

  T = ( T - 2451545.0 ) / 36525.; % Epoch 2000 January 1.5
  T = ( 23 + (26+21.45/60)/60 + T*(-46.815 +T*(-0.0006 + T*0.00181) )/3600 );
  T = T * pi / 180;

  ct = cos(T); st = sin(T);
  
  sl = sin(lon);

  rad = atan2( ( sl*ct - tan(lat)*st ) , cos(lon) );
  rad = mod2pi( rad );
  dec = asin( sin(lat)*ct + cos(lat)*st*sl );
  
end


function sun = SunPosition(T)
% Calculate coordinates for Sun
% Coordinates are accurate to about 10s (right ascension) 
% and a few minutes of arc (declination)

  deg = pi/180;

  D  = T - 2447891.5;
  
  eg = 279.403303 * deg;
  wg = 282.768422 * deg;
  e  =   0.016713;
  a  = 149598500;        % km
  d0 =   0.533128 * deg; % angular diameter of Moon at a distance
  
  anom = 2*pi / 365.242191 * D + eg - wg;    % MeanSunAnomaly

  nu   = anom + 2 * e * sin(anom);    % 360.*DEG/pi
  
  lon  =  mod2pi( nu + wg );
  lat  = 0;
  
  dist = ( 1 - sqr(e) ) / ( 1 + e*cos(nu) ); % distance in astronomical units
  diam = d0 / dist;                          % angular diameter in radians
  dist =  a * dist;                          % distance in km
  prlx = 6378.137 / dist;                    % horizonal parallax

  [rad,dec] = ecl2equ( lon , lat , T );
  

  sun = struct( 'lon'  , { lon  } , ...
                'lat'  , { lat  } , ...
            'anomaly'  , { anom } , ...
            'distance' , { dist } , ...
            'diameter' , { diam } , ...
            'parallax' , { prlx } , ...
            'rad'      , { rad  } , ...
            'dec'      , { dec  } , ...
            'sign'     , { StarSign(lon) } );


end


function moon = MoonPosition( sun , T )
% Calculate data and coordinates for the Moon
% Coordinates are accurate to about 1/5 degree (in ecliptic coordinates)

  deg = pi/180;

  D = deg * ( T - 2447891.5 );
  
  % Mean Moon orbit elements as of 1990.0
  l0  = 318.351648 * deg;
  P0  =  36.340410 * deg;
  N0  = 318.510107 * deg;
  i0  =   5.145396 * deg;
  e   =   0.054900;
  a   = 384401;            % km
  dm0 =   0.5181   * deg;  % angular diameter of Moon at a distance
  pr0 =   0.9507   * deg;  % parallax at distance a
  
  l    =     13.1763966 * D + l0;
  anom = l  - 0.1114041 * D - P0;            % Moon's mean anomaly M
  N    = N0 - 0.0529539 * D;                 % Moon's mean ascending node longitude
  C    = l  - sun.lon;
  Ev   = 1.2739 * deg * sin( 2*C - anom );
  Ae   = 0.1858 * deg * sin(sun.anomaly);
  A3   = 0.3700 * deg * sin(sun.anomaly);

  anom = anom + Ev - Ae - A3;                % corrected Moon anomaly

  Ec   = 6.2886 * deg * sin(anom);           % equation of centre
  A4   = 0.2140 * deg * sin(2*anom);

  l2   = l + Ev + Ec - Ae + A4;              % corrected Moon's longitude

  V    = 0.6583 * deg * sin(2*(l2-sun.lon));

  l3   = l2 + V;                             % true orbital longitude;

  N2   = N - 0.16 * deg * sin(sun.anomaly);
  
  lon = mod2pi( N2 + atan2( sin(l3-N2)* cos(i0), cos(l3-N2) ) );
  lat = asin( sin(l3-N2)* sin(i0) );
  orb = l3;
  
  dist = ( 1 - sqr(e) ) / ( 1 + e*cos(anom+Ec) ); % relative distance to semi mayor axis of lunar oribt
  diam = dm0 / dist;                         % angular diameter in radians
  prlx = pr0 / dist;                         % horizontal parallax in radians
  dist =   a * dist;                         % distance in km

  age = mod2pi( l3 - sun.lon );              % Age of Moon in radians since New Moon (0) - Full Moon (pi)
  phs = 0.5 * ( 1 - cos(age) );              % Moon phase, 0-1


  [ rad , dec ] = ecl2equ( lon , lat , T );

  phases = { 'New Moon' , 'Crescent Moon', 'First Quarter' , 'Zunnehmender Mond', 
  	     'Full Moon' , 'Decrescent Moon', 'Last Quarter', 'Abnehmende Sichel', 'New Moon' };

  main   = 1/29.53 * 2*pi; % show 'Newmoon, 'Quarter' for +/-1 day arond the actual event

  p = mod( age , pi/2 );

  ok = ( ( main <= p ) & ( p <= pi/2-main ) );

  p = 2*round( age / (pi/2) ) + ok + 1;


  moon = struct( 'lon'  , { lon  } , ...
                 'lat'  , { lat  } , ...
              'orblon'  , { orb } , ...
             'anomaly'  , { anom } , ...
             'distance' , { dist } , ...
             'diameter' , { diam } , ...
             'parallax' , { prlx } , ...
             'rad'      , { rad  } , ...
             'dec'      , { dec  } , ...
             'age'      , { age  } , ...
             'phase'    , { phs  } , ...
           'MoonPhases' , { phases{p} } , ...
             'sign'     , { StarSign(lon) } );
  

end


function riseset = GMSTRiseSet( rad, dec, lon, lat)
% returns Greenwich sidereal time (hours) of time of rise and set 
% for rad and dec
% at geographic position lon/lat (all values in radians)

 rd = 180/pi / 15;

 dec = acos(-tan(lat)*tan(dec));        % TagBogen
 rad = rad - lon;

 trans = rd*(rad);
 rise  = rd*(rad-dec) + 24; % calculate GMST of rise of object
 rset  = rd*(rad+dec);     % calculate GMST of set of object

  
 riseset = struct( 'transit' , { mod(trans,24) } , ...
                   'rise'    , { mod(rise ,24) } , ...
                   'set'     , { mod(rset ,24) }    );
 
end


function t = interpGMST(t0, t1, t2, tfac)
% Find GMST of rise/set of object from the two calculates (start)points (day 1 and 2) 
% and at midnight UT(0)
  t2 = ( t2 - t1 );
  t  = ( tfac*24.07*t1 - t0*t2 ) / ( tfac*24.07 - t2 );
end


function RiseSet( jd, c1, c2, lon, lat, tint)
% JD is the Julian Date of 0h UTC time (midnight)

  r1 = GMSTRiseSet( c1.rad , c1.dec , lon , lat);
  r2 = GMSTRiseSet( c2.rad , c2.dec , lon, lat);
    
  %alert( rise1.set  +"  "+ rise2.set );
  % unwrap GMST in case we move across 24h -> 0h
  if ( ( r1.transit > r2.transit ) & ( abs(r1.transit-r2.transit)>18 ) ), r2.transit = r2.transit + 24;
  if ( ( r1.rise    > r2.rise    ) & ( abs(r1.rise   -r2.rise)   >18 ) ), r2.rise    = r2.rise + 24;
  if ( ( r1.set     > r2.set     ) & ( abs(r1.set    -r2.set)    >18 ) ), r2.set     = r2.set + 24;

  T0 = GMST(jd);
%  var T02 = T0-zone*1.002738; % Greenwich sidereal time at 0h time zone (zone: hours)

  % Greenwich sidereal time for 0h at selected longitude
  T1 = T0-lon*180/pi/15*1.002738; if (T1 < 0) T1 = T1 + 24: end

  if (r1.transit < T1), r1.transit = r1.transit + 24; r2.transit = r2.transit + 24; end
  if (r1.rise    < T1), r1.rise    = r1.rise + 24;    r2.rise    = r2.rise + 24; end
  if (r1.set     < T1), r1.set     = r1.set + 24;     r2.set     = r2.set + 24; end
  
  %alert("after="+ rise1.set  +"  "+ rise2.set+ " T0="+ T0 );
 
  % Refraction and Parallax correction
  
  decMean = 0.5 * ( c1.dec + c2.dec );
  psi     = acos( sin(lat) / cos(decMean) );

  % altitude of sun center: semi-diameter, horizontal parallax and (standard) refraction of 34'
  alt     = 0.5 * c1.diameter - c1.parallax + 34/60 * pi/180;
  y       = asin( sin(alt) / sin(psi) );

  dt      = 240*180/pi * y/cos(decMean)/3600;  % time correction due to refraction, parallax
  %alert("T02="+T02+"  "+rise1.rise+"  " +rise1.transit+"  "+rise1.set + "  "+rise2.rise+"  " +rise2.transit+"  "+rise2.set);

  rise = struct( 'transit' , {} , ...
                 'rise'    , {} , ...
                 'set'     , {}        );

  rise.transit = GMST2UT( jd, interpGMST( T0, r1.transit, r2.transit, tint)     );
  rise.rise    = GMST2UT( jd, interpGMST( T0, r1.rise,    r2.rise,    tint) -dt );
  rise.set     = GMST2UT( jd, interpGMST( T0, r1.set,     r2.set,     tint) +dt );
  
  %rise.transit = Mod(rise.transit, 24.);
  %rise.rise    = Mod(rise.rise, 24.);
  %rise.set     = Mod(rise.set,  24.);
 
}


% Find (local) time of sunrise and sunset
% JD is the Julian Date of 0h local time (midnight)
% Accurate to about 1-2 minutes
% recursive: 1 - calculate rise/set in UTC
% recursive: 0 - find rise/set on the current local day (set could also be first)
function SunRise(JD, deltaT, lon, lat, zone, recursive)
{
  var jd0UT = Math.floor(JD-0.5)+0.5;   % JD at 0 hours UT
%  alert("jd0UT="+jd0UT+"  JD="+JD);
  var coor1 = SunPosition(jd0UT+0*deltaT/24/3600);
  var coor2 = SunPosition(jd0UT+1+0*deltaT/24/3600); % calculations for next day's UTC midnight
  
  var risetemp = new Object();
  var rise = new Object();
  rise = RiseSet(jd0UT, coor1, coor2, lon, lat, 1 ); % rise/set time in UTC
  if (!recursive) { % check and adjust to have rise/set time on local calendar day
    if (zone>0) {
	  if (rise.rise>=24-zone || rise.transit>=24-zone || rise.set>=24-zone) {% rise time was yesterday local time -> calculate rise time for next UTC day
		risetemp = SunRise(JD+1, deltaT, lon, lat, zone, 1);
		if (rise.rise>=24-zone) rise.rise = risetemp.rise;
		if (rise.transit >=24-zone) rise.transit = risetemp.transit;
		if (rise.set >=24-zone) rise.set  = risetemp.set;
	  }
	}
    else if (zone<0) {
	  if (rise.rise<-zone || rise.transit<-zone || rise.set<-zone) {% rise time was yesterday local time -> calculate rise time for next UTC day
		risetemp = SunRise(JD-1, deltaT, lon, lat, zone, 1);
		if (rise.rise<-zone) rise.rise = risetemp.rise;
		if (rise.transit<-zone) rise.transit = risetemp.transit;
		if (rise.set <-zone) rise.set  = risetemp.set;
	  }
	}
	rise.transit = Mod(rise.transit+zone, 24);
	rise.rise    = Mod(rise.rise   +zone, 24);
	rise.set     = Mod(rise.set    +zone, 24);
  }
  return( rise );  
}



% Find local time of moonrise and moonset
% JD is the Julian Date of 0h local time (midnight)
% Accurate to about 5 minutes or better
% recursive: 1 - calculate rise/set in UTC
% recursive: 0 - find rise/set on the current local day (set could also be first)
% returns '' for moonrise/set does not occur on selected day
function MoonRise(JD, deltaT, lon, lat, zone, recursive)
{
  var timeinterval = 0.5;
  
  var jd0UT = Math.floor(JD-0.5)+0.5;   % JD at 0 hours UT
  var suncoor1 = SunPosition(jd0UT+0*deltaT/24/3600);
  var coor1 = MoonPosition(suncoor1, jd0UT+0*deltaT/24/3600);

  var suncoor2 = SunPosition(jd0UT +timeinterval +0*deltaT/24/3600); % calculations for noon
  % calculations for next day's midnight
  var coor2 = MoonPosition(suncoor2, jd0UT +timeinterval +0*deltaT/24/3600); 
  
  var risetemp = new Object();
  var rise = new Object();
  
  rise = RiseSet(jd0UT, coor1, coor2, lon, lat, timeinterval); % rise/set time in UTC, time zone corrected later
  
  if (!recursive) { % check and adjust to have rise/set time on local calendar day
    if (zone>0) {
	  riseprev = MoonRise(JD-1, deltaT, lon, lat, zone, 1); % recursive call to MoonRise returns events in UTC
	  
	  %risenext = MoonRise(JD+1, deltaT, lon, lat, zone, 1); % recursive call to MoonRise returns events in UTC
  	  %alert("yesterday="+riseprev.transit+"  today="+rise.transit+" tomorrow="+risenext.transit);
  	  %alert("yesterday="+riseprev.rise+"  today="+rise.rise+" tomorrow="+risenext.rise);
  	  %alert("yesterday="+riseprev.set+"  today="+rise.set+" tomorrow="+risenext.set);

	  if (rise.transit >= 24-zone || rise.transit < -zone) { % transit time is tomorrow local time
		if (riseprev.transit < 24-zone) rise.transit = ''; % there is no moontransit today
		else rise.transit  = riseprev.transit;
	  }

	  if (rise.rise >= 24-zone || rise.rise < -zone) { % transit time is tomorrow local time
		if (riseprev.rise < 24-zone) rise.rise = ''; % there is no moontransit today
		else rise.rise  = riseprev.rise;
	  }

	  if (rise.set >= 24-zone || rise.set < -zone) { % transit time is tomorrow local time
		if (riseprev.set < 24-zone) rise.set = ''; % there is no moontransit today
		else rise.set  = riseprev.set;
	  }

	}
    else if (zone<0) {
	  % rise/set time was tomorrow local time -> calculate rise time for former UTC day
	  if (rise.rise<-zone || rise.set<-zone || rise.transit<-zone) { 
		risetemp = MoonRise(JD+1, deltaT, lon, lat, zone, 1);
		
		if (rise.rise < -zone) {
		  if (risetemp.rise > -zone) rise.rise = ''; % there is no moonrise today
		  else rise.rise = risetemp.rise;
		}
		
		if (rise.transit < -zone)
		{
		  if (risetemp.transit > -zone)  rise.transit = ''; % there is no moonset today
		  else rise.transit  = risetemp.transit;
		}
		
		if (rise.set < -zone)
		{
		  if (risetemp.set > -zone)  rise.set = ''; % there is no moonset today
		  else rise.set  = risetemp.set;
		}
		
	  }
	}
	
	if (rise.rise)    rise.rise = Mod(rise.rise+zone, 24);    % correct for time zone, if time is valid
	if (rise.transit) rise.transit  = Mod(rise.transit +zone, 24); % correct for time zone, if time is valid
	if (rise.set)     rise.set  = Mod(rise.set +zone, 24);    % correct for time zone, if time is valid
  }
  return( rise );  
}


