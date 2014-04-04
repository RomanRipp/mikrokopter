function [range,A12,A21] = distbear(lat,lon,arg1,arg2)
% DISTBEAR  Distance and bearing between geographical points.
%
%  [Y,AF,AR]=DISTBEAR(LAT,LON) computes the ranges y between points specified
%  in the lat and lon vectors (decimal degrees with positive indicating
%  north/east). Forward and reverse bearings (degrees) are returned in af, ar.
%
%  [Y,LAT,LON]=DISTBEAR(LAT,LON,N) computes n-point geodesics between successive
%  points. Each successive geodesic occupies it's own row (n >= 2).
%
%  [..] = DISTBEAR(...,'ellipsoid') uses the specified ellipsoid to get distances
%  and bearing. Available ellipsoids are:
%
%  'clarke66'  Clarke 1866
%  'iau73'     IAU 1973
%  'wgs84'     WGS 1984
%  'sphere'    Sphere of radius 6371.0 km
%
%  The default is 'wgs84'.
%
%  Ellipsoid formulas are recommended for distance d < 2000km, but can be
%  used for longer distances.

%  Christian Mertens, IfM Kiel
%  $Revision: 2.0 $ $Date: 1997/08/04 13:45:59 $
%  adapted from `dist.m' by R. Pawlowicz
%
%  renamed to distbear.m to distingiush from simple dist.m
%  by D. Kindler, 97/12/18

%  Notes:
%  RP (WHOI) 3/Dec/91
%    Mostly copied from BDC "dist.f" routine (copied from ....?), but
%    then wildly modified to bring it in line with Matlab vectorization.
%  RP (WHOI) 6/Dec/91
%    Creaturism! - added geodesic computations. This turned
%    out to be pretty hairy since there were a lot of branch problems
%    with asin, atan when computing geodesics subtending > 90 degrees
%    that were ignored in the original code!
%  RP (WHOI) 15/Jan/91
%    Fixed some bothersome special cases, like when computing geodesics
%    and n = 2, or lat = 0...

%  GIVEN THE LATITUDES AND LONGITUDES (IN DEG.) IT ASSUMES THE IAU SPHERO
%  DEFINED IN THE NOTES ON PAGE 523 OF THE EXPLANATORY SUPPLEMENT TO THE
%  AMERICAN EPHEMERIS.
%
%  THIS PROGRAM COMPUTES THE DISTANCE ALONG THE NORMAL
%  SECTION (IN M.) OF A SPECIFIED REFERENCE SPHEROID GIVEN
%  THE GEODETIC LATITUDES AND LONGITUDES OF THE END POINTS
%  *** IN DECIMAL DEGREES ***
%
%  IT USES ROBBIN'S FORMULA, AS GIVEN BY BOMFORD, GEODESY,
%  FOURTH EDITION, P. 122.  CORRECT TO ONE PART IN 10**8
%  AT 1600 KM.  ERRORS OF 20 M AT 5000 KM.
%
%  CHECK:  SMITHSONIAN METEOROLOGICAL TABLES, PP. 483 AND 484,
%  GIVES LENGTHS OF ONE DEGREE OF LATITUDE AND LONGITUDE
%  AS A FUNCTION OF LATITUDE. (SO DOES THE EPHEMERIS ABOVE)
%
%  PETER WORCESTER, AS TOLD TO BRUCE CORNUELLE...1983 MAY 27


lon0 = lon;
lat0 = lat;

spheroid = 'wgs84'; 
geodes = 0;
if (nargin >= 3),
  if (isstr(arg1)),
    spheroid = arg1;
  else
    geodes = 1;
    Ngeodes = arg1;
    if (Ngeodes < 2)
      error('Must have at least 2 points in a goedesic!')
    end
    if (nargin == 4)
      spheroid = arg2;
    end
  end
end

if (spheroid(1:3)=='sph')
  A = 6371000.0;
  B = A;
  E = sqrt(A*A-B*B)/A;
  EPS= E*E/(1 - E*E);
elseif (spheroid(1:3)=='cla')
  A = 6378206.4e0;
  B = 6356583.8e0;
  E= sqrt(A*A - B*B)/A;
  EPS = E*E/(1 - E*E);
elseif (spheroid(1:3)=='iau')
  A = 6378160.e0;
  B = 6356774.516e0;
  E = sqrt(A*A - B*B)/A; 
  EPS = E*E/(1 - E*E);
elseif (spheroid(1:3)=='wgs')
  % on 9/11/88, Peter Worcester gave me (B. Cornuelle) the constants for
  % the WGS84 spheroid, and he gave A (semi-major axis), F = (A-B)/A
  % (flattening) (where B is the semi-minor axis), and E is the
  % eccentricity, E = ((A**2 - B**2)**.5)/A
  % the numbers from Peter are: A = 6378137.; 1/F = 298.257223563
  % E = 0.081819191
  A = 6378137.;
  B = 6356752.314;        %       A * sqrt(1-E*E);
  E = sqrt(A*A - B*B)/A;  %   E = 0.081819191;
  EPS = E*E/(1 - E*E);
else
   error('Unknown spheroid specified!')
end


NN = max(size(lat));
if (NN ~= max(size(lon)))
   error('lat, lon vectors of different sizes!')
end

% It is easier if things are column vectors,
% but we have to fix things before returning!
if (NN == size(lat,1))
  rowvec=0;
else
  rowvec=1;
end 

% convert to radians
lat = lat(:)*pi/180;
lon = lon(:)*pi/180;

% Fixes some nasty 0/0 cases in the geodesics stuff
lat(find(lat==0)) = eps+0*find(lat==0);
                                      
% endpoints of each segment
PHI1 = lat(1:NN-1);    
XLAM1 = lon(1:NN-1);
PHI2 = lat(2:NN);
XLAM2 = lon(2:NN);

% wiggle lines of constant lat to prevent numerical probs.
if (any(PHI1==PHI2))
  for ii = 1:NN-1
    if (PHI1(ii) == PHI2(ii))
      PHI2(ii) = PHI2(ii) + 1e-14;
    end
  end
end

% wiggle lines of constant long to prevent numerical probs.
if (any(XLAM1==XLAM2))
  for ii=1:NN-1
    if (XLAM1(ii)==XLAM2(ii))
      XLAM2(ii)=XLAM2(ii)+ 1e-14;
    end
  end
end


% COMPUTE THE RADIUS OF CURVATURE IN THE PRIME VERTICAL FOR
% EACH POINT
xnu=A./sqrt(1.0-(E*sin(lat)).^2);
xnu1=xnu(1:NN-1);
xnu2=xnu(2:NN);

% COMPUTE THE AZIMUTHS.  A12 (A21) IS THE AZIMUTH AT POINT 1 (2)
% OF THE NORMAL SECTION CONTAINING THE POINT 2 (1)
TPSI2=(1.-E*E)*tan(PHI2) + E*E*xnu1.*sin(PHI1)./(xnu2.*cos(PHI2));
PSI2=atan(TPSI2);

% SOME FORM OF ANGLE DIFFERENCE COMPUTED HERE??
DPHI2=PHI2-PSI2;
DLAM=XLAM2-XLAM1;
CTA12=(cos(PHI1).*TPSI2 - sin(PHI1).*cos(DLAM))./sin(DLAM);
A12=atan((1.)./CTA12);
CTA21P=(sin(PSI2).*cos(DLAM) - cos(PSI2).*tan(PHI1))./sin(DLAM);
A21P=atan((1.)./CTA21P);

%C    GET THE QUADRANT RIGHT
DLAM2=(abs(DLAM)<pi).*DLAM + (DLAM>=pi).*(-2*pi+DLAM) + ...
        (DLAM<=-pi).*(2*pi+DLAM);
A12=A12+(A12<-pi)*2*pi-(A12>=pi)*2*pi;
A12=A12+pi*sign(-A12).*( sign(A12) ~= sign(DLAM2) );
A21P=A21P+(A21P<-pi)*2*pi-(A21P>=pi)*2*pi;
A21P=A21P+pi*sign(-A21P).*( sign(A21P) ~= sign(-DLAM2) );
%%A12*180/pi
%%A21P*180/pi


SSIG=sin(DLAM).*cos(PSI2)./sin(A12);
% At this point we are OK if the angle < 90...but otherwise
% we get the wrong branch of asin! 
% This fudge will correct every case on a sphere, and *almost*
% every case on an ellipsoid (wrong hnadling will be when
% angle is almost exactly 90 degrees)
dd2=[cos(lon).*cos(lat) sin(lon).*cos(lat) sin(lat)];
dd2=sum((diff(dd2).*diff(dd2))')';
if ( any(abs(dd2-2) < 2*((B-A)/A))^2 ),
   disp('dist: Warning...point(s) too close to 90 degrees apart');
end;
bigbrnch=dd2>2;
 
SIG=asin(SSIG).*(bigbrnch==0) + (pi-asin(SSIG)).*bigbrnch;
     
SSIGC=-sin(DLAM).*cos(PHI1)./sin(A21P);
SIGC=asin(SSIGC);
A21 = A21P - DPHI2.*sin(A21P).*tan(SIG/2.0);

%C   COMPUTE RANGE

G2=EPS*(sin(PHI1)).^2;
G=sqrt(G2);
H2=EPS*(cos(PHI1).*cos(A12)).^2;
H=sqrt(H2);
TERM1=-SIG.*SIG.*H2.*(1.0-H2)/6.0;
TERM2=(SIG.^3).*G.*H.*(1.0-2.0*H2)/8.0;
TERM3=(SIG.^4).*(H2.*(4.0-7.0*H2)-3.0*G2.*(1.0-7.0*H2))/120.0;
TERM4=-(SIG.^5).*G.*H/48.0;

range=xnu1.*SIG.*(1.0+TERM1+TERM2+TERM3+TERM4);



try

  range0 = range;

  range = sodanocean(lon0(1),lat0(1),lon0(2),lat0(2));

  nl = char(10);

  fprintf([ nl  'DISTBEAR.M: use Range of SODANO,  SODANO - DISTBEAR = '        ...
                 sprintf('%.0f m',range-range0)  nl   ]);

end


 

if (geodes),

%c now calculate the locations along the ray path. (for extra accuracy, could
%c do it from start to halfway, then from end for the rest, switching from A12
%c to A21...
%c started to use Rudoe's formula, page 117 in Bomford...(1980, fourth edition)
%c but then went to Clarke's best formula (pg 118)

%RP I am doing this twice because this formula doesn't work when we go
%past 90 degrees! 
    Ngd1=round(Ngeodes/2);

% First time...away from point 1
    if (Ngd1>1),
      wns=ones(1,Ngd1);
      CP1CA12 = (cos(PHI1).*cos(A12)).^2;
      R2PRM = -EPS.*CP1CA12;
      R3PRM = 3.0*EPS.*(1.0-R2PRM).*cos(PHI1).*sin(PHI1).*cos(A12);
      C1 = R2PRM.*(1.0+R2PRM)/6.0*wns;
      C2 = R3PRM.*(1.0+3.0*R2PRM)/24.0*wns;
      R2PRM=R2PRM*wns;
      R3PRM=R3PRM*wns;

%c  now have to loop over positions
      RLRAT = (range./xnu1)*([0:Ngd1-1]/(Ngeodes-1));

      THETA = RLRAT.*(1 - (RLRAT.^2).*(C1 - C2.*RLRAT));
      C3 = 1.0 - (R2PRM.*(THETA.^2))/2.0 - (R3PRM.*(THETA.^3))/6.0;
      DSINPSI =(sin(PHI1)*wns).*cos(THETA) + ...
               ((cos(PHI1).*cos(A12))*wns).*sin(THETA);
%try to identify the branch...got to other branch if range> 1/4 circle
       PSI = asin(DSINPSI);

      DCOSPSI = cos(PSI);
      DSINDLA = (sin(A12)*wns).*sin(THETA)./DCOSPSI;
      DTANPHI=(1.0+EPS)*(1.0 - (E^2)*C3.*(sin(PHI1)*wns)./DSINPSI).*tan(PSI);
%C compute output latitude (phi) and long (xla) in radians
%c I believe these are absolute, and don't need source coords added
      PHI = atan(DTANPHI);
%  fix branch cut stuff - 
      otherbrcnh= sign(DLAM2*wns) ~= sign([sign(DLAM2) diff(DSINDLA')'] );
      XLA = XLAM1*wns + asin(DSINDLA).*(otherbrcnh==0) + ...
                       (pi-asin(DSINDLA)).*(otherbrcnh);
    else
      PHI=PHI1;
      XLA=XLAM1;
    end;
      
% Now we do the same thing, but in the reverse direction from the receiver!
    if (Ngeodes-Ngd1>1),
      wns=ones(1,Ngeodes-Ngd1);
      CP2CA21 = (cos(PHI2).*cos(A21)).^2;
      R2PRM = -EPS.*CP2CA21;
      R3PRM = 3.0*EPS.*(1.0-R2PRM).*cos(PHI2).*sin(PHI2).*cos(A21);
      C1 = R2PRM.*(1.0+R2PRM)/6.0*wns;
      C2 = R3PRM.*(1.0+3.0*R2PRM)/24.0*wns;
      R2PRM=R2PRM*wns;
      R3PRM=R3PRM*wns;

%c  now have to loop over positions
      RLRAT = (range./xnu2)*([0:Ngeodes-Ngd1-1]/(Ngeodes-1));

      THETA = RLRAT.*(1 - (RLRAT.^2).*(C1 - C2.*RLRAT));
      C3 = 1.0 - (R2PRM.*(THETA.^2))/2.0 - (R3PRM.*(THETA.^3))/6.0;
      DSINPSI =(sin(PHI2)*wns).*cos(THETA) + ...
               ((cos(PHI2).*cos(A21))*wns).*sin(THETA);
%try to identify the branch...got to other branch if range> 1/4 circle
      PSI = asin(DSINPSI);

      DCOSPSI = cos(PSI);
      DSINDLA = (sin(A21)*wns).*sin(THETA)./DCOSPSI;
      DTANPHI=(1.0+EPS)*(1.0 - (E^2)*C3.*(sin(PHI2)*wns)./DSINPSI).*tan(PSI);
%C compute output latitude (phi) and long (xla) in radians
%c I believe these are absolute, and don't need source coords added
      PHI = [PHI fliplr(atan(DTANPHI))];
% fix branch cut stuff
      otherbrcnh= sign(-DLAM2*wns) ~= sign( [sign(-DLAM2) diff(DSINDLA')'] );
      XLA = [XLA fliplr(XLAM2*wns + asin(DSINDLA).*(otherbrcnh==0) + ...
                       (pi-asin(DSINDLA)).*(otherbrcnh))];
    else
      PHI = [PHI PHI2];
      XLA = [XLA XLAM2];
    end;

%c convert to degrees
      A12 = PHI*180/pi;
      A21 = XLA*180/pi;
      range=range*([0:Ngeodes-1]/(Ngeodes-1));


else

%C*** CONVERT TO DECIMAL DEGREES
   A12=A12*180/pi;
   A21=A21*180/pi;
   if (rowvec),
      range=range';
      A12=A12';
      A21=A21';
   end;
end;




