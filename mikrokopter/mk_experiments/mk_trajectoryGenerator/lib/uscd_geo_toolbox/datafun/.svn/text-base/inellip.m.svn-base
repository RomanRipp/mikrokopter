function [y,x] = inellip(c,r,x,y,mode);

% INELLIP   True for points inside elliptical region
%
% [ Ok , Distance ] = INELLIP( Center , Radius , X , Y )
%
%------------------------------------------------------------
%
% Elliptical region defined by Center and Radius:
%
%  Center = [ CenterX  CenterY ]
%
%  Radius =   Radius                          (Circle)
%  Radius = [ RadiusX  RadiusY ]              (Elliptic)
%  Radius = [ RadiusX  RadiusY   Rotation ]   (rotated)
%
%  Rotation = Angle of RadiusX to X-Achses, default: 0
%
% Note: Ellipse in PolarCoordinates: rad(phi)
%
%  rad(phi) = rx * ry / sqrt( ry² + (rx²-ry²) * sin²(phi-rot) )
% 
%  x = rad(phi) * cos(phi)
%  y = rad(phi) * sin(phi)
%
%------------------------------------------------------------
%
% Ok(p,q) = 1.0 if the point (X(p,q), Y(p,q)) is strictly inside, 
%           0.5 if the point (X(p,q), Y(p,q)) on, 
%           0.0 if the point (X(p,q), Y(p,q)) is outside
%
% the elliptical region specified by Radius arround Center
%
% Distance = Distances of points (X,Y) to Center
%
%------------------------------------------------------------
%
% INELLIP( ... , Mode )  use Geographical Coordinates
%
%  Mode = {'sodano'} | 'robbins'  
%         for Calculation of Distance using GEODIST
%
%  X == Longitude
%  Y == Latitude
%
%  Radius:   [deg]  
%  Rotation == Angle of RadiusX to North, default: 90
%
%  Distance: [deg]
%
%------------------------------------------------------------
%
% see also: INPOLYGON, ONLINE, GEODIST, DST_PROJ
%

if nargin < 4
   error('Not enough InputArguments.');
end

%*******************************************************

if ~( prod(size(c)) == 2 )
    error('Center must be [Lon Lat]');
end

p = prod(size(r));

if ~( ( 0 < p ) & ( p <= 3 ) )
    error('Radius have 1 .. 3 Elements');
end

w = [];

if     p == 1
   r = r([1 1]);
elseif p == 3
   w = r(3);
   r = r(1:2);
end

if nargin == 5
   md = { 's'  'r' };
   ok = chkstr(mode,1);
   if ok
      ok = ( lower(mode(1)) == md{2} );
   end
   mode = md{1+ok};
else
   mode = '';
end

is_mode = ~isempty(mode);

if isempty(w)
   w = 0 + 90 * is_mode;
end

%--------------------------------------------------------

sx = size(x);  px = prod(sx);  vx = ( px == max(sx) );
sy = size(y);  py = prod(sy);  vy = ( py == max(sy) );

if ~isequal(sx,sy)
    if     vx & vy
       x = ones(py,1) * x(:)';
       y = y(:) * ones(1,px);
    elseif vx & ( px == sy(2) ) & ( py == sy(1)*sy(2) )
       x = ones(sy(1),1) * x(:)';
    elseif vy & ( py == sx(1) ) & ( px == sx(1)*sx(2) )
       y = y(:) * ones(1,sx(2));
    else
       error('Size of X and Y must be agree.');
    end
    sx = size(x);
    sy = size(y);
end

if ( px == 0 ) | ( py == 0 )
   return
end

%*******************************************************

r = abs(r);

x = x(:); 
y = y(:);

c = c(:)';

if is_mode
   [x,y] = geodist(c([2 1]),[y x],mode);  % [Distance,Bearing]
   x = x / 1852 / 60;                     % [m] --> [deg]
   y = ( w - y ) * pi/180;                % Bearing to Direction W
else
   x = x - c(1);
   y = y - c(2);
   d = sqrt( (x).^2 + (y).^2 );
   y = atan2(y,x) - w*pi/180;
   x = d;
end

y = sin(y) .^ 2;

y = ( r(2)^2 + ( r(1)^2 - r(2)^2 ) * y );

y = r(1)*r(2) ./ sqrt(y);

y = ( x < y ) + 0.5 * ( x == y );

%--------------------------------------------
% OutPut

x = reshape(x,sx);
y = reshape(y,sx);

%--------------------------------------------


%*******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [DIST,A12,A21] = geodist(p1,p2,varargin)

% GEODIST  Distance and bearing between geographical points.
%
%  [ Distance , AF , AR ] = GEODIST( P1 , P2 ) 
%
%  computes the ranges Distance [m], between points P1 --> P2, 
%   using the WGS84-ellipsoid,
%  AF and AR are the forward and reverse bearings [deg].
% 
%  The Points P1, P2 are specified in [ Lat lon ], i.e. 2 Columns or 2 Rows.
%  P1 and P2 must be either Matrices of same Size or define single Points.
%
%
%  [ D , LAT , LON ] = GEODIST( P1 , P2 , N ) 
%
%  computes N-point geodesics between P1 and P2, with N >= 2
%
%
%  [ D , ... ] = GEODIST( P1 , P2 , ... , Method )
%
%  specify the Method to use for range: {'sodano'}  |  'robbins'
%
%
% 

Nin = nargin;

N   = 0;
fcn = { 'sodano'  'robbins' };

%****************************************************************
% Check Inputs

if Nin < 2
   error('Not enough InputArguments.')
end

[msg,p1,p2,fcn,N,flip] = checkin(p1,p2,N,fcn,varargin{:});

if ~isempty(msg)
   error(msg);
end

%****************************************************************
% convert to Radiant

p1(:,2) = p1(:,2) - 360 * floor( (p1(:,2)+180) / 360 );  % [ -180 .. 180 )
p2(:,2) = p2(:,2) - 360 * floor( (p2(:,2)+180) / 360 );  % [ -180 .. 180 )

p180 = pi/180;

r1 = p1 * p180;
r2 = p2 * p180;

%****************************************************************
% Distance & Bearings

[DIST,A12,A21] = feval(fcn,r1(:,1),r1(:,2),r2(:,1),r2(:,2));

%****************************************************************
% Check for equal Points

eq = ( sum( ( p1 == p2 ) , 2 ) == 2 );

DIST = DIST .* (~eq);

%****************************************************************
% Check for equal Longitude

eq = ( p1(:,2) == p2(:,2) );

A12 = A12 .* (~eq) + pi * ( p1(:,1) >  p2(:,1) ) .* eq;
A21 = A21 .* (~eq) + pi * ( p1(:,1) <= p2(:,1) ) .* eq;

%****************************************************************
if ~( N == 0 )

  %  now calculate the locations along the ray path. (for extra accuracy, could
  %  do it from start to halfway, then from end for the rest, switching from A12
  %  to A21...

  %  RP I am doing this twice because this formula doesn't work when we go
  %  past 90 degrees! 

  N2 = floor( (N+1) / 2 );

  DLAM = quad( r2(:,2) - r1(:,2) );

  [ PHI1 , LAM1 ] = geoline(DIST,r1(:,1),r1(:,2),DLAM,A12,N,  N2, 1);
  [ PHI2 , LAM2 ] = geoline(DIST,r2(:,1),r2(:,2),DLAM,A21,N,N-N2,-1);

  A12 = cat( 2 , PHI1 , PHI2 ) / p180;
  A21 = cat( 2 , LAM1 , LAM2 ) / p180;

  %---------------------------------------------------------------
  % Equal Distance

 acc = 1/2 * 1e3 * eps;

  ok = ( DIST > acc );

  A12 = A12 .* ok(:,ones(1,N)) +  (p1(:,1).*(~ok)) * ones(1,N);
  A21 = A21 .* ok(:,ones(1,N)) +  (p1(:,2).*(~ok)) * ones(1,N);

  %---------------------------------------------------------------
  % Equal Longitude

  ok = ( abs( p1(:,2) - p2(:,2) ) > acc );

  A21 = A21 .* ok(:,ones(1,N)) +  (p1(:,2).*(~ok)) * ones(1,N);
  
  %---------------------------------------------------------------

  A12(:,1) = p1(:,1);
  A12(:,N) = p2(:,1);
 
  A21(:,1) = p1(:,2);
  A21(:,N) = p2(:,2);
 
  DIST = DIST * ( ( 0 : (N-1) ) / (N-1) );

%****************************************************************
else

  % convert to Degree

   A12 = A12 / p180;
   A21 = A21 / p180;

   A12 = A12 .* (~eq) + 180 * ( p1(:,1) >  p2(:,1) ) .* eq;
   A21 = A21 .* (~eq) + 180 * ( p1(:,1) <= p2(:,1) ) .* eq;

end

%****************************************************************
% Flip back

if flip
   DIST = permute(DIST,[2 1]);
   A12  = permute(A12,[2 1]);
   A21  = permute(A21,[2 1]);
end

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [A,B,E,EPS] = wgs84

% [A,B,E,EPS] = wgs84

   A = 6378137.;
   B = 6356752.314;        %       A * sqrt(1-E*E);
   E = sqrt(A*A - B*B)/A;  %   E = 0.081819191;
 EPS = E*E/(1 - E*E);


%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [PHI,LAM] = geoline(DIST,PHI,LAM,DLAM,ANG,N,N2,mode)

% [PHI,LAM] = geoline(DIST,PHI,LAM,ANG,N,N2,mode)
%
%  started to use Rudoe's formula, page 117 in Bomford...(1980, fourth edition)
%  but then went to Clarke's best formula (pg 118)
%

if N2 == 1
   return
end

%****************************************************************
% WGS84

[A,B,E,EPS] = wgs84;

%****************************************************************
% First time...away from point 1

xnu = A ./ sqrt( 1 - (E*sin(PHI)).^2 );

wns = ones(1,N2);

CP1CA = ( cos(PHI) .* cos(ANG)) .^ 2;
R2PRM =  -EPS * CP1CA;
R3PRM = 3*EPS * (1-R2PRM) .* cos(PHI) .* sin(PHI) .* cos(ANG);

C1 = ( R2PRM .* ( 1 +   R2PRM )/06 ) * wns;
C2 = ( R3PRM .* ( 1 + 3*R2PRM )/24 ) * wns;
   
R2PRM = R2PRM * wns;
R3PRM = R3PRM * wns;

%   now have to loop over positions
   
RLRAT = ( DIST ./ xnu ) * ( ( 0 : N2-1 ) / (N-1) );

THETA = RLRAT .* ( 1 - (RLRAT.^2) .* (C1 - C2.*RLRAT) );

C3 = 1 - (R2PRM.*(THETA.^2))/2 - (R3PRM.*(THETA.^3))/6;
  
DSINPSI = (sin(PHI)*wns) .* cos(THETA) + ( (cos(PHI).*cos(ANG)) * wns ).*sin(THETA);

%try to identify the branch...got to other branch if range > 1/4 circle

    PSI = asin(DSINPSI);

DCOSPSI = cos(PSI);
DSINDLA = (sin(ANG)*wns) .* sin(THETA) ./ DCOSPSI;

DTANPHI = (1+EPS) * ( 1 - (E^2) * C3.*(sin(PHI)*wns) ./ DSINPSI ) .* tan(PSI);

%  compute output latitude (phi) and long (xla) in radians
%  I believe these are absolute, and don't need source coords added

PHI = atan(DTANPHI);

%  fix branch cut stuff - 
   
obr = ~( sign(mode*DLAM*wns) == sign([ sign(mode*DLAM) diff(DSINDLA,1,2) ]) );

LAM = LAM*wns + asin(DSINDLA).*(~obr) + (pi-asin(DSINDLA)).*obr;

%****************************************************************
% ( mode == -1 ) ==> FLIP

if ( mode == -1 ) 

   ind = ( N2 : -1 : 1 );

   PHI = PHI(:,ind);
   LAM = LAM(:,ind);

end

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [DIST,A12,A21] = robbins(PHI1,LAM1,PHI2,LAM2)

%  [DIST,ang1,ang2] = robbins( lat1 , lon1 , lat2 , lon2 )
%
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
%
%  Notes:
%  RP (WHOI) 3/Dec/91
%    Mostly copied from BDC "dist.f" routine (copied from ....?), but
%    then wildly modified to bring it in line with Matlab vectorization.
%  RP (WHOI) 6/Dec/91
%    Creaturism! - added geodesic computations. This turned
%    out to be pretty hairy since there were a lot of branch problems
%    with asin, atan when computing geodesics subtending > 90 degrees
%    that were ignored in the original code!
%
%  Christian Mertens, IfM Kiel
%  $Revision: 2.0 $ $Date: 1997/08/04 13:45:59 $
%  adapted from `dist.m' by R. Pawlowicz
%

%****************************************************************
% WGS84

[A,B,E,EPS] = wgs84;

 acc = 1/2 * 1e3 * eps;

%****************************************************************

% COMPUTE THE RADIUS OF CURVATURE IN THE PRIME VERTICAL FOR
% EACH POINT

xnu1 = A ./ sqrt( 1 - (E*sin(PHI1)).^2 );
xnu2 = A ./ sqrt( 1 - (E*sin(PHI2)).^2 );

% We are OK if the angle < 90...but otherwise
% we get the wrong branch of asin! 
% This fudge will correct every case on a sphere, and *almost*
% every case on an ellipsoid (wrong hnadling will be when
% angle is almost exactly 90 degrees)

dd2 =   cat( 2 , cos(LAM1).*cos(PHI1) , sin(LAM1).*cos(PHI1) , sin(PHI1) ) ...
      - cat( 2 , cos(LAM2).*cos(PHI2) , sin(LAM2).*cos(PHI2) , sin(PHI2) );

dd2 = sum( dd2.^2 , 2 );

if any( abs(dd2-2) < 2*((B-A)/A) )^2
   warning('point(s) too close to 90 degrees apart');
end

%****************************************************************
% Fixes some nasty 0/0 cases in the geodesics stuff

PHI1 = sign(PHI1) .* max(abs(PHI1),acc);
PHI1 = PHI1 + acc * ( PHI1 == 0 );

PHI2 = sign(PHI2) .* max(abs(PHI2),acc);
PHI2 = PHI2 + acc * ( PHI2 == 0 );
                                      
%****************************************************************
% wiggle lines of constant lat to prevent numerical probs.

PHI2 = PHI2 - 2*acc * sign(PHI2) .* ( PHI1 == PHI2 );
LAM2 = LAM2 - 2*acc * sign(LAM2) .* ( LAM1 == LAM2 );

%****************************************************************

% COMPUTE THE AZIMUTHS.  A12 (A21) IS THE AZIMUTH AT POINT 1 (2)
% OF THE NORMAL SECTION CONTAINING THE POINT 2 (1)

TPSI2 = (1-E*E)*tan(PHI2) + E*E * xnu1.*sin(PHI1) ./ (xnu2.*cos(PHI2));
PSI2  = atan(TPSI2);

% SOME FORM OF ANGLE DIFFERENCE COMPUTED HERE??
DPHI = PHI2 - PSI2;
DLAM = LAM2 - LAM1;

CTA12 = ( cos(PHI1).*TPSI2     - sin(PHI1).*cos(DLAM) ) ./ sin(DLAM);
CTA21 = ( sin(PSI2).*cos(DLAM) - cos(PSI2).*tan(PHI1) ) ./ sin(DLAM);

A12   = atan( 1 ./ CTA12 );
A21   = atan( 1 ./ CTA21 );

%     GET THE QUADRANT RIGHT

[DLAM2,A12,A21] = quad(DLAM,A12,A21);

SSIG = sin(DLAM) .* cos(PSI2) ./ sin(A12);

% At this point we are OK if the angle < 90...but otherwise
% we get the wrong branch of asin! 
% This fudge will correct every case on a sphere, and *almost*
% every case on an ellipsoid (wrong hnadling will be when
% angle is almost exactly 90 degrees)

bigbrnch = ( dd2 > 2 );
 
SIG   = asin(SSIG).*(bigbrnch==0) + (pi-asin(SSIG)).*bigbrnch;
     
SSIGC = -sin(DLAM).*cos(PHI1)./sin(A21);
SIGC  = asin(SSIGC);

A21 = A21 - DPHI.*sin(A21).*tan(SIG/2.0);

%    COMPUTE RANGE

G2 = EPS * sin(PHI1).^2;
G  = sqrt(G2);

H2 = EPS * ( cos(PHI1) .* cos(A12) ).^2;
H  = sqrt(H2);

TERM1 = -(SIG.^2) .* H2 .*      ( 1 -   H2 ) / 6;
TERM2 =  (SIG.^3) .* G  .* H .* ( 1 - 2*H2 ) / 8;
TERM3 =  (SIG.^4) .*  (   H2 .* ( 4 - 7*H2 ) - ...
                        3*G2 .* ( 1 - 7*H2 ) ) / 120;
TERM4 = -(SIG.^5) .* G .* H / 48;

DIST = xnu1 .* SIG .* (1.0+TERM1+TERM2+TERM3+TERM4);

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [DLAM2,A12,A21] = quad(DLAM,A12,A21);

% Get Quadrant right

Nin  = nargin;
Nout = nargout;

DLAM2 = ( abs(DLAM) <   pi ) .*   DLAM + ...
        (     DLAM  >=  pi ) .* ( DLAM - 2*pi ) + ...
        (     DLAM  <= -pi ) .* ( DLAM + 2*pi );


if ( Nin == 1 ) | ( Nout == 1 )
   A12 = 0;
   A21 = 0;
   return
end

A12 = A12 + 2*pi * ( ( A12 < -pi ) - ( A12 >= pi ) );
A21 = A21 + 2*pi * ( ( A21 < -pi ) - ( A21 >= pi ) );

A12 = A12 +   pi * sign(-A12) .* ( sign(A12) ~= sign( DLAM2) );
A21 = A21 +   pi * sign(-A21) .* ( sign(A21) ~= sign(-DLAM2) );


%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [dist,ang1,ang2]=sodano(lat1,lon1,lat2,lon2)

% function [dist,ang1,ang2]=sodano(lat1,lon1,lat2,lon2)
%
% calculates the distance between two positions in m
% for an ellipsoidal body
%
% second position may be array of positions

% translated from FORTRAN, IFREMER-Brest / EPSHOM-Brest
% Gerd Krahmann, IfM Kiel, Jun 1993
% last change 14. Jun 1994 G. Krahmann


%
%  values of half large and short axis of the ellipsoid WGS72: dgax72, dpax72
%  new values WGS84: ea84, eb84
%

%****************************************************************
% WGS84

[a,b] = wgs84;

%-----------------------------------------------------------------------

Nout = nargout;

%=======================================================================
%
% sodan8 (lat1,lon1,lat2,lon2,dist,a,b,azab,azba)
%
%=======================================================================
%
%  Version CDC.1.01            28 mars 1989     source EPSHOM
%  -------
%
%  Modifications
%  -------------
%  Correction calcul des azimuths                            T.Terre
%  Correction latitudes identiques              ../01/88     F.Jaulgey
%
%  Objet :  Calcul de la distance geodesique (formule de sodano) et des
%  -----    azimuths.
%
%  Entree:    lat1, lon1 : latitude et longitude du point 1 en radians
%  ------     lat2, lon2 : idem pour le point 2
%             a, b         : demi-grand axe et demi petit axe de l'ellip-
%                            soide de reference
%=======================================================================
%

flat  = 1 - b/a;
flat2 = flat * flat;

f1 = flat2 * 1.25;
f2 = flat2 * 0.50;
f3 = flat2 * 0.25;
f4 = flat2 * 0.125;
f5 = flat2 * 0.0625;
f6 = flat2 + flat;

f7 = f6 + 1.0;
f8 = f6 * 0.5;

beta1  = atan( (1-flat) * (sin(lat1)./cos(lat1)) );
sbeta1 = sin(beta1);
cbeta1 = cos(beta1);

beta2  = atan( (1-flat) * (sin(lat2)./cos(lat2)) );
sbeta2 = sin(beta2);
cbeta2 = cos(beta2);

dlat   = lat1 - lat2;
dlon   = lon1 - lon2;

adell  = abs(dlon);

if ( adell >= pi )
   adell = 2*pi - adell;
end

sidel = sin(adell);
codel = cos(adell);

a1 = sbeta1 .* sbeta2;
b1 = cbeta1 .* cbeta2;

cophi = a1 + b1.*codel;

tmp   = ( sbeta2.*cbeta1 - sbeta1.*cbeta2.*codel );
tmp0  = ( sidel.*cbeta2 );

siphi = sqrt( tmp0.*tmp0 + tmp.*tmp );

   i0 = ( siphi == 0 );

siphi = siphi + i0;

  c   = b1 .* ( sidel ./ siphi );
  em  = 1 - c.*c;

  phi = asin(siphi);
  if(cophi<0.)
    phi=pi-phi;
  end

  phisq =   phi .*   phi;
  csphi =     1 ./ siphi;
  ctphi = cophi ./ siphi;

  psyco = siphi .* cophi;

  term1 = f7*phi;
  term2 = a1 .* ( f6*siphi - f2*(phisq.*csphi) );
  term3 = em .* ( f2*(phisq.*ctphi) - f8*(phi+psyco) );
  term4 = f2 * (a1.*a1) .* psyco;
  term5 = (em.*em) .* ( f5*(phi+psyco) - f2*(phisq.*ctphi) - ...
                        f4*(psyco.*(cophi.*cophi)) );
  term6 = f2 * (a1.*em) .* ( phisq.*csphi + psyco.*cophi );

  dist = b*(term1+term2+term3-term4+term5+term6);

  dist = dist .* (1-i0);

if Nout < 2
   return
end

ctaz1 =  ( (sbeta2./cbeta2).*cbeta1 - sbeta1.*codel );
ctaz2 = -( (sbeta1./cbeta1).*cbeta2 - sbeta2.*codel );

ang1 = atan2(sidel,ctaz1);  % azab
ang2 = atan2(sidel,ctaz2);  % azba

ang1 = ang1 .* sign(-dlon);

ang2 = ang2 - pi;
ang2 = ang2 .* sign(-dlon);

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,p1,p2,fcn,N,flip] = checkin(p1,p2,N,fcn,varargin);

% CHECKIN  Check Inputs
%

Nin = nargin - 4;

msg = '';
nl  = char(10);

method = '';

  flip = 0;

%****************************************************************
% Check Coordinates

s1 = size(p1);
s2 = size(p2);

ps1 = prod(s1);
ps2 = prod(s2);

ok = ( ( s1(1)*s1(2) == ps1 )  &  any( s1 == 2 )  &  ...
       ( s2(1)*s2(2) == ps2 )  &  any( s2 == 2 )  &  ...
       ( isequal(s1,s2) | ( ps1 == 2 ) | ( ps2 == 2 ) )  );

if ~ok

   msg = 'Positions must define a single or same Number of Points in [ Lat Lon ].';

else

  f1 = ( ( s1(1) == 2 ) & ~( s1(2) == 2 ) );  % Flip p1
  f2 = ( ( s2(1) == 2 ) & ~( s2(2) == 2 ) );  % Flip p2

  if f1
     p1 = permute(p1,[2 1]);
  end

  if f2
     p2 = permute(p2,[2 1]);
  end

  flip = ( ( f1 & ~( ps1 == 2 ) ) | ( f2 & ~( ps2 == 2 ) ) | ( f1 & f2 ) );

  if  ( ps1 == 2 ) & ~( ps2 == 2 )
     p1 = p1(ones(1,size(p2,1)),:);
  end

  if ~( ps1 == 2 ) & ( ps2 == 2 )
     p2 = p2(ones(1,size(p1,1)),:);
  end

end


%****************************************************************
% Check Npoints and Method

if Nin > 0

  Nin = min( Nin , 2 );

  ok = zeros(1,Nin);

  for ii = 1 : Nin
      v = varargin{ii};
      ok(ii) = chkstr(v,1);
      if ok(ii)
         v = lower(v(1));
         ok(ii) = any( strncmp(fcn,v,1) );
         if ok(ii)
            method = v;
         end
      else
         ok(ii) = isnumeric(v) & ( prod(size(v)) == 1 );
         if ok(ii)
            ok(ii) = 2 * ( isfinite(v) & ( mod(v,1) == 0 ) & ( v >= 2 ) );
            if ok(ii)
               N = v;
            end
         end
      end
  end

  if ( any( ok == 0 ) | ( sum(ok==1) > 1 ) | ( sum(ok==2) > 1 ) )

     msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                'Following Inputs must be any of PointNumber and Method.' , nl , ...
                'PointNumber must be an Integer larger 1, Method must be a String.'  );

  end

end

if isempty(method)  

   fcn = fcn{1};

else

   fcn = fcn{ find(strncmp(fcn,method,1)) };

end


%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = chkstr(str,opt)


% CHKSTR  Checks Input for String
%
%  ok = chkstr(str,Option)
%
%  Option ~= 0 ==> only true for nonempty Strings
%
%   default:   Option == 0  (true for empty Strings)
%

 
if nargin < 2
   opt = 0;
end

ok = ( strcmp( class(str) , 'char' )      & ...
       ( prod(size(str)) == size(str,2) ) & ...
       ( isequal(opt,0) | ~isempty(str) )         );

