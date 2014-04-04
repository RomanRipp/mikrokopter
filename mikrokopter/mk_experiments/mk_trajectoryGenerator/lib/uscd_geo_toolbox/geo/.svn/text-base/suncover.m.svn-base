function x = suncover(c,x,y,mode)

% SUNCOVER  Covering of EarthCoordinates by Sun
%
% C = SUNCOVER( Zenit , Lon , Lat , [Mode] )
% C = SUNCOVER( Date  , Lon , Lat , [Mode] )
% C = SUNCOVER( Day   , Lon , Lat , [Mode] )
%
% C = SUNCOVER(  []   , Lon , Lat , [Mode] )
%
% Zenit = [ Lon  Lat ]    Position of Zenit of Sun
%
% Date  = [ YY MM DD hh mm ss ]
%
% Day   =  DATENUM( Date)
%
% C = SUNCOVER(  []   , Lon , Lat , [Mode] )
% 
%  use actual Date by CLOCK
%
%
% Mode  = Overlap + i*Scale
%
% Overlap  [degree]  smooth Gradient between Night and Day
%
% Scale = [ -2 .. 2 ]
%
% Scale = -2    ZERO
%         -1    1-COS
%          0    SIN^2  |  1-COS^2
%          1    SIN
%          2    ONE
%
% default: Mode = 5 + 0*i  ( 5deg == 20min )
%
%----------------------------------------------------------------------------
%
% see also: SUNZENIT, SUNCYCLE, SPH_PROJ
%
%----------------------------------------------------------------------------
% Code adapted from: AIR_SEA TOOLBOX (version 2.0: 8/9/99)
%                    Rich Pawlowicz 
%
% It is put together from expressions taken from Appendix E in the
% 1978 edition of Almanac for Computers, Nautical Almanac Office, U.S.
% Naval Observatory. They are reduced accuracy expressions valid for the
% years 1800-2100. Solar declination computed from these expressions is
% accurate to at least 1'.
%


Nin = nargin;

if Nin < 3
   error('Not enough InputArguments.');
end

%***************************************************************
% Check Zenit

if isempty(c)
   c = clock;
end

c = c(:)';  s2 = size(c,2);

if ~any( s2 == [ 1  2  6 ] )
    error(['First Inputs must be a single DayNumber, ' ...
           '6-Element-DateVector or 2-Element-ZenitPosition.' ]);
end

if ~( s2 == 2 )
   if ( s2 == 1 )
      c = datevec(c);
   end
   [c(1),c(2)] = zenit(c);
   c = c([1 2]);
end

%*******************************************************
% Check Size of LON and LAT

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

%***************************************************************
% Check for Overlap

if Nin < 4
   mode = 5;
elseif isempty(mode)
   mode = 0;
end

%***************************************************************
% Projection

x = proj(c,x,y);

%***************************************************************
% Overlay

int = real(mode);
scl = imag(mode);

if int == 0
   x = ( x >= 0 );
   return
end

int = int * pi/180;

ii = find( abs(x) < asin(int/2) ); %  Half Z-Intervall : ( -i2 .. i2 

y  = x(ii);

y = 2 * pi * ( sin(y) + int/2 ) / (4*int); % Angle ( 0 .. pi/2 )


sgn = sign(scl);
sgn = sgn + ( sgn == 0 );  %  -1    1 
scl = 2 - min(2,abs(scl));

neg = ( sgn == -1 );

y = neg + sgn * sin( pi/2*neg + sgn*y ) .^ scl;

x     = ( x >= 0 );
x(ii) = y;

%***************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function z = proj(c,x,y,z)

% PROJ  Projection of geodesic SphereCoordinates to Cartesian
%
%  Z = PROJ( CENTER , LON , LAT )
%
%  CENTER = [ CenterLongitude  CenterLatitude ]
%
% The Direction of the Cartesian System are:
%
%    X-Axis: East  + Rotation
%    Y-Axis: North + Rotation
%    Z-Axis: Zenit ==  [ 0  0  0 ] --> CENTER
%
% Units of LON, LAT are in degree.
% Units of X, Y and Z are normalized: [ -1 .. 1 ]
%
% Note:  Values for LON will returned in [ -180 .. 180 ],
%                   NaN-Values for LAT at [ -90  90 ]
%


sx = size(x);

x = x(:); 
y = y(:);

%-------------------------------------
% Basis CoordinateSystem

% ex0 = [ 1 ; 0 ; 0 ];
% ey0 = [ 0 ; 1 ; 0 ];
% ez0 = [ 0 ; 0 ; 1 ];

%-------------------------------------
% New CoordinateAxes in Basis System

c = c * pi/180;

sc = sin(c);
cc = cos(c);

%%% ey = [ -sc(1)*sc(2) ;  cc(1)*sc(2) ; cc(2) ];

    ez = [  sc(1)*cc(2) ; -cc(1)*cc(2) ; sc(2) ];

%%% ex = cross(ey,ez);
%%%  T = [ ex ey ez ] * eye(3,3);

%********************************************************
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%  Length > 100000  ==>  Use parts off Length == 100000  

nn = prod(sx);
mm = 1e5;

m = nn + ( mm - nn ) * ( nn > mm );

n = ceil( nn / m );

z = zeros(nn,1);

for ii = 1 : n

 ende = m*ii + ( nn - m*ii ) * ( nn < m*ii );

 ind = ( m*(ii-1)+1 :  ende );

 z(ind) = doproj(x(ind),y(ind),ez);   % Z-Komponent ONLY

end

%--------------------------------------------
% OutPut

z = reshape(z,sx);


%*******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function z = doproj(x,y,ez)

% function z = doproj(x,y,T)
%
% Input: x,y    ...  [ N by 1 ];  % Column
% Output x,y,z  ...  [ N by 1 ];  % Column
%
% T = [ ex ey ez ]
%

p180 = pi/180;

%--------------------------------------------
% Cartesic Coordinates of Grid in Basis System

z = ~cat( 2 , ( round((x+00)/180) == (x+00)/180 ) , ...
              ( round((x+90)/180) == (x+90)/180 ) , ...
              ( round((y+90)/180) == (y+90)/180 )        );

z(:,1) =  ( z(:,1) & z(:,3) );
z(:,2) =  ( z(:,2) & z(:,3) );
z(:,3) = ~( round((y+00)/180) == (y+00)/180 );

z(:,1) =  sin(x*p180) .* cos(y*p180) .* z(:,1);
z(:,2) = -cos(x*p180) .* cos(y*p180) .* z(:,2);
z(:,3) =  sin(y*p180) .* z(:,3);

%--------------------------------------------
% Projection %  xyz * T == (e0*inv(e1)) * xyz'

z  = ( z * ez );   % Z-Komponent ONLY

is_nan = find(isnan(z));

z = min(max(z,-1),1);

z(is_nan) = NaN;


%***************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [JD,DEC] = zenit(JD);

% SUN_POS   computes Position of SUN
%
% [ ZenitLon  ZenitLat ] = SUN_POS( [ YY MM DD hh mm ss ] )
%
%----------------------------------------------------------------------------
% Code adapted from: AIR_SEA TOOLBOX (version 2.0: 8/9/99)
%                    Rich Pawlowicz 
%----------------------------------------------------------------------------
%
% It is put together from expressions taken from Appendix E in the
% 1978 edition of Almanac for Computers, Nautical Almanac Office, U.S.
% Naval Observatory. They are reduced accuracy expressions valid for the
% years 1800-2100. Solar declination computed from these expressions is
% accurate to at least 1'.
%
% The solar constant (1368.0 W/m^2) represents a mean of satellite measurements
% made over the last sunspot cycle (1979-1995)  taken from 
%  Coffey et al (1995), Earth System Monitor, 6, 6-10.  
%


  d2r = pi/180;   % deg --> rad

% compute Universal Time in hours
   UT = JD(4) + JD(5) / 60 + JD(6) / 3600;

% compute Julian ephemeris date in days (Day 1 is 1 Jan 4713 B.C.=-4712 Jan 1)
  JD = 367 * JD(1) - fix( 7 * ( JD(1) + fix( (JD(2)+9) / 12 ) ) / 4 ) + ...
        fix( 275 * JD(2) / 9 ) + JD(3) + 1721013 + UT/24;

% compute interval in Julian centuries since 1900
  JD = ( JD - 2415020 ) / 36525;

% compute mean anomaly of the sun
   G = 358.475833 + 35999.049750 * JD - 0.000150 * JD.^2;

% compute mean longitude of sun
   L = 279.696678 + 36000.768920 * JD + 0.000303 * JD.^2;

% compute mean anomaly of Jupiter: 225.444651 + 2880 * JD + 154.906654 * JD;
  JP = 225.444651 + 3034.906654 * JD;

% compute mean anomaly of Venus
  VN = 212.603219 + 58517.803875 * JD + 0.001286 * JD.^2;

% compute longitude of the ascending node of the moon's orbit
  NM = 259.183275 - 1934.142008 * JD + 0.002078 * JD.^2;

   G = (  G - 360 * fix(  G / 360 ) ) * d2r;
   L = (  L - 360 * fix(  L / 360 ) ) * d2r;
  JP = ( JP - 360 * fix( JP / 360 ) ) * d2r;
  VN = ( VN - 360 * fix( VN / 360 ) ) * d2r;
  NM = ( NM - 360 * fix( NM / 360 ) + 360 ) * d2r;

% compute sun theta (THETA)
  DEC = +0.397930 * sin(L)       - 0.000040 * cos(L)       ...
        +0.009999 * sin(G-L)     + 0.003334 * sin(G+L)     ...
        +0.000042 * sin(2*G+L)  - 0.000014 * sin(2*G-L)   ...
        -0.000030 * JD.*sin(G-L) - 0.000010 * JD.*sin(G+L) ...
        -0.000208 * JD.*sin(L)   - 0.000039 * sin(NM-L)    ...
        -0.000010 * cos(G-L-JP);

% compute sun rho
  RHO = 1.000421 - 0.033503 * cos(G) - 0.000140 * cos(2*G) + ...
        0.000084 * JD.*cos(G) - 0.000033 * sin(G-JP) + 0.000027 * sin(2*G-2*VN);

% compute declination: DEC = asin( THETA ./ sqrt(RHO) );

   DEC = asin( DEC ./ sqrt(RHO) ) / d2r;

% compute equation of time (in seconds of time)

   JD = 276.697 + (0.98564734*36525) * JD;    % [deg]
   JD = ( JD - 360 * fix( JD / 360 ) ) * d2r;

   JD =   -97.8 * sin(  JD) - 431.3 * cos(  JD) ...
         +596.6 * sin(2*JD) -   1.9 * cos(2*JD) ...
           +4.0 * sin(3*JD) +  19.3 * cos(3*JD) - 12.7 * sin(4*JD);

   JD = JD / 3600;

% compute local hour angle (LHA)
   JD = -15 * ( JD + UT - 12 );
