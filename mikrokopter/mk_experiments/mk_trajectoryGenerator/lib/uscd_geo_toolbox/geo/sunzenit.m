function [JD,DEC,ALT] = sunzenit(JD)

% SUNZENIT  Calculates ZenitPosition of the Sun on Earth
%
% [ Lon , Lat , Alt ] = SUNZENIT( Day )
% [ Lon , Lat , Alt ] = SUNZENIT( [ YY MM DD ] )
% [ Lon , Lat , Alt ] = SUNZENIT( [ YY MM DD  hh mm ss ] )
%
% Day is the decimal Day of DateTime, using DATENUM / DATEVEC.
%
% SUNZENIT without any InputArgument returns the ZenitPosition 
%  at the actual DateTime, using CLOCK.
%
%----------------------------------------------------------------------------
%
% see also: SUNCYCLE, SUNCOVER
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


Nout = nargout;

is_DEC = ( Nout >= 2 );
is_ALT = ( Nout >= 3);

d2r = pi/180;   % deg --> rad

if nargin < 1
   JD = clock;
end

si = size(JD);

is_date = ( ( ndims(JD) == 2 ) & any( si(2) == [3 6] ) );

if isempty(JD)
   if is_date
      si(2) = 1;
      JD= zeros(si);
   end
   DEC = zeros(si);
   ALT = zeros(si);
   return
end

if is_date

   if si(2) == 3
      JD = cat( 2 , JD , zeros(si(1),3) );
   end

   si(2) = 1;

else 

   JD = datevec(JD(:));

end

flip = ~( si(1) == prod(si) );

%*********************************************************************

% compute Universal Time in hours
   UT = JD(:,4) + JD(:,5) / 60 + JD(:,6) / 3600;

% compute Julian ephemeris date in days (Day 1 is 1 Jan 4713 B.C.=-4712 Jan 1)
  JD = 367 * JD(:,1) - fix( 7 * ( JD(:,1) + fix( (JD(:,2)+9) / 12 ) ) / 4 ) + ...
        fix( 275 * JD(:,2) / 9 ) + JD(:,3) + 1721013 + UT/24;

% compute interval in Julian centuries since 1900
  JD = ( JD - 2415020 ) / 36525;

  if flip
     JD = reshape(JD,si);
     UT = reshape(UT,si);
  end

%*********************************************************************
if is_DEC
%*********************************************************************

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
   DEC = DEC ./ sqrt(RHO);

%*********************************************************************
end % DEC
%*********************************************************************

% compute equation of time (in seconds of time)

   JD = 276.697 + (0.98564734*36525) * JD;    % [deg]
   JD = ( JD - 360 * fix( JD / 360 ) ) * d2r;

   JD =   -97.8 * sin(  JD) - 431.3 * cos(  JD) ...
         +596.6 * sin(2*JD) -   1.9 * cos(2*JD) ...
           +4.0 * sin(3*JD) +  19.3 * cos(3*JD) - 12.7 * sin(4*JD);

% compute Greenwich hour angle (LHA)

   JD =  -15 * ( JD/3600 + UT - 12 );

   JD = JD + 360 * ( JD < -180 ) - 360 * ( JD >= 180 );

%*********************************************************************

if ~is_DEC
    return
end

%*********************************************************************
if is_ALT
%*********************************************************************

      ALT = sqrt(1-DEC.^2) .* cos(JD*d2r);
      ALT = asin( ALT ) / d2r;

%*********************************************************************
end
%*********************************************************************

   DEC = asin( DEC ) / d2r;
