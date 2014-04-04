function [ssp,info] = snd_spd(D,T,S,equation,lat);

% SND_SPD Various sound-speed equations.
%
%          1)  SND_SPD(D,T,S) returns the sound speed (m/sec) given vectors
%             of salinity (ppt), temperature (deg C) and DEPTH !!!  (m) using
%             the formula of Mackenzie:
%
%           Mackenzie, K.V. "Nine-term Equation for Sound Speed in the Oceans",
%           J. Acoust. Soc. Am. 70 (1981), 807-812.
%
%
%          2) SND_SPD(P,T,S,'del grosso',lat) returns the sound speed (m/sec) 
%            given vectors of salinity (ppt), temperature (deg C), and 
%            PRESSURE (dbar) using the Del Grosso equation:
%
%            Del Grosso, "A New Equation for the speed of sound in Natural
%            Waters", J. Acoust. Soc. Am. 56#4 (1974).
%
%          3) SND_SPD(P,T,S,'chen') returns the sound speed (m/sec) given
%            given vectors of salinity (ppt), temperature (deg C), and 
%            PRESSURE (dbar) using the Chen and Millero equation:
%
%            Chen and Millero, "The Sound Speed in Seawater", J. Acoust. 
%            Soc. Am. 62 (1977), 1129-1135
%
%          4) SND_SPD(P,T,S,'state') returns the sound speed (m/sec) given
%            given vectors of salinity (ppt), temperature (deg C), and 
%            PRESSURE (dbar) by using derivatives of the EOS80 equation of
%            state for seawater and the adiabatic lapse rate.
% 
%  [SSP,Info] = SND_SPD( ... ) returns the sound speed SSP and 
%                               the [ 3 by 1 ] CellStringArray Info with:
%                               { Author ; Title ; Reference }
%
%  requires m-files:  ADIABATT.M,
%                     SWSTATE.M,
%

%	converted to matlab v 5.1  AN


if (nargin<4), equation='mackenzie'; end;
if (nargin<5), lat=45; end;

if ~ischar(equation) | isempty(equation)
 equation='mackenzie';
end

equation = [ upper(equation(1,:))  '   ' ];

switch equation(1:3)

%********************************************************************
case 'MAC'

% Fixed a small typo with the salinity stuff in the T*(S-35) term
%  15/11/91

     c= 1.44896e3;     t= 4.591e0; 
    t2=-5.304e-2;     t3= 2.374e-4;
     s= 1.340e0;       d= 1.630e-2;
    d2= 1.675e-7;     ts=-1.025e-2;
   td3=-7.139e-13;  


   ssp=c+t*T+t2*T.*T+t3*T.*T.*T+s*(S-35.0)+d*D+d2*D.*D+ts*T.*(S-35.0) ...
        +td3*T.*D.*D.*D;

   info = { 'Mackenzie'
            'Nine-term Equation for Sound Speed in the Oceans'
            'J. Acoust. Soc. Am. 70 (1981), 807-812'  };



%********************************************************************
case 'DEL'

% Del grosso uses pressure in kg/cm^2. To get to this from dbars we must 
% divide by "g". From the UNESCO algorithms (referring to ANON (1970) 
% BULLETIN GEODESIQUE) we have this formula for g as a function of latitude 
% and pressure. We set latitude to 45 degrees for convenience!

      XX=sin(lat*pi/180).^2;

      GR = 9.780318*(1.0+(5.2788E-3+2.36E-5*XX).*XX) + 1.092E-6*D;

      P=D./GR;

% This is from VSOUND.f:
    

      C000 = 1402.392;
      DCT = (0.501109398873e1-(0.550946843172e-1 - 0.221535969240e-3*T).*T).*T;
      DCS = (0.132952290781e1 + 0.128955756844e-3*S).*S;
      DCP = (0.156059257041e0 + (0.244998688441e-4 - 0.883392332513e-8*P).*P).*P;
      DCSTP = -0.127562783426e-1*T.*S + 0.635191613389e-2*T.*P +0.265484716608e-7*T.*T.*P.*P ...
 - 0.159349479045e-5*T.*P.*P+0.522116437235e-9*T.*P.*P.*P - 0.438031096213e-6*T.*T.*T.*P;
      DCSTP=DCSTP - 0.161674495909e-8*S.*S.*P.*P + 0.968403156410e-4*T.*T.*S+ ...
   0.485639620015e-5*T.*S.*S.*P - 0.340597039004e-3*T.*S.*P;
   ssp= C000 + DCT + DCS + DCP + DCSTP;


   info = {'Del Grosso'
           'A New Equation for the speed of sound in Natural Waters'
           'J. Acoust. Soc. Am. 56#4 (1974)' };


%********************************************************************
case 'CHE'


      P0=D;
% This is copied directly from the UNESCO algorithmms, with some minor changes (like adding
% ";" and changing "*" to ".*") for Matlab.

% CHECKVALUE: SVEL=1731.995 M/S, S=40 (IPSS-78),T=40 DEG C,P=10000 DBAR

%   SCALE PRESSURE TO BARS
      P=P0/10.;
%**************************
      SR = sqrt(abs(S));
% S**2 TERM
      D = 1.727E-3 - 7.9836E-6*P;
% S**3/2 TERM
      B1 = 7.3637E-5 +1.7945E-7*T;
      B0 = -1.922E-2 -4.42E-5*T;
      B = B0 + B1.*P;
% S**1 TERM
      A3 = (-3.389E-13*T+6.649E-12).*T+1.100E-10;
      A2 = ((7.988E-12*T-1.6002E-10).*T+9.1041E-9).*T-3.9064E-7;
      A1 = (((-2.0122E-10*T+1.0507E-8).*T-6.4885E-8).*T-1.2580E-5).*T+9.4742E-5;
      A0 = (((-3.21E-8*T+2.006E-6).*T+7.164E-5).*T-1.262E-2).*T+1.389;
      A = ((A3.*P+A2).*P+A1).*P+A0;
% S**0 TERM
      C3 = (-2.3643E-12*T+3.8504E-10).*T-9.7729E-9;
      C2 = (((1.0405E-12*T-2.5335E-10).*T+2.5974E-8).*T-1.7107E-6).*T +3.1260E-5;
      C1 = (((-6.1185E-10*T+1.3621E-7).*T-8.1788E-6).*T+6.8982E-4).*T +0.153563;
      C0 = ((((3.1464E-9*T-1.47800E-6).*T+3.3420E-4).*T-5.80852E-2).*T+5.03711).*T+1402.388;
      C = ((C3.*P+C2).*P+C1).*P+C0;
% SOUND SPEED RETURN
      ssp = C + (A+B.*SR+D.*S).*S;


      info = {'Chen and Millero'
              'The Sound Speed in Seawater'
              'J. Acoust. Soc. Am. 62 (1977), 1129-1135'    };



%********************************************************************
case 'STA'

     P=D;
% (Copied somewhat from program EOSSPEED.F)
     [svan,sigma]=swstate(S,T,P);
     VOL=(1.)./(1000.+sigma);
%     DV/DP|ADIA = (DV/DP) AT CONSTANT T + ADIA.LAPSE RATE *
%                  (DV/DT) AT CONSTANT P
%     Note: factor of 10 is convert pressure from dB to Bars
     dVdP=swstate(S,T,P,'dP');
     dVdT=swstate(S,T,P,'dT');
     dVdPad=(dVdP+adiabatt(S,T,P).*dVdT)*10;
%     C = V * SQRT ( 1/DV/DP| ADIA)
     ssp=VOL.*sqrt(abs( (1.e5)./dVdPad ));

     info = {'EOS80'
             'EOS80 equation of state for seawater and the adiabatic lapse rate'
             '' };

%********************************************************************
otherwise

   ssp = [];

end;

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                    
function [SVAN,SIGMA]=swstate(S,T,P0,ftn);
% SWSTATE State equation for seawater
%
%        [SVAN,SIGMA]=SWSTATE(S,T,P) returns the specific volume
%        anomaly SVAN (m^3/kg*1e-8) and the density anomaly SIGMA (kg/m^3)
%        given the salinity S (ppt), temperature T (deg C) and pressure
%        P (dbars).
%
%        [dVdT,dRdT]=SWSTATE(S,T,P,'dT') returns derivatives w.r.t.
%        temperature of the volume and density.
%
%        [dVdS,dRdS]=SWSTATE(S,T,P,'dS') returns derivatives w.r.t.
%        salinity.
%
%        [dVdP,dRdP]=SWSTATE(S,T,P,'dP') returns derivatives w.r.t.
%        pressure.
%
%        All elements can be scalars, vectors, or matrices but should be
%        the same size.

%Notes: RP (WHOI 2/dec/91)
%
%  This stuff is directly copied from the UNESCO algorithms, with some
%  minor changes to make it Matlab compatible (like adding ";" and changing
%  "*" to ".*" when necessary.
%
%      RP  (WHOI 3/dec/91)
%
%  Added first derivative calculations.

derivT=0;
derivS=0;
derivP=0;

if (nargin==4),
   if     (ftn=='dT'), derivT=1;
   elseif (ftn=='dS'), derivS=1;
   elseif (ftn=='dP'), derivP=1;
   else error('swstate: Unrecognized option!');
   end;
end;

% ******************************************************
% SPECIFIC VOLUME ANOMALY (STERIC ANOMALY) BASED ON 1980 EQUATION
% OF STATE FOR SEAWATER AND 1978 PRACTICAL SALINITY SCALE.
% REFERENCES
% MILLERO, ET AL (1980) DEEP-SEA RES.,27A,255-264
% MILLERO AND POISSON 1981,DEEP-SEA RES.,28A PP 625-629.
% BOTH ABOVE REFERENCES ARE ALSO FOUND IN UNESCO REPORT 38 (1981)
% UNITS:      
%       PRESSURE        P0       DECIBARS
%       TEMPERATURE     T        DEG CELSIUS (IPTS-68)
%       SALINITY        S        (IPSS-78)
%       SPEC. VOL. ANA. SVAN     M**3/KG *1.0E-8
%       DENSITY ANA.    SIGMA    KG/M**3
% ******************************************************************
% CHECK VALUE: SVAN=981.3021 E-8 M**3/KG.  FOR S = 40 (IPSS-78) ,
% T = 40 DEG C, P0= 10000 DECIBARS.
% CHECK VALUE: SIGMA = 59.82037  KG/M**3 FOR S = 40 (IPSS-78) ,
% T = 40 DEG C, P0= 10000 DECIBARS.
%HECK VALUE: FOR S = 40 (IPSS-78) , T = 40 DEG C, P0= 10000 DECIBARS.
%        DR/DP                  DR/DT                 DR/DS
%       DRV(1,7)              DRV(2,3)             DRV(1,8)
%
% FINITE DIFFERENCE WITH 3RD ORDER CORRECTION DONE IN DOUBLE PRECSION
%
%       3.46969238E-3       -.43311722           .705110777
%
% EXPLICIT DIFFERENTIATION SINGLE PRECISION FORMULATION EOS80 
% 
%       3.4696929E-3        -.4331173            .7051107
%
% (RP...I think this ---------^^^^^^ should be -.4431173!);


% *******************************************************
% DATA
   R3500=1028.1063;
   R4=4.8314E-4;
   DR350=28.106331;

% CONVERT PRESSURE TO BARS AND TAKE SQUARE ROOT SALINITY.
      P=P0/10.;
      SAL=S;
      SR = sqrt(abs(S));
% *********************************************************
% PURE WATER DENSITY AT ATMOSPHERIC PRESSURE
%   BIGG P.H.,(1967) BR. J. APPLIED PHYSICS 8 PP 521-537.
%
      R1 = ((((6.536332E-9*T-1.120083E-6).*T+1.001685E-4).*T ...
            -9.095290E-3).*T+6.793952E-2).*T-28.263737;
% SEAWATER DENSITY ATM PRESS. 
%  COEFFICIENTS INVOLVING SALINITY
      R2 = (((5.3875E-9*T-8.2467E-7).*T+7.6438E-5).*T-4.0899E-3).*T+8.24493E-1; 
      R3 = (-1.6546E-6*T+1.0227E-4).*T-5.72466E-3;
%  INTERNATIONAL ONE-ATMOSPHERE EQUATION OF STATE OF SEAWATER
      SIG = (R4*S + R3.*SR + R2).*S + R1;
% SPECIFIC VOLUME AT ATMOSPHERIC PRESSURE
      V350P = 1.0/R3500;
      SVA = -SIG*V350P./(R3500+SIG);
      SIGMA=SIG+DR350;
      V0 = 1.0./(1000.0 + SIGMA);
%  SCALE SPECIFIC VOL. ANAMOLY TO NORMALLY REPORTED UNITS
      SVAN=SVA*1.0E+8;

if (derivS),               % These are derivatives for (S,T,0).
      R4S=9.6628E-4;
      RHO1 = 1000.0 + SIGMA;

      RHOS=R4S*SAL+1.5.*R3.*SR+R2;  
      V0S=-RHOS./(RHO1.*RHO1);  
elseif (derivT),
      R1 =(((3.268166E-8*T-4.480332E-6).*T+3.005055E-4).*T...
          -1.819058E-2).*T+6.793952E-2;
      R2 = ((2.155E-8*T-2.47401E-6).*T+1.52876E-4).*T-4.0899E-3;
      R3 = -3.3092E-6*T+1.0227E-4;
      RHO1 = 1000.0 + SIGMA;

      RHOT = (R3.*SR + R2).*SAL + R1;    
      V0T = -RHOT./(RHO1.*RHO1);
end;
      
% ******************************************************************
% ******  NEW HIGH PRESSURE EQUATION OF STATE FOR SEAWATER ********
% ******************************************************************
%        MILLERO, ET AL , 1980 DSR 27A, PP 255-264
%               CONSTANT NOTATION FOLLOWS ARTICLE
%********************************************************
% COMPUTE COMPRESSION TERMS
      E = (9.1697E-10*T+2.0816E-8).*T-9.9348E-7;
      BW = (5.2787E-8*T-6.12293E-6).*T+3.47718E-5;
      B = BW + E.*S;    % Bulk Modulus (almost)
%  CORRECT B FOR ANAMOLY BIAS CHANGE
      Bout = B + 5.03217E-5;

if (derivS),
      DBDS=E;
elseif (derivT),
      BW = 1.05574E-7*T-6.12293E-6;
      E = 1.83394E-9*T +2.0816E-8;
      BT = BW + E.*SAL;
end;
%             
      D = 1.91075E-4;
      C = (-1.6078E-6*T-1.0981E-5).*T+2.2838E-3;
      AW = ((-5.77905E-7*T+1.16092E-4).*T+1.43713E-3).*T-0.1194975;
      A = (D*SR + C).*S + AW;    
%  CORRECT A FOR ANAMOLY BIAS CHANGE
      Aout = A + 3.3594055;

if (derivS),
      DADS=2.866125E-4*SR+C;
elseif (derivT),
      C = -3.2156E-6*T -1.0981E-5;
      AW = (-1.733715E-6*T+2.32184E-4).*T+1.43713E-3;
      AT = C.*SAL + AW;
end;
           
      B1 = (-5.3009E-4*T+1.6483E-2).*T+7.944E-2;
      A1 = ((-6.1670E-5*T+1.09987E-2).*T-0.603459).*T+54.6746;
      KW = (((-5.155288E-5*T+1.360477E-2).*T-2.327105).*T+148.4206).*T-1930.06;
      K0 = (B1.*SR + A1).*S + KW;

if (derivS),
      K0S=1.5*B1.*SR+A1;
      KS=(DBDS.*P+DADS).*P+K0S;
elseif (derivT),
      B1 = -1.06018E-3*T+1.6483E-2;
      % APRIL 9 1984 CORRECT A1 BIAS FROM -.603457 !!!
      A1 = (-1.8501E-4*T+2.19974E-2).*T-0.603459;
      KW = ((-2.0621152E-4*T+4.081431E-2).*T-4.65421).*T+148.4206;
      K0T = (B1.*SR+A1).*SAL + KW;
      KT = (BT.*P + AT).*P + K0T;
end;


% EVALUATE PRESSURE POLYNOMIAL 
% ***********************************************
%   K EQUALS THE SECANT BULK MODULUS OF SEAWATER
%   DK=K(S,T,P)-K(35,0,P)
%  K35=K(35,0,P)
% ***********************************************
      DK = (B.*P + A).*P + K0;
      K35  = (5.03217E-5*P+3.359406).*P+21582.27;
      GAM=P./K35;
      PK = 1.0 - GAM;
      SVA = SVA.*PK + (V350P+SVA).*P.*DK./(K35.*(K35+DK));
%  SCALE SPECIFIC VOL. ANAMOLY TO NORMALLY REPORTED UNITS
      SVAN=SVA*1.0E+8;      % Volume anomaly
      V350P = V350P.*PK;
%  ****************************************************
% COMPUTE DENSITY ANAMOLY WITH RESPECT TO 1000.0 KG/M**3
%  1) DR350: DENSITY ANAMOLY AT 35 (IPSS-78), 0 DEG. C AND 0 DECIBARS
%  2) DR35P: DENSITY ANAMOLY 35 (IPSS-78), 0 DEG. C ,  PRES. VARIATION
%  3) DVAN : DENSITY ANAMOLY VARIATIONS INVOLVING SPECFIC VOL. ANAMOLY
% ********************************************************************
% CHECK VALUE: SIGMA = 59.82037  KG/M**3 FOR S = 40 (IPSS-78),
% T = 40 DEG C, P0= 10000 DECIBARS.
% *******************************************************
      DR35P=GAM./V350P;
      DVAN=SVA./(V350P.*(V350P+SVA));
      SIGMA=DR350+DR35P-DVAN;  % Density anomaly

      K=K35+DK;
      VP=1.0-P./K;
      V = (1.) ./(SIGMA+1000.0);

if (derivS),
      VS=V0S.*VP+V0.*P.*KS./(K.*K);

      SVAN=VS;              % dVdS
      SIGMA=-VS./(V.*V);    % dRdS
elseif (derivT),
      VT = V0T.*VP + V0.*P.*KT./(K.*K);

      SVAN=VT;              % dVdT
      SIGMA=-VT./(V.*V);    % dRdT
elseif (derivP),
      DKDP = 2.0*Bout.*P + Aout;   
% CORRECT DVDP TO PER DECIBAR BY MULTIPLE *.1
      DVDP = -.1*V0.*(1.0 - P.*DKDP./K)./K;

      SVAN=DVDP;            % dVdP
      SIGMA=-DVDP./(V.*V);  % dRdP
end;

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ATG=adiabatt(S,T,P);
% ADIABATT Computes the adiabatic temperature gradient (called
%                 by POTTEMP).
%
%                 ADIABATT(S,T,P) returns the gradient. Units are:
%   
%       PRESSURE        P        DECIBARS
%       TEMPERATURE     T        DEG CELSIUS (IPTS-68)
%       SALINITY        S        (IPSS-78)
%       ADIABATIC       ATG      DEG. C/DECIBAR
%
% REF: BRYDEN,H.,1973,DEEP-SEA RES.,20,401-408

% Notes: RP 29/Nov/91
%
% I have modified the FORTRAN code to make it Matlab compatible, but
% no numbers have been changed. In certain places "*" has been replaced
% with ".*" to allow vectorization.
%
% Modified for PC 8 character length and Matlab version 5.1  (AN)
%
%
% This routine is called by pottemp.m, and is called ATG in the
% UNESCO algorithms.

%C CHECKVALUE: ATG=3.255976E-4 C/DBAR FOR S=40 (IPSS-78),
%C T=40 DEG C,P0=10000 DECIBARS
%      IMPLICIT REAL*8 (A-H,O-Z)

      DS = S - 35.0 ;
      ATG = (((-2.1687E-16*T+1.8676E-14).*T-4.6206E-13).*P+...
((2.7759E-12*T-1.1351E-10).*DS+((-5.4481E-14*T+8.733E-12).*T-6.7795E-10).*T+...
1.8741E-8)).*P+(-4.2393E-8*T+1.8932E-6).*DS+((6.6228E-10*T-6.836E-8).*T+...
8.5258E-6).*T+3.5803E-5;
