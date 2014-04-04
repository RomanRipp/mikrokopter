function A = bulk(ur,zr,Ta,zt,rh,zq,Pa,Ts,sal,dlw,dsw,nsw)

% BULK  computes sensible and latent heat fluxes and other variables.
%
% A = BULK(ur,zr,Ta,zt,rh,zq,Pa,Ts,sal,dlw,dsw,nsw) computes the following:
%
% A = [ Hs Hl Hl_webb stress U_star T_star Q_star L zetu CD CT CQ RI ]
%
%         Hs      = sensible heat flux INTO ocean [W/m^2]
%         Hl      = latent heat flux INTO ocean [W/m^2]
%         Hl_webb = Webb correction to latent heat flux INTO ocean [W/m^2]
%         stress  = wind stress [N/m^2]
%         U_star  = velocity friction scale [m/s]
%         T_star  = temperature scale [deg C]
%         Q_star  = humidity scale [kg/kg]
%         L       = Monin-Obukhov length [m]
%         zetu    = zr/L
%         CD      = drag coefficient
%         CT      = temperature transfer coefficient (Stanton number)
%         CQ      = moisture transfer coefficient (Dalton number)
%         RI      = bulk Richardson number
%         Dter    = cool-skin temperature difference (optional output) [C]; 
%                   positive if surface is cooler than bulk (presently no 
%                   warm skin permitted by model)
%                    
% Based on the following buoy input data:
%
%           ur     = wind speed [m/s] measured at height zr [m] 
%           Ta     = air temperature [C] measured at height zt [m]
%           rh     = relative humidity [%] measured at height zq [m]
%         i*sh     = specific humidity [kg/kg] measured at height zq [m]
%           Pa     = air pressure [mb]
%           Ts     = sea surface temperature [C]
%           sal    = salinity [psu (PSS-78)]
%                    (optional - only needed for cool-skin)
%           dlw    = downwelling (INTO water) longwave radiation [W/m^2]
%                    (optional - only needed for cool-skin)
%           dsw    = measured insolation [W/m^2]
%                    (optional - only needed for cool-skin)
%           nsw    = net shortwave radiation INTO the water [W/m^2]
%                    (optional - only needed for cool-skin)
%
% where ur, Ta, rh, Pa, Ts, zr, zt, and zq (and optional sal, dlw,
% dsw, and nsw) may be either row or column vectors; and rh, Pa, 
% zr, zt, and zq (and optional sal) may also be fixed scalars.
%
% Output variables are given as column vectors in A:
%
% 1) without cool-skin correction:
%
%   A = [Hs Hl Hl_webb stress U_star T_star Q_star L zetu CD CT CQ RI]
%
% 2) with cool-skin correction: 
%
%   A = [Hs Hl Hl_webb stress U_star T_star Q_star L zetu CD CT CQ RI Dter];
%

% Code follows Edson and Fairall TOGA COARE code (version 2.0), modified 
% to include Rogers' weighting factor for unstable conditions.  Code does
% include gustiness, and assumes that the marine boundary layer height is
% known and constant over time for simiplicity. zr/L is limited to 
% be <=3.0 to ensure that the code converges to nonzero stress and heat 
% flux values for strongly stable conditions.  The bulk Richardson number
% is computed between the sea surface and zr as a diagnostic about whether
% turbulent boundary layer theory is applicable.  Code does not include 
% warm layer effects to modify Ts.  See Fairall et al (1996), J. Geophys. 
% Res., 101, 3747-3764, for description of full TOGA COARE code and 
% comparison with data. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8/19/98: version 1.1 (rewritten by RP to remove inconsistencies in 
%          virtual and real temperatures, improve loop structure, 
%          correct gustiness component of stress computation) 
% 4/9/99: version 1.2 (rewritten by AA to clarify some variable names
%         and include cool-skin effect and Webb correction to latent 
%         heat flux added to output matrix)
% 8/5/99: version 2.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% change to column vectors


ur = ur(:);
Ta = Ta(:);
rh = rh(:);
Ts = Ts(:);
Pa = Pa(:);
zr = zr(:);
zt = zt(:);
zq = zq(:);


M = size(ur,1);

% create vectors for rh, Pa, zr, zt, and zq, if scalars are input
if length(rh)==1 & M>1
  rh=rh*ones(M,1);
end
if length(Pa)==1 & M>1
  Pa=Pa*ones(M,1);
end
if length(zr)==1 & M>1
  zr=zr*ones(M,1);
end
if length(zt)==1 & M>1
  zt=zt*ones(M,1);
end
if length(zq)==1 & M>1
  zq=zq*ones(M,1);
end


if nargin > 8
  % optional cool-skin stuff
  sal = sal(:);
  dlw = dlw(:);
  dsw = dsw(:);
  nsw = nsw(:);
  % create vector for sal if scalar is input
  if length(sal)==1 & M>1
    sal=sal*ones(M,1);
  end
end

% initialize various constants
as_consts;

tol = .001;                         % tolerance on Re changes to make sure soln has converged.

onethird = 1./3;
o61 = 1/eps_air-1;                  % 0.61 (moisture correction for temperature)

visc  = viscair(Ta);                % viscosity

Qsats = Qsat_coeff * qsat(Ts,Pa);   % saturation specific humidity; the Qsat_coeff
                                    % value is set in routine as_consts.m
Q = imag(rh);
ii = find( Q == 0 );
if ~isempty(ii)
    Q(ii) = 0.01 * real(rh(ii));
    Q(ii) = Q(ii).*qsat(Ta,Pa);     % specific humidity of air [kg/kg]
end

T   = Ta + CtoK;                    % convert to K
Tv  = T .* ( 1 + o61*Q );           % air virtual temperature
rho = 100 * Pa ./ (gas_const_R*Tv); % air density
Dt  = ( Ta + 0.0098.*zt ) - Ts;     % adiabatic temperature difference
Dq  = Q - Qsats;                    % humidity difference

% compute initial neutral scaling coefficients
S     = sqrt( ur.^2 + min_gustiness.^2 );
cdnhf = sqrt(cdntc(S,zr,Ta));       % Smith's neutral cd as first guess

z0t   = 7.5*10^(-5);
ctnhf = kappa ./ log(zt./z0t);

z0q   = z0t;
cqnhf = kappa ./ log(zq./z0q);

U_star = cdnhf .* S;      % (includes gustiness)
T_star = ctnhf .* Dt;     % 
Q_star = cqnhf .* Dq;     %

Dter   = 0;
Dqer   = 0;

if nargin > 8
% initial cool-skin thickness guess  
  delta = 0.001*ones(size(Ts));
end

Reu = 0;
Ret = 0;
Req = 0;

% begin iteration loop to compute best U_star, T_star, and Q_star

for iter1 = 1 : 80

    ReuO = Reu; 
    RetO = Ret; 
    ReqO = Req; % Save old values
    
    % Compute Monin-Obukov length (NB - definition given as eqn (7)
    % of Fairall et al (1996) probably wrong, following, e.g.
    % Godfrey and Bellars (1991), JGR, 96, 22043-22048 and original code)

    bs = g * ( T_star .* ( 1 + o61*Q ) + o61*T.*Q_star ) ./ Tv; 
    L  = (U_star.^2) ./ (kappa*bs);

    % set upper limit on zr/L = 3.0 to force convergence under 
    % very stable conditions. Assume that zr, zt and zq comparable.

      index_limit  = ( ( L < zr/3 )  &  ( L > 0 ) );
    L(index_limit) = zr(index_limit) / 3;
    
    zetu = zr ./ L;  % nondimensionalized heights
    zett = zt ./ L;
    zetq = zq ./ L;

    % surface roughness
    z0 = (Charnock_alpha/g) .* U_star.^2 + R_roughness .* visc ./ U_star;

    % compute U_star correction for non-neutral conditions
    cdnhf  = kappa ./ ( log(zr./z0) - psiutc(zetu) );
    U_star = cdnhf.*S;
  
     Reu      = z0 .* U_star ./ visc;   % roughness Reynolds #
    [Ret,Req] = LKB(Reu);               % compute other roughness Reynolds #s

    % compute t and q roughness scales from roughness R#s
    z0t = visc .* Ret ./ U_star;
    z0q = visc .* Req ./ U_star;

    % compute new transfer coefficients at measurement heights
    cthf = kappa ./ ( log(zt./z0t) - psittc(zett) );
    cqhf = kappa ./ ( log(zq./z0q) - psittc(zetq) );

    % compute new values of T_star, Q_star
    T_star = cthf .* ( Dt + Dter );
    Q_star = cqhf .* ( Dq + Dqer );

    % estimate new gustiness
    Ws     = U_star .* (-CVB_depth./(kappa*L)).^onethird;
    wg     = min_gustiness * ones(M,1);
       jj  = find( zetu < 0 );                     % convection in unstable conditions only
    wg(jj) = max(min_gustiness,beta_conv.*Ws(jj)); % set minimum gustiness

    S = sqrt( ur.^2 + wg.^2 );

    if nargin > 8
    % compute cool-skin parameters
      [delta,Dter,Dqer] = cool_skin(sal,Ts-Dter,rho,cpp,Pa, ...
                                    U_star,T_star,Q_star, ...
                                    dlw,dsw,nsw,delta,g,gas_const_R, ...
                                    CtoK,Qsat_coeff);
    end

end % end of iteration loop

ii = ( ( abs(Reu-ReuO) > tol ) | ...
       ( abs(Ret-RetO) > tol ) | ...
       ( abs(Req-ReqO) > tol )       );

if any(ii),
 disp(['Algorithm did not converge for ' int2str(sum(ii)) ' values. Indices are:']);
 disp(find(ii)');
 warning('Not converged!');
end;


% compute latent heat
Le = ( 2.501 - 0.00237*(Ts-Dter) ) * 10^6;

% compute fluxes into ocean
Hs = rho .* cpp .* U_star .* T_star;
Hl = rho .* Le  .* U_star .* Q_star;

% compute transfer coefficients at measurement heights
CD = (U_star ./ S ).^2;
CT =  U_star .* T_star ./ (S.*(Dt + Dter)); % Stanton number
CQ =  U_star .* Q_star ./ (S.*(Dq + Dqer)); % Dalton number

% to compute mean stress, we don't want to include the effects
% of gustiness which average out (in a vector sense).
stress = rho .* CD .* S .* ur;

% compute bulk Richardson number (as a diagnostic) - the "T"
% is probably not quite right - assumes T \ approx. Ts (good enough though)
RI = g * zr .* ( Dt + Dter + o61*T.*(Dq + Dqer) ) ./ (Tv.*S.^2);

% compute Webb correction to latent heat flux into ocean
W = 1.61 .* U_star .* Q_star + ( 1 + 1.61*Q ) .* U_star .* T_star./T; % eqn. 21

Hl_webb = rho .* Le .* W .* Q; % eqn. 22, Fairall et al. (1996), JGR, 101, p3751.

% output array
if nargin > 8
  % output additional cool-skin parameters 
  A=[Hs Hl Hl_webb stress U_star T_star Q_star L zetu CD CT CQ RI Dter];
else
  % otherwise
  A=[Hs Hl Hl_webb stress U_star T_star Q_star L zetu CD CT CQ RI];
end

%*****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function y = psiutc(zet)

% PSIUTC: computes velocity profile function following TOGA/COARE.
% y=PSIUTC(zet) computes the turbulent velocity profile function given 
% zet = (z/L), L the Monin-Obukoff length scale, following Edson and
% Fairall TOGA COARE code (version 2.0) as modified to include Rogers' 
% weighting factor to combine the Dyer and free convection forms for 
% unstable conditions. 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8/28/98: version 1.1
% 8/5/99: version 1.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c13 = 1.0 / 3.0;
sq3 = sqrt(3.0);

% stable conditions
y = -4.7 * zet;

% unstable conditions
  jj = find( zet < 0 );
zneg = zet(jj);

% nearly stable (standard functions)
 x  = ( 1 - 16*zneg ).^0.25;
 y1 = 2.0 * log((1+x)/2) + log((1+x.^2)/2) - 2*atan(x) + pi/2;

% free convective limit
 x  = ( 1 - 12.87*zneg ).^c13;
 y2 = 1.5 * log(c13*(x.^2+x+1)) - sq3*atan((2*x+1)/sq3) + pi/sq3;
	
% weighted sum of the two
 F     = 1 ./ ( 1 + zneg.^2 );
 y(jj) = F .* y1 + (1-F) .* y2;


%*****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function y = psittc(zet)

% PSITTC: computes potential temperature profile following TOGA/COARE.
% y=PSITTC(zet) computes the turbulent potential temperature profile 
% function given zet = (z/L), L the Monin-Obukoff length scale, following 
% Edson and Fairall TOGA COARE code (version 2.0), as modified to use
% Rogers' weighting factor to combine the Dyer and free convective 
% forms for unstable conditions. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8/28/98: version 1.1
% 8/5/99: version 1.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c13 = 1.0 / 3.0;
sq3 = sqrt(3.0);

% stable conditions
y = -4.7 * zet;

% unstable conditions
  jj = find( zet < 0 );
zneg = zet(jj);

% nearly stable (standard functions)
 x  = ( 1 - 16*zneg ).^0.25;
 y1 = 2 * log((1+x.^2)/2);

% free convective limit
 x  = ( 1 - 12.87*zneg ).^c13;
 y2 = 1.5 * log(c13*(x.^2+x+1)) - sq3.*atan((2*x+1)/sq3) + pi/sq3;

% weighted sum of the two
 F     = 1 ./ ( 1 + zneg.^2 );
 y(jj) = F .* y1 + (1-F) .* y2;


%*****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Ret,Req] = LKB(Reu);

% LKB: computes rougness Reynolds numbers for temperature and humidity
% [Ret,Req]=LKB(Reu) computes the roughness Reynolds for temperature
% and humidity following Liu, Katsaros and Businger (1979), J. Atmos.
% Sci., 36, 1722-1735.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8/28/98: version 1.1
% 8/5/99: version 1.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   Ret = .177 * ones(size(Reu));
   Req = .292 * ones(size(Reu));

       jj  = find( ( Reu > .11 ) & ( Reu <= .825 ) );
   Ret(jj) = 1.376 .* Reu(jj).^0.929;
   Req(jj) = 1.808 .* Reu(jj).^0.826;

       jj  = find( ( Reu > .825 ) & ( Reu <= 3 ) );
   Ret(jj) = 1.026 ./ Reu(jj).^0.599;
   Req(jj) = 1.393 ./ Reu(jj).^0.528;

       jj  = find( ( Reu > 3 ) & ( Reu <= 10 ) );
   Ret(jj) = 1.625 ./ Reu(jj).^1.018;
   Req(jj) = 1.956 ./ Reu(jj).^0.870;

        jj = find( ( Reu > 10 ) & ( Reu <= 30 ) );
   Ret(jj) = 4.661 ./ Reu(jj).^1.475;
   Req(jj) = 4.994 ./ Reu(jj).^1.297;

       jj  = find( Reu > 30 );
   Ret(jj) = 34.904 ./ Reu(jj).^2.067;
   Req(jj) = 30.790 ./ Reu(jj).^1.845;

%*****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function q = qsat(Ta,Pa)

% QSAT: computes specific humidity at saturation. 
% q=QSAT(Ta) computes the specific humidity (kg/kg) at satuation at
% air temperature Ta (deg C). Dependence on air pressure, Pa, is small,
% but is included as an optional input.
%
%    INPUT:   Ta - air temperature  [C]
%             Pa - (optional) pressure [mb]
%
%    OUTPUT:  q  - saturation specific humidity  [kg/kg]

% Version 1.0 used Tetens' formula for saturation vapor pressure 
% from Buck (1981), J. App. Meteor., 1527-1532.  This version 
% follows the saturation specific humidity computation in the COARE
% Fortran code v2.5b.  This results in an increase of ~5% in 
% latent heat flux compared to the calculation with version 1.0.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/8/97: version 1.0
% 4/7/99: version 1.2 (revised as above by AA)
% 8/5/99: version 2.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==1,
  as_consts;
  Pa = P_default; % pressure in mb
end;

% original code
% a=(1.004.*6.112*0.6220)./Pa;
% q=a.*exp((17.502.*Ta)./(240.97+Ta))

% as in Fortran code v2.5b for COARE
ew = 6.11210 * ( 1.0007 + 3.46e-6 * Pa ) .* exp(17.502*Ta./(240.97+Ta)); % in mb
q  = 0.62197 * ew ./ ( Pa - 0.378*ew );                         % mb -> kg/kg


%*****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [cdrg,u10] = cdntc(sp,z,Ta)

% CTDTC: computes the neutral drag coefficient following Smith (1988).
% [cdrg,u10]=CDNTC(sp,z,Ta) computes the neutral drag coefficient and 
% wind speed at 10m given the wind speed and air temperature at height z 
% following Smith (1988), J. Geophys. Res., 93, 311-326. 
%
% INPUT:   sp - wind speed  [m/s]
%          z - measurement height [m]
%          Ta - air temperature (optional)  [C] 
%
% OUTPUT:  cdrg - neutral drag coefficient at 10m
%          u10  - wind speed at 10m  [m/s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/8/97: version 1.0
% 8/26/98: version 1.1 (vectorized by RP)
% 8/5/99: version 2.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get constants
as_consts; 

if nargin==2,
  Ta = Ta_default;
end

% iteration endpoint
tol  = .00001;

visc = viscair(Ta);

% remove any sp==0 to prevent division by zero
sp = sp + 0.1 * ( sp == 0 );

% initial guess
ustaro = zeros(size(sp));
ustarn = 0.036 .* sp;

% iterate to find z0 and ustar
ii = ( abs(ustarn-ustaro) > tol );

while any(ii(:))

  ustaro = ustarn;

  z0 = Charnock_alpha/g * ustaro.^2 + R_roughness * visc ./ ustaro;

  ustarn = sp .* kappa ./ log(z./z0);

  ii = ( abs(ustarn-ustaro) > tol );

end

sqrcd = kappa ./ log(10./z0);
cdrg  = sqrcd.^2;

u10   = ustarn ./ sqrcd;

%*****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  vis = viscair(Ta)

% VISCAIR: computes viscosity of air 
% vis=VISCAIR(Ta) computes the kinematic viscosity of dry air as a 
% function of air temperature following Andreas (1989), CRREL Report 
% 89-11.
%
% INPUT:   Ta  -  air temperature  [C]
%
% OUTPUT:  vis  -  air viscosity  [m^2/s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/8/97: version 1.0
% 8/5/99: version 2.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vis = 1.326e-5 * ( 1 + 6.542e-3 * Ta + 8.301e-6 * Ta.^2 - 4.84e-9 * Ta.^3 );


%*****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function as_consts

% AS_CONSTS: returns values of many constants used in AIR_SEA TOOLBOX
% AS_CONSTS: returns values of many constants used in the AIR_SEA
% TOOLBOX.  At end of this file are values of constants from COARE
% to be used when running the test program t_hfbulktc.m 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8/19/98: version 1.1 (contributed by RP)
% 4/7/99: version 1.2 (revised to include COARE test values by AA)
% 8/5/99: version 2.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------- physical constants

g           = 9.8;       % acceleration due to gravity [m/s^2]
sigmaSB     = 5.6697e-8; % Stefan-Boltzmann constant [W/m^2/K^4]
eps_air     = 0.62197;   % molecular weight ratio (water/air)
gas_const_R = 287.04;    % gas constant for dry air [J/kg/K]
CtoK        = 273.16;    % conversion factor for [C] to [K]


% ------- meteorological constants

kappa          = 0.4;    % von Karman's constant
Charnock_alpha = 0.011;  % Charnock constant (for determining roughness length
                         % at sea given friction velocity), used in Smith
                         % formulas for drag coefficient and also in Fairall
                         % and Edson.  use alpha=0.011 for open-ocean and 
                         % alpha=0.018 for fetch-limited (coastal) regions. 
R_roughness   = 0.11;    % limiting roughness Reynolds # for aerodynamically 
                         % smooth flow         
                         
                        
% ------ defaults suitable for boundary-layer studies

%cp           = 1004.7;   % heat capacity of air [J/kg/K]
cpp           = 1004.7;   % heat capacity of air [J/kg/K]
rho_air       = 1.22;     % air density (when required as constant) [kg/m^2]
Ta_default    = 10;       % default air temperature [C]
P_default     = 1020;     % default air pressure for Kinneret [mbars]
psych_default = 'screen'; % default psychmometer type (see relhumid.m)
Qsat_coeff    = 0.98;     % satur. specific humidity coefficient reduced 
                          % by 2% over salt water


% the following are useful in hfbulktc.m 
%     (and are the default values used in Fairall et al, 1996)

CVB_depth     = 600; % depth of convective boundary layer in atmosphere [m]
min_gustiness = 0.5; % min. "gustiness" (i.e., unresolved fluctuations) [m/s]
                     % should keep this strictly >0, otherwise bad stuff
                     % might happen (divide by zero errors)
beta_conv     = 1.25;% scaling constant for gustiness


% ------ short-wave flux calculations

Solar_const = 1368.0; % the solar constant [W/m^2] represents a 
                      % mean of satellite measurements made over the 
                      % last sunspot cycle (1979-1995) taken from 
                      % Coffey et al (1995), Earth System Monitor, 6, 6-10.
                      
                      
% ------ long-wave flux calculations

emiss_lw = 0.985;     % long-wave emissivity of ocean from Dickey et al
                      % (1994), J. Atmos. Oceanic Tech., 11, 1057-1076.

bulkf_default = 'berliand';  % default bulk formula when downward long-wave
                             % measurements are not made.

%*************************************************

n = whos;
n = {n.name};

for nn = n(:)'
    assignin('caller',nn{1},eval(nn{1}));
end

%*************************************************
 
% ------ constants used for COARE; to use simply delete the %     
% g           = 9.7803;   % acceleration due to gravity [m/s^2]
% sigmaSB     = 5.67e-8;  % Stefan-Boltzmann constant [m^2/K^4] 
% gas_const_R = 287.1;    % gas constant for dry air [J/kg/K]
% cp          = 1004.67;  % heat capacity of air [J/kg/K]
% beta_conv   = 1.20;     % scaling constant for gustiness  
% emiss_lw    = 0.97;     % long-wave emissivity                           
