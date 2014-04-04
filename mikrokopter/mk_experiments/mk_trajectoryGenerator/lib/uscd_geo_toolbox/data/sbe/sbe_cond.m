function c = sbe_cond(c,t,p,pref,coef);

% SBE_COND Compressibility Compensation of Sea-Bird Conductivity Sensors
%
% C = SBE_COND(C,T,P,PREF,COEF)
%
% Correct conductivity of Sea-Bird-Instruments (MicroCat/SeaCat) 
%  without build-in pressure sensor.
%
% If Sea-Bird-Instrument are deployed without build-in pressure sensors
% a reference pressure 'pref' has to be chosen according to the
% depth where the pressure sensor is moored. Conductivity itself
% is pressure dependend for it depends on the geometric layout
% of the conductivity cell. The latter is changed due to changing
% pressures. 
%
%  C  Conductivity
%  T  Water Temperature [degC]
%  P  Pressure          [dbar]
%
%  PREF  setup reference pressure [dbar]
%  COEF  calibration coefficients [ CTCOR  CPCOR ]
%
% defaults:  PREF = 0dbar
%            COEF =  [ 3.25e-06  -9.57e-08 ] = [ CTCOR  CPCOR ]
%
% Correction:
%
% C = C * ( 1 + CTCOR * T + CPCOR  * PREF ) / ( 1 + CTCOR * T + CPCOR  * P )

%
% Reference: 
%
% Seabird Electonics Application Note 10,
% "Compressibility Compensation of Sea-Bird Conductivity Sensors",
% http://www.seabird.com/application_notes/AN10.htm
%
% f = v / 256 / 1000;  % [Volt] ---> [kHz]
%
%           ( G + H * f^2 + I * f^3 + J * f^4 )
% c = 10 *  -----------------------------------  [mS/cm]
%              ( 1 + CTCOR * T + CPCOR  * P )
%

Nin = nargin;

if Nin < 4
   pref = 0;
end

if Nin < 5
   coef = [ 3.250000e-06  -9.570000e-08 ];  % [ CTCOR CPCOR ]
end

sc = size(c); pc = prod(sc);
st = size(t); pt = prod(st);
sp = size(p); pp = prod(sp);

if ~( ( isequal(sc,st) | ( pt == 1 ) | ( pc == 1 ) ) & ...
      ( isequal(sc,sp) | ( pp == 1 ) | ( pc == 1 ) ) & ...
      ( isequal(st,sp) | ( pt == 1 ) | ( pp == 1 ) )       )
    error('Size of C, T and P must be agree.');
end

if ~( prod(size(pref)) == 1 )
    error('PREF must be a single Value');
end

if ~( prod(size(coef)) == 2 )
    error('COEF must have 2 Elements: [ CTCOR CPCOR ]');
end

c = c .* ( 1 + coef(1) * t + coef(2) * pref ) ./  ...
         ( 1 + coef(1) * t + coef(2) * p );
