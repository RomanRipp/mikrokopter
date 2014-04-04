function c = oxsat_umol(t,s)

% OXSAT_UMOL   Oxygen equilibrium concentration.
%
%              The equation for the Henry coefficent as a function of
%              temperature and salinity is used to calculate values for
%              unit standard atmospheric concentrations (USAC) in freshwater
%              and seawater in equilibrium with air at total pressure 
%              of 1 atmosphere. 
%
%              usage  : osat = oxsat_umol(pt,s);
%
%              input  : pt    potential temperature  [deg C]
%                       s     salinity               [psu]
%
%              output : osat  oxygen equilibrium concentration [umol/kg]
%
%              uses   : ---
%
%              check values : oxsat_umol(40,40) = 159.53
%
%              validity     : 0 < t < 40 deg C, 0 < s < 40 
%
%              reference    : Benson, B. B. and D. Krause, Jr.,
%                             "The concentration and isotopic fractionation of
%                             oxygen dissolved in freshwater and seawater in
%                             equilibrium with the atmosphere",
%                             Limnol. Oceanogr., 29(3), 1984, p. 620-632 
%
%              
%              version 1.0, Wed Jun 21 19:15:56 MET DST 2000, d.kieke

t  = t + 273.15;

a0 = -135.29996;
a1 =  1.572288e5;    % ... 10^5
a2 = -6.637149e7;    % ... 10^7
a3 =  1.243678e10;   % ... 10^10
a4 = -8.621061e11;   % ... 10^11

b0 = 0.020573;
b1 = -12.142;
b2 = 2363.1;

c = a0 + a1 ./ t + a2 ./ (t.^2) + a3 ./ (t.^3) + a4 ./ (t.^4);

c = c - s .*( b0 + b1./t + b2./(t.^2) );

c = exp(c);



