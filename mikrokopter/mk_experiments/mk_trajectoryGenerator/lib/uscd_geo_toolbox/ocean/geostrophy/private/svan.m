
function [sva,sigma] = svan(s,t,p)
%SVAN Specific volume anomaly.
%  [SVA,SIGMA] = SVAN(S,T,P) returns the specific volume anomaly 
%  (steric anomaly) and density anomaly sigma (kg/m³) as a function of 
%  salinity s (PSS-78), temperature t (øC, IPTS-68) and pressure p (dbar),
%  based on the equation of state for seawater (1980) and the practical
%  salinity scale (1978).
%
%  Checkvalue:
%  svan = 981.30210 (e-8 m³/kg) and sigma = 59.82037 kg/m³ for s = 40, 
%  t = 40 øC and p = 10000 dbar
% 
%  References:
%  Millero, F. J., C.-T. Chen, A. Bradshaw and K. Schleicher (1980).
%    A new high-pressure equation of state for seawater. 
%    Deep-Sea Res., 27(a), 255-264.
%  Millero, F. J. and A. Poisson (1981). Summary of data treatment
%    for the Unesco one atmosphere equation of state for seawater. 
%    Deep-Sea Res., 28(a), 625-629.
%  Unesco (1981). Background papers and supporting data on the international
%    equation of state of seawater 1980. Unesco Tech. Pap. in Mar. Sci.,
%    No. 38, 192 pp.
%  Unesco (1983). Algorithms for computation of fundamental properties
%    of seawater. Unesco Tech. Pap. in Mar. Sci., No. 44, 54 pp.

%  Christian Mertens, IfM Kiel
%  $Revision: 1.1 $ $Date: 1996/03/14 16:39:19 $

dr350 = 28.106331;
% convert pressure to bars and take square root salinity
p = p/10;
sr = sqrt(abs(s));
% pure water density at atmospheric pressure
% Bigg, P. H. (1967). Br. J. Applied Physics, 8, 521-537
r1 = ((((6.536332e-9*t - 1.120083e-6).*t + 1.001685e-4).*t ...
  - 9.095290e-3).*t + 6.793952e-2).*t - 28.263737;
% seawater density at atmospheric pressure,
% coefficients involving salinity
r2 = (((5.3875e-9*t - 8.2467e-7).*t + 7.6438e-5).*t - 4.0899e-3).*t + 8.24493e-1;
r3 = (-1.6546e-6.*t + 1.0227e-4).*t - 5.72466e-3;
r4 = 4.8314e-4;
% international one-atmosphere equation of state of seawater
sig = (r4.*s + r3.*sr + r2).*s + r1;
% specific volume at atmospheric pressure 
r3500 = 1028.1063;
v350p = 1/r3500;
sva = -sig.*v350p./(r3500 + sig);
sigma = sig + dr350;

%
% New high pressure equation of state for seawater.
%

% compute compression terms 
e = (9.1697e-10*t + 2.0816e-8).*t - 9.9348e-7;
bw = (5.2787e-8*t - 6.12293e-6).*t + 3.47718e-5;
b = bw + e.*s;

d = 1.91075e-4;
c = (-1.6078e-6*t - 1.0981e-5).*t + 2.2838e-3;
aw = ((-5.77905e-7*t + 1.16092e-4).*t + 1.43713e-3).*t - 0.1194975;
a = (d.*sr + c).*s + aw;

b1 = (-5.3009e-4*t + 1.6483e-2).*t + 7.944e-2;
a1 = ((-6.1670e-5*t + 1.09987e-2).*t - 0.603459).*t + 54.6746;
kw = (((-5.155288e-5*t + 1.360477e-2).*t - 2.327105).*t + 148.4206).*t - 1930.06;
k0 = (b1.*sr + a1).*s + kw;

% evaluate pressure polynomial, k equals the secant bulk modulus of seawater
% dk = k(s,t,p) - k(35,0,p), k35 = k(35,0,p) 
dk = (b.*p + a).*p + k0;
k35 = (5.03217e-5*p + 3.359406).*p + 21582.27;
gam = p./k35;
pk = 1 - gam;
sva = sva.*pk + (v350p + sva).*p.*dk./(k35.*(k35 + dk));
v350p = v350p.*pk;
% compute density anomaly with respect to 1000 kg/m³
% dr350 = density anomaly at 35 (PSS-78), 0 deg. C and 0 decibars
% dr35p = density anomaly at 35 (PSS-78), 0 deg. C, press. variation
% dvan = density anomaly variations involving spec. vol. anomaly
dr35p = gam./v350p;
dvan = sva./(v350p.*(v350p + sva));
i = find(p ~= 0);
sigma(i) = dr350 + dr35p(i) - dvan(i);

% scale specific volume anomaly to normally reported units
sva = sva*1e8;






