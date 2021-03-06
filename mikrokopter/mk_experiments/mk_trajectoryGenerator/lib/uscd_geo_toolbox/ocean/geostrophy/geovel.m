
function [v,pav,latv,lonv] = geovel(s,t,p,lat,lon)

% GEOVEL          v = geovel(s,t,p,lat,lon)
%
%         Calculation of geostrophic velocity. GEOVEL returns geostrophic 
%         velocity as a function of salinity S, temperature T (�C, IPTS-68) and
%         pressure P (dbar), relative to the surface. The columns of S, T, and P
%         are assumed to be individual profiles located at the positions given 
%         in the vectors lat and lon.
%
%         WARNING:
%         v - v_surf is calculated 
%          
%
%         [V,P,LAT,LON] = GEOVEL(S,T,P,LAT,LON) returns pressure matrix 
%         corresponding to v and positions of velocity profiles.
%
%         input  : s       salinity                    [psu]
%                  t       in situ temperature IPTS-68 [deg C]
%                  p       pressure                    [dbar]  
%                  lat     latitude                    [deg]
%                  lon     longitude                   [deg]
%
%         output : v       geostrophic velocity relative sea surface [m/sec]
%                  pav     pressure on velocity grid [dbar], optional
%                  latv    latitude between profiles [deg], optional
%                  lonv    longitude betwenn profiles [deg], optional
%  
%
%         uses   : svan.m, dist.m
%
%
%
%  Christian Mertens, IfM Kiel
%  $Revision$ $Date$
%  04.03.1998, d.kieke, improved header
%  15.10.1998, kanzow, warning in header 
%
%

%--- SET OUTPUT LATITUDE AND LONGITUDE -----------------------------------------

m    = length(lat);
latv = 0.5*(lat(1:m-1) + lat(2:m));
lonv = 0.5*(lon(1:m-1) + lon(2:m));

[m,n] = size(s);
lat   = lat(:);
lon   = lon(:);

%--- SPECIFIC VOLUME ANOMALY ---------------------------------------------------

sva = svan(s,t,p);

%--- GEOPOTENTIAL ANOMALY ------------------------------------------------------

gpa = [zeros(1,n); cumsum(0.5*(sva(1:m-1,:) + sva(2:m,:)).*diff(p))];

%--- CORIOLIS PARAMETER TIMES DISTANCE BETWEEN STATION PAIRS -------------------

omega = 7.27e-5 ;
%keyboard
lf    = 2*omega*sin(pi/180*0.5*(lat(1:n-1) + lat(2:n))).*dist2(lat,lon);

%--- GEOSTROPHIC VELOCITY (M/S) ------------------------------------------------

v = 1e-4*(diff(gpa') ./ lf(:,ones(1,m)))';

%--- SET OUTPUT PRESSURE MATRIX ------------------------------------------------

pav = 0.5*(p(:,1:n-1) + p(:,2:n));


