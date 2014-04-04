
function [v,pav,latv,lonv] = geovelsva(sva,p,lat,lon)

% GEOVELSVA          v = geovelsva(sva,p,lat,lon)
%
%         Calculation of geostrophic velocity. GEOVEL returns geostrophic 
%         velocity as a function of in spec. volume anomaly SVA relative 
%         to the surface.
%         The columns of SVA and P
%         are assumed to be individual profiles located at the positions given 
%         in the vectors lat and lon.
%
%         [V,P,LAT,LON] = GEOVELSVA(SVA,P,LAT,LON) returns pressure matrix 
%         corresponding to v and positions of velocity profiles.
%
%         input  : sva     spec. volume anomaly        [m^3/kg]
%   
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
%  Torsten Kanzow 23.9.98, adapted from function geovel.m from: 
%  Christian Mertens, IfM Kiel
%  $Revision$ $Date$
%  04.03.1998, d.kieke, improved header
%
%

%--- SET OUTPUT LATITUDE AND LONGITUDE -----------------------------------------

m    = length(lat);
latv = 0.5*(lat(1:m-1) + lat(2:m));
lonv = 0.5*(lon(1:m-1) + lon(2:m));

%%%[m,n] = size(s);

[m,n] = size(sva);  
lat   = lat(:);
lon   = lon(:);

%--- SPECIFIC VOLUME ANOMALY ---------------------------------------------------

%spv0 = alpha(p,zeros(m,n),35*ones(m,n)); %spez. Volume of Watercolumn: T=0C,S=35PSU
sva  = sva;   

%--- GEOPOTENTIAL ANOMALY ------------------------------------------------------

gpa = [zeros(1,n); cumsum(0.5*(sva(1:m-1,:) + sva(2:m,:)).*diff(p))];

%--- CORIOLIS PARAMETER TIMES DISTANCE BETWEEN STATION PAIRS -------------------

omega = 7.27e-5 ;
lf    = 2*omega*sin(pi/180*0.5*(lat(1:n-1) + lat(2:n))).*dist(lat,lon);

%--- GEOSTROPHIC VELOCITY (M/S) ------------------------------------------------

v = 1e-4*(diff(gpa') ./ lf(:,ones(1,m)))';

%--- SET OUTPUT PRESSURE MATRIX ------------------------------------------------

pav = 0.5*(p(:,1:n-1) + p(:,2:n));



