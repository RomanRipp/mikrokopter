function depth=p2z(p,lat)
% function  depth=p2z(p,lat)
% !!!!!! USES Z=0 AT P=0  (I.E. NOT 1ATM AT SEA SURFACE)
%	pressure to depth conversion using
%	saunders&fofonoff's method (deep sea res.,
%	1976,23,109-111)
%	formula refitted for alpha(p,t,s) = eos80
%	units:
%		depth         z        meter
%       pressure      p        dbars                                   
%	latitude      lat      deg
%
% lat could be a single value or
%  Matrice with size of depth.
%  also a vector of stationpositions
%
%	checkvalue:
%		depth =       9712.654  m
%	for
%		p     =        10000.     dbars
%		lat   =           30.     deg
%
% version 1.0.0         last change ??
% version 1.0.1         last change by M Hamann, 14.08.2000


      

 p     = p / 10.;

 if size(lat,1) < size(p) | size(lat,2) < size(p)
 if size(lat,1) ==1 & size(lat,2) ==1
  lat = lat;
 else
  lat = lat(:);
  lat = ones(length(p),1)*lat';
 end  
 end  
 x     = sin(lat/57.29578).^2;

 gr    = 9.780318*(1.0+(5.2788e-3+2.36e-5*x).*x)+1.092e-5*p;
 depth = (((-1.82e-11*p+2.279e-7).*p-2.2512e-3).*p+97.2659).*p;
 depth = depth./gr;
