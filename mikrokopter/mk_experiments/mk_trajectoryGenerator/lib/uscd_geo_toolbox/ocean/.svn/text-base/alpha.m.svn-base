function [spv,k,steran] = alpha(p,t,s)

% ALPHA  EQUATION OF STATE FOR SEAWATER PROPOSED BY JPOTS 1980
%
% [SPV,K,SterAn] = ALPHA(P,T,S)
%
% input  :	P		: pressure [dbar]
%		T		: temperature [degC] IPTS-78
%		S		: salinity [psu]     IPSS-78
% 
% output :	SPV		: specific volume [m^3/kg]
%		K		: secant bulk modulus [dbar]
%		SterAn		: steric anomaly [m^3/kg]
%
% check values :	alpha(10000,40,40) = 9.435561e-4 m^3/kg
%
%      P could be a [ M by 1 ] Vector if T,S 2-dimensional,                    
%           or a [ 1 x 1 x N ] Vector if T,S 3-dimensional,                    
%                   (the last Matlab 5.# only)                                 
%
% version 1.1.0		last change 06.09.1995

% reference : Landolt-Boernstein V/3a pp 237-242
%       IFM KIEL        T.MUELLER
%       18/02/92, C. Mertens, IfM Kiel, changed to Matlab
% revised header, added bulk-modulus, steric anomaly
%	G.Krahmann, IfM Kiel, Sep 1995


if length(p) == 1
 p = p + 0*t;
end

if length(size(t)) == 3                                                        
                                                                               
 p_si = size(p);                                                               
 if p_si(1:2) == [1 1]                                                         
  ok=1;
  try
    [hilf1,hilf2,p]=meshgrid(t(1,:,1),t(:,1,1),p);
  catch
   ok=0;
  end
  if ~ok
   error(lasterr)
  end
    clear hilf1 hilf2
 end                                                                           
                                                                               
elseif length(size(t)) == 2                                                     
                                                                               
 [m,n] = size(p) ;                                                             
 if n == 1 ,                                                                   
         [m,n] = size(t) ;                                                     
         p = p*ones(1,n) ;                                                     
 end                                                                           
                                                                               
end  


p = p/10 ;
sr = sqrt(abs(s)) ;
%pure water density at atm pressure
rhow = ((((6.536332E-9*t - 1.120083E-6).*t +1.001685E-4).*t - 9.095290E-3).*t ...
         + 6.793952E-2).*t + 999.842594 ;

%seawater density at atm pressure
r1 = (((5.3875E-9*t - 8.2467E-7).*t + 7.6438E-5).*t - 4.0899E-3).*t ...
      + 8.24493E-1 ;
r2 = (-1.6546E-6*t + 1.0227E-4).*t - 5.72466E-3 ;
r3 = 4.8314E-4 ;
rho0 = (r3.*s + r2.*sr + r1).*s + rhow ;
%specific volume at atm pressure
spv = 1 ./ rho0 ;

%compute secant bulk modulus k(p,t,s)
e = (9.1697E-10*t + 2.0816E-8).*t -9.9348E-7 ;
bw = (5.2787E-8*t - 6.12293E-6).*t + 8.50935E-5 ;
b = bw + e.*s ;
d = 1.91075E-4 ;
c = (-1.6078E-6*t - 1.0981E-5).*t + 2.2838E-3 ;
aw = ((-5.77905E-7*t + 1.16092E-4).*t + 1.43713E-3).*t + 3.239908 ;
a = (d.*sr + c).*s + aw ;
b1 = (-5.3009E-4*t + 1.6483E-2).*t + 7.944E-2 ;
a1 = ((-6.1670E-5*t + 1.09987E-2).*t -0.603459).*t + 54.6746 ;
kw = (((-5.155288E-5*t + 1.360477E-2).*t - 2.327105).*t + 148.4206).*t ...
       + 19652.21 ;

%compute k(0,t,s)
k0 = (b1.*sr + a1).*s + kw ;
%evaluate k(p,t,s)
k = (b.*p + a).*p + k0 ;
spv = spv.*(1-p./k) ;

% convert k to dbar for output
k = k*10;

if nargout < 3
   return
end

% compute steric anomaly
if ~all( isnan(t(:)) | isnan(s(:)) )
  steran=spv-alpha(p,0,35);
else
  steran=zeros(size(p));
end
