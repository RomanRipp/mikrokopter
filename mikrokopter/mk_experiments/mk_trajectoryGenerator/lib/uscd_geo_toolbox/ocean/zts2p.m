function  [p,d] = zts2p(z,t,s,lat,N)

% ZTS2P  Calculates Pressure from Depth, Temp, Salinity
% 
%  P = ZTS2P( Z , T , S , LAT );
%
% Formula:  dP = (1000+sigma(P,T,S)) * g * dZ   
%            P = cumsum( P0 ; dP );  
%
%  The Gravity    g  calculates from GRAVIT(LAT)
%  The StartValue P0 calculates from  Z2P80(Z,LAT)
%  The Pressure P in sigma iterates, Start: P = Z2P80(Z,LAT)
%
%  [ P , STD ] = ZTS2P( Z , T , S , Lat , N );
%
%   gives the Number N of Loops for the Iteration (default: 4 )
%   and returns the Standard deviation STD of the Difference 
%    of the Pressure P  to the Pressure in the Step before.
%
%  If  T and S are a 2-dimensional Array; [ z by x ] 
%      Z can be a Vector with Length of 1. Dimension of T,
%    LAT can be a Vector with Length of 2. Dimension of T.
%
%  If  T and S are a 3-dimensional Array: [ y by x by z ] 
%      Z can be a Vector with Length of 3. Dimension of T,
%    LAT can be a Vector with Length of 1. Dimension of T.
%
%

p = [];
d = [];


Nin = nargin;
Nout = nargout;

if Nin < 4
  lat = 45;
end

if Nin < 5
  N = 4;
end

%---------------------------------------------------
% Check Size of T and S

if ~isequal(size(t),size(s))
    error('Size of T and S must be the same.');
end

sd = size(t);
nd = size(sd,2);  % NDIMS

if nd > 3
   error('NDIMS of T and S must be 2 or 3.');
end

%---------------------------------------------------
% Check Size of Z

sz = size(z);
pz = prod(sz);

if ~isequal(sz,sd) 

    if any( sz == 1 )

       z = z(:);

       if     ( nd == 2 )  &  ( pz  == sd(1) )
          z = z(:,ones(1,sd(2)));
       elseif ( nd == 3 )  &  ( pz  == sd(3) )
          z = permute(z,[2 3 1]);
          z = z(ones(1,sd(1)),ones(1,sd(2)),:);
       else
          error('Length of Vector Z must match Length of 1. or 3. Dimension of T.');
       end

    else

        error('Size of Matrice Z must be equal to Size of T.');

    end

end
   
%---------------------------------------------------
% Check LAT

sl = size(lat);
pl = prod(sl);

if ~( isequal(sl,sd)  |  ( pl == 1 ) )

    if any( sl == 1 ) 

       lat = lat(:);

       if     ( nd == 2 )  &  ( pl  == sd(2) )
          lat = permute(lat,[2 1]);
       elseif ( nd == 3 )  &  ( pl  == sd(1) )
          lat = lat(:,ones(1,sd(2)));
       else
          error('Length of Vector LAT must match Length of 2. or 1. Dimension of T.');
       end

    else

        error('Size of Matrice LAT must be equal to Size of T.');

    end

end

%***************************************************

%---------------------------------------------------
% Permute and Reshape if 3D

if nd == 3

   perm = [ 3 1 2 ];

   z  = reshape(permute(z,perm),sd(3),sd(1)*sd(2));
   t  = reshape(permute(t,perm),sd(3),sd(1)*sd(2));
   s  = reshape(permute(s,perm),sd(3),sd(1)*sd(2));

   if ~( prod(size(lat)) == 1 )
       lat  = reshape(permute(lat,perm),1,sd(1)*sd(2));
   end

end

if ( prod(size(lat)) == 1 )
   lat = lat * ones(1,size(t,2));
end

%---------------------------------------------------

is_nan = find( isnan(z) | isnan(t) | isnan(s) );

%---------------------------------------------------
% Initialisation

% dP = Rho * g * dZ 

   g = gravit(lat);

  n1 = size(t,1);

  ind = ( 2 : n1 );

 % First Initialisation of Pressure

  p = z2p80(z,lat(ones(1,n1),:));

 % Mean Values of T , S 

  t = ( t(ind-1,:) + t(ind,:) ) / 2;
  s = ( s(ind-1,:) + s(ind,:) ) / 2;

 % DepthDifferential

  z =  ( z(ind,:) - z(ind-1,:) );


%---------------------------------------------------
% Loop

for ii = 1 : N
 
 if ( Nout == 2 ) & ( ii == N )
    p1 = p;
 end
  
 % Mean Values of P 
 p(ind,:) = ( p(ind-1,:) + p(ind,:) ) / 2;

 % Pressure,  N/m^2  -->  dbar

 p(ind,:) = g(ones(1,n1-1),:).*(1000+sigma(p(ind,:),t,s)).*z / 1e4;

 p = cumsum( p , 1 );

end

%---------------------------------------------------

if Nout < 2
   p(is_nan) = NaN;
   if nd == 3
      p = ipermute(reshape(p,sd(3),sd(1),sd(2)),perm);
   end
   return
end

%---------------------------------------------------
% STD

d = NaN * ones(1,size(p,2));

for ii = 1 : size(p,2)

    dp = p(:,ii)-p1(:,ii);

    jj = find(~isnan(dp));
    if ~isempty(jj)
         d(ii) = std(dp(jj));
    end
end

p(is_nan) = NaN;

if nd == 3
   p = ipermute(reshape(p,sd(3),sd(1),sd(2)),perm);
   d = reshape(d,sd(1),sd(2));
end
 

%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function p = z2p80(z,l);

% Z2P80 Computes pressure given the depth at some latitude
%
%          P = Z2P80(D,LAT) gives the pressure P (dbars) at a depth D (m)
%
%          at some latitude LAT (degrees).
%
%          This probably works best in mid-latitude oceans, if anywhere!
%
%          Ref: Saunders, "Practical Conversion of Pressure to Depth",
%              J. Phys. Oceanog., April 1981.
% 

%         I copied this directly from the UNESCO algorithms.


% CHECK VALUE: P80=7500.004 DBARS;FOR LAT=30 DEG., DEPTH=7321.45 METERS

      l = 5.92e-3 + 5.25e-3 * sin( abs(l*pi/180) ).^2;
      z = sqrt( (1-l).^2 - 8.84e-6 * z  );
      p = (  1 - l - z ) / 4.42e-6;

%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function g = gravit(lat);

% GRAVIT  Returns the Gravitation Acceleration from Latitude
%
% GRAVIT( Latitude )
%
%-------------------------------------------------------------------------------
%
% gravit        - Calcul de l'acceleration de gravite en fonction
%                 de la latitude en degres decimaux. (formule GRS-80)
%
%-------------------------------------------------------------------------------
% Version:
% -------
%  1.01 Criation (d'aprhs gravi8, chaine hydro)          14/06/94 F. Gaillard
% 
%-------------------------------------------------------------------------------

lat   = lat * pi/180;

lat2  = sin( 2 * lat );
lat   = sin( 1 * lat );

ge     =  9.780318;
b2     =  0.530244e-02;
b4     = -0.585000e-05;

g = ge * ( 1 + b2 * lat.*lat + b4 * lat2.*lat2 );

%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function sig = sigma(p,t,s)

%	function sig=sigma(p,t,s)
%	
% computation of density of seawater 
% referenced to arbitrary pressures 
% based on 'alpha.m'
%
% input  :	p		: pressure [dbar]
%		t		: in situ temperature [degrees Celsius]
%		s		: salinity [psu]
%
% output :	sig		: density of seawater at pressure P (adiabatic)
%				  [kg/m^3]
%
% check values :	sigma(0,40,40) = 21.6788   kg/m^3
%               	sigma(0, 0,35) = 28.106331 kg/m^3
%
%      P could be a [ M by 1 ] Vector if T,S 2-dimensional,                    
%           or a [ 1 x 1 x N ] Vector if T,S 3-dimensional,                    
%                   (the last Matlab 5.# only)                                 
%
%
% version 1.1.0		last change 01.09.1995

% modified from SIGMATH, Uwe Send, March 1995
% optional without Pref		G.Krahmann, IfM Kiel, Sep 1995


xats=alpha(p,t,s);

sig=1./xats-1000;

%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [spv,k,steran] = alpha(p,t,s)

% function [spv,k,steran] = alpha(p,t,s)
%
% EQUATION OF STATE FOR SEAWATER PROPOSED BY JPOTS 1980
%
% input  :	p		: pressure [dbar]
%		t		: temperature [degrees Celsius]
%		s		: salinity [psu]
% 
% output :	spv		: specific volume [m^3/kg]
%		k		: secant bulk modulus [dbar]
%		steran		: steric anomaly [m^3/kg]
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
