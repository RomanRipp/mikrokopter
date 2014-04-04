function [ct,pv,lav,gpa1,gpa2] = geo_transp(p1,t1,s1,la1,p2,t2,s2,la2,pr)

% function   [ct,pv,lav,gpa1,gpa2] = geo_transp(p1,t1,s1,la1,p2,t2,s2,la2,pr)
%
%         Calculation of geostrophic transport. GT returns geostrophic
%         transport as a function of salinity S1/2, temperature T1/2 (0C, IPTS-68) 
%         and  pressure P1/2 (dbar), relative to given pressure level PR. 
%         The columns of S1/2, T1/2, and P1/2 are assumed to be individual profiles 
%         located at the positions given in the points la1/2 and lo1/2.
%
%         WARNING:
%         v - v(pr) is calculated
%
%
%
%         input  : s1,s2     salinity time series (mooring1,mooring2)     [psu]
%                            row index = pressure, column index = time
%                  t1,t2     temperature time series (mooring1,mooring2)  [deg C]
%                  p1,p2     pressure grid (mooring1,2) , can be time  
%                            series or simply column vector, if t/s time  
%                            series are on a static pressure grid         [dbar]
%                  la1,la2   latitude of mooring 1,2                      [deg]
%                  pr        reference pressure 
%             
%         output : ct        cumulative transport [Sv = 1e6 m^3/s] from velocities
%                            rel. to pr. Then integrated from bottom to top
%                  pv        pressure on velocity grid [dbar], optional
%                  lav       latitude between profiles [deg], optional
%                  gpa1,gpa2 geopotential anomaly [m^2/s^2] rel. to surface 
%                 
%         uses   : svan.m
%
%
%
%
% T.Kanzow

loop = 'y';

%--- SET OUTPUT LATITUDE AND LONGITUDE -----------------------------------------

lav = 0.5*(la1 + la2);

[m,n] = size(s1);

if size(p1,2) == 1
   p1 = p1*ones(1,n);
   p2 = p2*ones(1,n);
   loop = 'n'; 
end
%--- SPECIFIC VOLUME ANOMALY ---------------------------------------------------

sva1 = svan(s1,t1,p1) * 1e-8;  %unit m^3/s
sva2 = svan(s2,t2,p2) * 1e-8;

%--- GEOPOTENTIAL ANOMALY ------------------------------------------------------

if loop == 'n'

  gpa1 = cumtrapz(p1(:,1)*1e4,sva1);
  gpa2 = cumtrapz(p2(:,1)*1e4,sva2);

  %----------- reference level ------------------------

  I    = find(abs(p1(:,1)-pr)==min(abs(p1(:,1)-pr)));
  I    = I(1);
  gpa1r = gpa1 - ones(m,1)*gpa1(I,:);    

  I    = find(abs(p2(:,1)-pr)==min(abs(p2(:,1)-pr)));
  I    = I(1);
  gpa2r = gpa2 - ones(m,1)*gpa2(I,:);

else

  for  i = 1 : n,
 
    gpa1(1:m,i) = cumtrapz( p1(:,i)*1e4,sva1(:,i) );
    gpa2(1:m,i) = cumtrapz( p2(:,i)*1e4,sva2(:,i) );

    %---------- reference level ---------------------

    I           = find(abs(p1(:,i)-pr)==min(abs(p1(:,i)-pr)));
    I           = I(1);
    gpa1r(1:m,i) = gpa1(:,i) - ones(m,1)*gpa1(I,i);

    I           = find(abs(p2(:,i)-pr)==min(abs(p2(:,i)-pr)));
    I           = I(1);
    gpa2r(1:m,i) = gpa2(:,i) - ones(m,1)*gpa2(I,i);

  end

end  % end loop if 

%--- CORIOLIS PARAMETER -------------------

omega = 7.27e-5 ;
f    = 2*omega*sin(pi/180*lav);

%--- CUMULATIVE TRANSPORT ------------------------------------------------

z  = - 0.5 *  ( p2z(p1,la1) + p2z(p2,la2) ); % depth 

if loop == 'n' 
  ct = - flipud( cumtrapz(flipud(z(:,1)), flipud(gpa2r - gpa1r)) ) / f;
else
  for i = 1 : n,
    ct(1:m,i) = - flipud( cumtrapz(flipud(z(:,i)), flipud(gpa2r(:,i) - gpa1r(:,i))) ) / f;
  end
end
ct = ct ./ 1e6; % Unit: Sverdrup
  
%--- SET OUTPUT PRESSURE MATRIX ------------------------------------------------

pv = 0.5*(p1+p2);
 


