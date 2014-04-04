function [U,V,GP] = gp2uv2d(GP,lon,lat,pref,P)                   
%GP2UV2D Calculation of geostrophic velocities.               
%       (works 2-dimensional)                                 
%       [U,V,GP]=GP2UV2d(GP,Lon,Lat,pref,P)                      
%                                                             
%       Units:                                                
%                                                             
%       GeoPotAnomaly    GP     m^2/s^2  dyn cm  M x N               
%       Longitude        Lon    degrees  M x 1
%       Latitude         Lat    degrees  M x 1 
%       reference level  pref   dBar     (optional)           
%       Pressure         P      dBar     (optional)           
%       U-velocity       U      m/s                           
%       V-velocity       V      m/s                           
%                                                             
%       Velocities with absolute > 5 (m/s)  will be setted to NaN.
%                                                             
%  If Pref negative ,it will be calculated                    
%  from maximum, finite Element of each Column in P.    
%
%  The OutPut GP is GP with Correction of ReferenceLevel at pref.
%                                                             

%  Christian Begler

if nargin <= 3                                                
 pref = [nan];                                                   
else                                                          
  if nargin <= 4                                              
   error(['Require P if Pref is given.'])                     
  else                                                        
   pref = pref(1);                                            
  end                                                         
end                                                           
                                                              
[d1,d2] = size(GP); 




if ~isnan(pref) 

 if any(size(P)==1)
         P = P(:)*ones(1,d2) ;
 end   

if pref >= 0
 [hilf,k] = minnan(abs(P-pref)); 
       k  = ((1:d2)-1)*d1 + k;

 pref = round(min(P(k)));
  disp(['    GP2UV2D: Reference level set to ',num2str(pref),' dBar.'])

 [hilf,k] = minnan(abs(P-pref)) ;
       k  = ((1:d2)-1)*d1 + k;

  GP = GP - ones(d1,1)*GP(k);

  clear P k 

else
 [hilf,k]=maxnan(P); 
       k  = ((1:d2)-1)*d1 + k;

    Pini = P(k) + pref;
  Pini = ones(d1,1)*Pini(:)';
[hilf,k] = minnan(abs(P-Pini));
      kk  = ((1:d2)-1)*d1 + k;

   GP = GP - ones(d1,1)*GP(kk);        

  clear P kk
  
 no = 1-min(k); % RowMoving  upward
 nu = m-max(k); % RowMoving downward                                           
 [k,nn]= meshgrid(k,[no:nu]');                                               
 
 k = nn + k; clear nn

 k  = ((1:d2)-1)*d1 + k;

end  
end 

                                                              
lat = lat(:)';
lon = lon(:)';

L = 60*1852;  % Grad --> m                                    
if pref>= 0  |  isnan(pref)
 [V]=gradnan(GP,L*lon.*cos(lat*pi/180),1);             
 [U]=gradnan(GP,L*lat,1);
else                                                          
 U=nan*zeros(d1,d2);                                       
 V=nan*zeros(d1,d2);
 [V(k)]=gradnan(GP(k),L*lon(k).*cos(lat(k)*pi/180),1);     
 [U(k)]=gradnan(GP(k),L*lat(k),1);
end                                                           
 om = 2*pi/24/3600;                                           
 U = -U ./ (2*om*sin([ones(size(U,1),1)*lat]*pi/180));                            
 V = +V ./ (2*om*sin([ones(size(V,1),1)*lat]*pi/180));
                                                              
                                                              
k = find(abs(U.^2+V.^2) > 25) ;
U(k) = nan*k; 
V(k) = nan*k;

