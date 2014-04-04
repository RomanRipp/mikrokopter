function [U,V,GP] = gp2uv3d(GP,lon,lat,pref,P)
%GP2UV3D Calculation of geostrophic velocities.
%       (works 3-dimensional)
%       [U,V]=GP2UV3d(GP,Lon,Lat,pref,P)
%     
%       Units:
%
%       GeoPotAnomaly    GP     m^2/s^2  M x N     dyn cm 
%       Longitude        Lon    degrees  M x N x 1 or Length N 
%       Latitude         Lat    degrees  M x N x 1 or Length M
%       reference level  pref   dBar     (optional)
%       Pressure         P      dBar     (optional)
%       U-velocity       U      m/s
%       V-velocity       V      m/s
%
%       Velocities with absolute > 5 (m/s)  will be setted to NaN.
% 
%  If Pref negative ,it will be calculated 
%  from maximum, finite Element of each 3. Dimension in P.
%
%  The OutPut GP is GP with Correction of ReferenceLevel at pref.
%

%  Christian Begler




if nargin <= 3
 pref = nan;
else
  if nargin <= 4
   error(['Require P if Pref is given.'])
  else
   pref = pref(1);
  end
end

[d1,d2,d3] = size(GP);
 
dim = 3;




if ~isnan(pref)
 p_si = size(P); 
 if p_si(1:2) == [1 1]
  [hilf1,hilf2,P]=meshgrid(GP(1,:,1),GP(:,1,1),P);                              
    clear hilf1 hilf2
 end 

 d12   = d1*d2;
 d123  = d1*d2*d3;
 ind1 = reshape( (1:d1)'*ones(1,d2) , d12 , 1 );
 ind2 = reshape( ones(d1,1)*(1:d2) , d12 , 1 );                                  
 
if pref >= 0                                                                    
 [hilf,k] = min(abs(P-pref),[],dim);                                            
        k = sub2ind([d1 d2 d3],ind1,ind2,reshape(k,d12,1));                     
                                                                                
 pref = round(min(P(k)));                                                       
  disp(['    GP2UV3D: Reference level set to ',num2str(pref),' dBar.'])        
                                                                                
 [hilf,k] = min(abs(P-pref),[],dim) ; clear hilf              
                  
  kk = sub2ind([d1 d2 d3],ind1,ind2,reshape(k,d12,1)); clear k
  kk = reshape(kk , d1 , d2);
  GP = GP - GP(kk(:,:,ones(1,d3))); clear P kk
 
elseif pref < 0 
% pref < 0 

 [Pmax,k]=max(P,[],dim);
    Pmax = Pmax(:,:,ones(1,d3));
[hilf,k] = min(abs(P-(Pmax+pref)),[],dim); clear hilf

  kk = sub2ind([d1 d2 d3],ind1,ind2,reshape(k,d12,1));
  kk = reshape(kk , d1 , d2);
  GP = GP - GP(kk(:,:,ones(1,d3))); clear P kk

  no =  1 - min(min(k)); % RowMoving  upward
  nu = d3 - max(max(k)); % RowMoving downward 
  nn = [ no : nu ] ; 
 lnn = length(nn);
  k = k(:,:,ones(1,lnn));
 [hilf1,hilf2,nn] = meshgrid((1:d2),(1:d1),nn);
 clear hilf1 hilf2 
  k  = nn+k; clear nn 

 d12n = d12 * lnn;
 ind1 = reshape(ind1*ones(1,lnn),d12n,1);
 ind2 = reshape(ind2*ones(1,lnn),d12n,1);
 k   = sub2ind([d1   d2   d3],ind1,ind2,reshape(k ,d12n,1));
  
  clear ind1 ind2 pp ppn

  k  = reshape(k ,d1  ,d2  ,lnn);

end 

end



% Make Lat, Lon Grid's
 lat_si = size(lat);
 if any(lat_si(1:2)==1)
  lat=meshgrid(lat(:,:,1),GP(1,:,1))';
 end

 lon_si = size(lon);
 if any(lon_si(1:2)==1)
  lon=meshgrid(lon(:,:,1),GP(:,1,1));
 end

 om = 2*pi/24/3600;
 ff = 2*om*sin(lat*pi/180);
  

L = 60*1852;  % Grad --> m

lon = L*lon.*cos(lat*pi/180);
lat = L*lat;

 size(lon),size(lat),size(GP)
if pref >= 0  |  isnan(pref)
 [V,U]=gradnan(GP,lon(:,:,ones(d3,1)),lat(:,:,ones(d3,1)),ones(size(GP)));
else
 U=nan*zeros(d1,d2,d3);
 V=U;
 lon = lon(:,:,ones(d3,1));
 lat = lat(:,:,ones(d3,1));
 [V(k),U(k)]=gradnan(GP(k),lon(k),lat(k));
end

 U = -U ./ ff(:,:,ones(d3,1));
 V =  V ./ ff(:,:,ones(d3,1));


if 0
k = find(abs(U.^2+V.^2) > 25) ;
U(k) = nan*k;
V(k) = nan*k;
end
