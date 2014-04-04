function [out1,out2,out3,out4] = geovel2d(P,T,S,arg4,arg5,arg6,arg7)
%GEOVEL2D Calculation of geostrophic velocities.
%         (works 2-dimensional)
%       GEOVEL2d(P,T,S, ... )                                  
%                                                              
%       OptionalInputs are:  Lat , Lon , Pref , '[<OutPut>]'   
%          use them only in the order above.                   
%       For Velocities Lat anmd Lon are requiered.             
%                                                              
%       possible OutPutArguments are: [GP,GPa,U, V ]
%       '[<OutPut>]' give the Number of the OutPutArgument(s),
%         you want to get.
%       '[<OutPut>]' is default:     '[ 1  2  3  4 ]'
%
%       The order of the numbers in [<OutPut>]  give the order 
%        of the OutPutArguments.                               
%                                                              
%       Units:                                                 
%                                                              
%       Pressure         P      dBar     M x N  or Length M 
%       Temperature      T      deg C    M x N              
%       Salinity         S      psu      M x N              
%       Longitude        Lon    degrees  N x 1   
%       Latitude         Lat    degrees  N x 1  
%       reference level  pref   dBar                           
%       GeoPotential     GP     m^2/s^2  dyn cm                
%       GeoPotAnomaly    GPa    m^2/s^2                        
%       U-velocity       U      m/s                            
%       V-velocity       V      m/s                            
%
%       Velocities with absolute > 5 (m/s)  will be setted to NaN.
%
%  If Pref negative ,it will be calculated 
%    from maximum, finite Element of each Column in P.
%
%  Pref = NaN , will give only GP or GPa  originally,
%   no Difference to any Level.
%
 

Nout = nargout;
                                                    
Nin = nargin;                                       
if nargin < 3                                       
 error('P,T,S required')                            
end                                                 
                                                    
if length(size(T)) ~= 2                             
                                                    
 error('Use GEOVEL3D for 3-dimensional Fields.')    
                                                    
end                                                 
                                                    
                                                    
Out = [ 1 : Nout ];                                 
pref = 0;                                           
                                                    
for ii = 4:Nin                                      
 eval(['val = arg' int2str(ii) ';'])                
  if isstr(val)                                     
    Out = eval(['[' val ']']) ;                      
   OutNr = ii;                                      
  elseif length(val) == 1                           
   pref = val;                                      
   pNr  = ii;                                       
  end                                               
end                                         

Out = Out(:)';                                      
jj = find([ 1 <= Out  &  Out <= 4   ...             
            &   Out == round(Out) ]);               
Out = Out(jj);                                      
if length(Out) > 4, Out = Out(1:4); end             
if length(Out) > Nout; Out = Out(1:Nout); end       
                                                    
Out = [ Out (max(Out)+1:Nout) ];                           

if isnan(pref)                                 
 jj = find( Out>=3 );                          
 Out(jj) = [];                                 
end                                            
   
                                                    
                                                    
is_uv = [ any(Out==3)  |  any(Out==4) ];            
is_GPa = [ is_uv  |  any(Out==2) ];                  
is_GP = any(Out==1);                               
                                                    
ll = [];                                            
if is_uv                                            
% Lat Lon required                                  
  for ii = 4:Nin                                    
   if ~any(ii==[OutNr pNr])                         
    eval(['val = arg' int2str(ii) ';'])             
    if ~isstr(val)                                  
     ll = [ ll ii ];                                
    end                                             
   end                                              
  end                                               
  if length(ll) < 2                                 
   error(['Lon,Lat requiered.'])                    
  end                                               

  eval(['lat =  arg' int2str(ll(1)) ';' ])
  eval(['lon =  arg' int2str(ll(2)) ';' ])
  
end                                                 
       

[d1,d2]=size(T);

 if any(size(P)==1)
         P = P(:);
         P = P*ones(1,d2) ;
 end   

dim = 1;

is_nan = find( isnan(T) | isnan(S) | isnan(P) );

 dP = diff(P*1e4);                            
 jj = find(isnan(dP));                              
 dP(jj) = 0*jj;                                     
                                                    
                                                    
%specific volume(anomaly) [m^3/kg]                  
if is_GP | is_GPa                                   
 sva = alpha(P,T,S);                                
 sva(is_nan) = 0*is_nan;                            
end                                                 
if is_GPa                                           
 sva0 = sva-alpha(P,0*T,35+0*S);                    
 sva0(is_nan) = 0*is_nan;                           
end                                                 
                                                    
clear  T S 

%geopotential(Anomaly) [m^2/s^2] ( 1e4*P : pressure [dbar] --> Pa)
if is_GP 
 GP = -[zeros(1,d2); cumsum(0.5*(sva(1:d1-1,:) + sva(2:d1,:)) .* dP )] ;
end
if is_GPa
 GPa = -[zeros(1,d2); cumsum(0.5*(sva0(1:d1-1,:) + sva0(2:d1,:)) .* dP )] ;
end


clear sva sva0 dP jj 


if pref >= 0
 [hilf,k] = minnan(abs(P-pref)); 
       k  = ((1:d2)-1)*d1 + k;

  pref = round(min(P(k)));
   disp(['    GEOVEL2: Reference level set to ',num2str(pref),' dBar.'])

 [hilf,k] = minnan(abs(P-pref)) ;

elseif pref < 0

 [hilf,k]=maxnan(P); 
       k  = ((1:d2)-1)*d1 + k;

    Pini = P(k) + pref;
  Pini = ones(d1,1)*Pini(:)';
[hilf,k] = minnan(abs(P-Pini));

end



if ~isnan(pref)

 k  = ((1:d2)-1)*d1 + k;
if is_GP 
 GP = GP - ones(d1,1)*GP(k);
 GP(is_nan) = nan*is_nan;
end
if is_GPa
 GPa = GPa - ones(d1,1)*GPa(k);
 GPa(is_nan) = nan*is_nan;
end
end

clear is_nan k

if is_uv 
 if pref >= 0
  [U,V]=gp2uv2d(GPa,lon,lat);
 else
  [U,V]=gp2uv2d(GPa,lon,lat,pref,P);
 end
end
% is_uv 

vars = ['GP ';'GPa';'U  ';'V  '];

for ii = 1:Nout
 eval(['out' int2str(ii) ' = ' deblank(vars(Out(ii),:)) ';' ])
end

