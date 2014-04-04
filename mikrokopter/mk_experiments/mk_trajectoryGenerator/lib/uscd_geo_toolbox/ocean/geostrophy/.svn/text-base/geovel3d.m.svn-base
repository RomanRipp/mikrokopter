function [out1,out2,out3,out4] = geovel3d(P,T,S,arg4,arg5,arg6,arg7)
%GEOVEL3D Calculation of geostrophic velocities.
%       (works 3-dimensional)
%       GEOVEL3d(P,T,S, ... )
%
%       Optional Inputs are:  Lon , Lat , Pref , '[<OutPut>]'
%          use them only in the order above.
%       For Velocities Lat and Lon are requiered.
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
%       Pressure         P      dBar     M x N x O or Length O
%       Temperature      T      deg C    M x N x O
%       Salinity         S      psu      M x N x O
%       Longitude        Lon    degrees  M x N x 1 or Length N 
%       Latitude         Lat    degrees  M x N x 1 or Length M
%       reference level  pref   dBar
%       GeoPotential     GP     m^2/s^2  dyn cm
%       GeoPotAnomaly    GPa    m^2/s^2
%       U-velocity       U      m/s
%       V-velocity       V      m/s
%
%       Velocities with absolute > 5 (m/s)  will be setted to NaN.
% 
%  If Pref negative ,it will be calculated 
%  from maximum, finite Element of each 3. Dimension in P.
%   
%  Pref = NaN , will give only GP or GPa  originally, 
%   no Difference to any Level.
%

%09.07.92 Christian Mertens, IFM Kiel
% changes to 3D  Christian Begler


Nout = nargout;

Nin = nargin;
if nargin < 3
 error('P,T,S required')
end

if length(size(T)) ~= 3                                                          
 error('Use GEOVEL2D for 2-dimensional Fields.')                                 
end


Out = [ 1 : Nout ]; 
pref = 0;
pNr = nan;

for ii = 4:Nin
 eval(['val = arg' int2str(ii) ';'])
  if isstr(val)
   Out = eval(['[' val ']']);
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
is_GP = [  any(Out==1) ];
is_GPa = [ is_uv |  any(Out==2)  ];

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
end

lon = [];
lat = [];
if ~isempty(ll)
  eval(['lon = arg' int2str(ll(1)) ';']);
  eval(['lat = arg' int2str(ll(2)) ';']); 
  eval([' clear arg' int2str(ll(1)) ' arg' int2str(ll(2)) ])
end
 
[d1,d2,d3] = size(T);
 

 p_si = size(P);                                                               
 if p_si(1:2) == [1 1]                                                         
  [hilf1,hilf2,P]=meshgrid(T(1,:,1),T(:,1,1),P);                               
    clear hilf1 hilf2
 end                                                                           
                                                                               
dim = 3;

is_nan = find( isnan(T) | isnan(S) | isnan(P) );


 dP = diff(P*1e4,1,dim);
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
 GP = zeros(d1,d2,d3);
 GP(:,:,2:d3) =  ...
   -cumsum(0.5*(sva(:,:,1:d3-1) + sva(:,:,2:d3)) .* dP , dim ) ;
end
if is_GPa
 GPa = zeros(d1,d2,d3);
 GPa(:,:,2:d3) =  ...
   -cumsum(0.5*(sva0(:,:,1:d3-1) + sva0(:,:,2:d3)) .* dP , dim ) ;
end


clear sva sva0 dP jj 

 d12 = d1*d2;
 d123 = d12 * d3;
 ind1 = reshape( (1:d1)'*ones(1,d2) , d12 , 1 );
 ind2 = reshape( ones(d1,1)*(1:d2) , d12 , 1 );

if pref >= 0
 [hilf,k] = min(abs(P-pref),[],dim); 
% save last1
        k = sub2ind([d1 d2 d3],ind1,ind2,reshape(k,d12,1));

 pref = round(min(P(k)));
  disp(['    GEOPOT3D: Reference level set to ',num2str(pref),' dBar.'])

 [hilf,k] = min(abs(P-pref),[],dim) ; clear hilf
% save last2
       k = sub2ind([d1 d2 d3],ind1,ind2,reshape(k,d12,1));

elseif pref < 0 
% pref < 0
 [Pmax,k]=max(P,[],dim);  
    Pmax = Pmax(:,:,ones(1,d3));
[hilf,k] = min(abs(P-(Pmax+pref)),[],dim); clear hilf
       k = sub2ind([d1 d2 d3],ind1,ind2,reshape(k,d12,1));

end

clear  Pini d12 ind1 ind2

% save last 


if ~isnan(pref)
  k = reshape( reshape(k*ones(1,d3),d123,1) , d1 , d2 , d3);
if is_GP
 GP = GP - GP(k);  
 GP(is_nan) = nan*is_nan;
else
 GP = [];
end
if is_GPa 
 GPa = GPa - GPa(k);
 GPa(is_nan) = nan*is_nan;
else
 GPa = [];
end
end
clear is_nan k




if is_uv 

 if pref >= 0
  [U,V] = gp2uv3d(GPa,lon,lat);
 else
  [U,V] = gp2uv3d(GPa,lon,lat,pref,P);
 end
else
 U = [];
 V = [];
end
% is_uv

vars = ['GP ';'GPa';'U  ';'V  '];

for ii = 1:Nout
 eval(['out' int2str(ii) ' = ' deblank(vars(Out(ii),:)) ';' ])
end

  


