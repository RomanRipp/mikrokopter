function [uv,phi,mn,pc,e,l] = paxes(uv)

% PAXES  Principle Axes of 2-Component TimeSeries
%
% [XY,PHI,MN] = PAXES(UV)
%
% UV =  [ U V ] Components (2 Columns!)
%
% XY  = Principle Axes Components
% PHI = Orientation of Principle Axes
% MN  = MeanValues of UV-Components, removed before SVD
%
% [XY,PHI,MN,PC,E,L] = PAXES(UV)
%
% E  BasisVectors for XY
% L  EigenValues
%
% Use SVD-Analysis:
%
% [V,L,F] = SVD(UV-MN)
%
% UV-MN = V*L*F'
%
% XY = V*L; E = F';
%
% UV-MN = XY * E ; E = [ EX ; EY ]
%

phi = [];
mn  = [];
pc  = [];
e   = [];

if isempty(uv)
   return
end

if ~( ( ndims(uv) == 2 ) & ( size(uv,2) == 2 ) )
    error('Two Columns required.');
end

%----------------------------------------------------------

% Remove Mean

o = ones(size(uv,1),1);

mn = mean(uv,1);
uv = uv - mn(o,:); 

%----------------------------------------------------------
% SVD-Analysis
%
% [V,L,E]=svd(uv,0); 
%
% uv = V*L*E'
%
%      V  EigenVectors  sum(V.^2,1) == 1
% diag(L) EigenValues   for EigenVectors
%      E  BasisVectors for EigenVectors
%      E = [ ex  ey ], [det(E) = -1, E =] E' = inv(E)
%      New Coordinate-Axes for V=[Vx Vy] in BasisSystem 
%       

[uv,l,e]=svd(uv,0); 

if ( det(e) < 0 )  % Flip if negative Determinant!!!
    e(:,2) =  -e(:,2);
   uv(:,2) = -uv(:,2);
end

e = e';  % !!! [ ex ; ey ] !!!

pc = diag(l.^2); 
pc = (pc/sum(pc))*100;

%
%  V * L = V .* (ones(size(V,1),1)*diag(L)')
%

l = diag(l);
l = l(:)';

uv = uv .* l(o,:);

phi = atan2( e(1,2) , e(1,1) ) * 180/pi; % !!! ATAN2 !!!

