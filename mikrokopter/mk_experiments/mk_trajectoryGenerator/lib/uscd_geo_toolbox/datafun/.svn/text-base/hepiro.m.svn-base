function [x,y,z,T] = hepiro(x,y,z,v);

% HEPIRO   Rotates XYZ-Coordinates by Heading, Pitch and Roll
%
% [ X , Y , Z ] = HEPIRO( X , Y , Z , HPR )
%
% HPR = [ Heading Pitch Roll ]
%
% Heading  Rotation around Z-Axis
% Pitch    Rotation around Y-Axis
% Roll     Rotation around X-Axis
%
% A single Input for HPR or a fourth Output 
% returns the RotationMatrice:
%
% T = HEPIRO( HPR )  or  [ ... , T ] = HEPIRO( ... , HPR ) 
%
%
%--------------------------------------------------------
% Example with GLIDER:
%
% [v,f,d,m] = glider;
%
% [x,y,z] = hepiro(v(:,1),v(:,2),v(:,3),[-40 30 20]);
%
% figure('ColorMap',m);
%
% axes( 'dataaspectratio' , [ 1  1  1 ] , ...
%       'view'            , [ -30  30 ] , ...
%       'nextplot'        , 'add'  );
%
% patch(     'vertices' , [ x y z ] , ...
%               'faces' , f         , ...
%     'facevertexcdata' , d         , ...
%        'cdatamapping' , 'direct'  , ...
%           'facecolor' , 'flat'    , ...
%           'edgecolor' , 'none'          )
%
% h = light('style','infinite','color','w');
%
% lightangle(h,60,60)
%

T = [];

Nin = nargin;

if Nin == 1
   v = x;
   x = []; y = []; z = [];
elseif Nin < 3
   error('X, Y and Z required.');
elseif Nin < 4
   v = [];
end

%**********************************************************************
% Check Inputs

msg = cell(0,1);

%---------------------------------------------------------------------
% XYZ

if ~( isnumeric(x) & isnumeric(y) & isnumeric(z) )
    msg = cat(1,msg,{'X, Y and Z must be numeric.'});
end

sz = size(z);

if ~isequal(size(x),size(y),sz)
    msg = cat(1,msg,{'Size of X, Y and Z must be agree.'});
end

%---------------------------------------------------------------------
% HPR

ok = isnumeric(v);
if ok
   p  = prod(size(v));
   ok = ( p <= 3 );
   if ok & ~isempty(v)
      v = cat( 1 , v(:) , zeros(3-p,1) );
   elseif ~ok
      msg = cat(1,msg,{'HPR must be numeric with max. 3 Elements.'});
   end
end

%---------------------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    msg = sprintf('Invalid Inputs.\n%s',msg);
    error(msg);
end

if isempty(v)
   return
end

%**********************************************************************

ex = [ 1  0  0 ];
ey = [ 0  1  0 ];
ez = [ 0  0  1 ];

T = cat( 1 , ex , ey , ez );

v = v - 360 * floor( v / 360 );   % [ 0 .. 360 )

if all( v == 0 )
   return
end

v = v * pi/180;

sv = sin(v);
cv = cos(v);

%*********************************************
% Zyclic Permutation: ZXY / YZX / XYZ
%------------------------------------------
% Head around Z
%    w = [  cv(1)  sv(1)  0
%          -sv(1)  cv(1)  0
%            0      0     1 ];
%    T = w * T;
%------------------------------------------
% Pitch around Y
%    w = [ cv(2)  0  -sv(2)
%           0     1    0   
%          sv(2)  0   cv(2) ];
%    T = w * T;
%------------------------------------------
% Roll around X
%    w = [  1     0     0   
%           0   cv(3)  sv(3)
%           0  -sv(3)  cv(3) ];
%    T = w * T;
%*********************************************

for ii = 1 : 3

    if ~( v(ii) == 0 )

        w = [  cv(ii)  sv(ii)  0
              -sv(ii)  cv(ii)  0
                0       0      1 ];

        jj = ( 1 : 3 ) - 1 + ii;     % Zyclic Index
        jj = jj - 3 * floor(jj/3);
        jj = jj + 3 * ( jj == 0 );
       
        T = w(jj,jj) * T;

    end

end

% Check for valid RotationMatrice
% d = inv(T) - T'; max(abs(d(:)))  % must be ZERO

if isempty(x)
   if Nin == 1
      x = T;
      T = [];
   end
   return
end

%------------------------------------------
% Transformation

xyz = cat(2,x(:),y(:),z(:));

xyz = xyz * T;

%------------------------------------------
% Reshape back

x = reshape(xyz(:,1),sz);
y = reshape(xyz(:,2),sz);
z = reshape(xyz(:,3),sz);


