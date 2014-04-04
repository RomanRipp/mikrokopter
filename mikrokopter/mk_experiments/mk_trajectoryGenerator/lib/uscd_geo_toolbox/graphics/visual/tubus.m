function [v,f] = tubus(geo,pos,hpr)

% TUBUS  Returns Faces and Vertices of a Tubus
%
% [V,F] = TUBUS(GEO,POS,HPR)
%
%--------------------------------------------------------
%
% GEO = [ LEN+INCL*i DO+DI*i RES ]  Geometry
%
%        LEN  = Length of Tube along X-Axis, centered to ZERO
%        INCL = Inclination of along Length, 
%                a nonzero Value for INCL creates a Cone
%
%        DO = outher Diameter, DI = Inner Diameter
%
% POS = [ PX PY PZ ]  Position of Center
%
% POS = [ PX PY PZ ] + DEV * i   Deviation of Center in CUBUS
%
% DEV = [ DevLEN DevRAD DIR ] DevLEN and DevRad are normalized to
%             half LEN or DO: DevLen * LEN/2; DevRad * DO/2
%
% HPR = [ Heading  Pitch  Roll ]  Pitch: positive down
%
%--------------------------------------------------------
%
%  V    XYZ-Coordinates of Vertices
%
%  F    VerticeIndex for Faces
%
% TUBUS without an Output plots the Tube,
%       default: regular Cube
%
%--------------------------------------------------------
%
% see also: CUBUS, CORPUS, GLIDER, HEPIRO, PATCH
%

Nin = nargin;

Nout = nargout;

if Nin < 1, geo = []; end
if Nin < 2, pos = []; end
if Nin < 3, hpr = []; end

%--------------------------------------------------------------
% Basic Check of Inputs

for v = { geo pos hpr }
    ok = isempty(v{1});
    if ~ok
        ok = ( isnumeric(v{1}) & ( prod(size(v{1})) <= 3 ) );
        if ok
           ok = all(isfinite(v{1}));
        end
    end
    if ~ok
        break
    end
end

if ~ok
    error('Inputs must be finite numeric with max. 3 Elements.')
end

%--------------------------------------------------------------

def = [ 5  1+0.8*i  72 ];

n = prod(size(geo));
m = size(def,2);

if n == 0
   geo = def;
elseif n < m
   if n < 2
      geo = cat(2,geo,real(geo(1))*def(2)/def(1));
   end
   geo(3) = def(3);
end

if isempty(pos)
   pos = [ 0  0  0 ];
elseif prod(size(pos)) < 3
   pos = cat(2,pos(:)',0,0);
   pos = pos(1:3);
end

if isempty(hpr)
   hpr = [ 0  0  0 ];
elseif prod(size(pos)) < 3
   hpr = cat(1,hpr(:),0,0);
   hpr = hpr(1:3);
end

%--------------------------------------------------------------

v = zeros(0,3);  % XYZ-Coordinates
f = zeros(0,4);  % Vertices

len = geo(1);
rad = geo(2);
res = geo(3);

ang = imag(len);
len = real(len);

rad = [ real(rad)  imag(rad) ] / 2;

rad = rad([1 2]+[1 -1]*(rad(1)<rad(2)));

if ( len == 0 ) | all( rad == 0 )
   return
end

rad = rad( 1 : (2-(rad(1)==rad(2))) );

rd  = rad + len * sin(ang*pi/180);   % Apply Inclination

rf  = rd ./ ( rad + ( rad == 0 ) );

%********************************************************

n = res;

p = pi/180 * 360/n * (0:(n-1))';

y = cos(p) * rad;
z = sin(p) * rad;

x = cat( 1 , -1*ones(size(y)) , 1*ones(size(y)) ) * len/2;

y = cat( 1 , y , y.*rf(ones(n,1),:) );
z = cat( 1 , z , z.*rf(ones(n,1),:) );

v = cat( 2 , x(:) , y(:) , z(:) );

c = ( 1 : n )';
c = cat( 2 , c , c+1 );
c = c - n * ( c > n );

% Base CylinderBarrel
f = cat( 2 , c , c(:,[2 1])+n );

if size(rad,2) > 1

   % Inner and Outher barrel
   if rad(1) == 0
      f = f + 2*n;
   elseif ~( rad(2) == 0 )
      f = cat( 1 , f , f + 2*n );
   end

   w = cat( 2 , c , c(:,[2 1])+2*n );  % Walls

   f = cat( 1 , f , w , w+n );

end

%--------------------------------------------------------------
% Position, Deviation and Rotation

o = ones(size(v,1),1);

dev = imag(pos) .* [ len/2  max(abs(rad)) 1 ];
pos = real(pos);

dev(3) = pi/180 * dev(3);
dev([2 3]) = [ cos(dev(3)) sin(dev(3)) ] * dev(2);

v = v + dev(o,:);

if ~all( hpr == 0 )  % Subfunction see below
    [v(:,1),v(:,2),v(:,3)] = hepirot(v(:,1),v(:,2),v(:,3),hpr);
end

v = v + pos(o,:);

%--------------------------------------------------------------
% Plot

if Nout == 0

   lim = cat( 1 , min(v,[],1) , max(v,[],1) );
   lim = permute( lim  , [ 2 1 ] );
   dlm = diff(lim,1,2);
   dlm = dlm + ( dlm == 0 );
   lim = lim + dlm/4 * [ -1  1 ];

   figure
   axes( 'xlim' , lim(1,:) , ...
         'ylim' , lim(2,:) , ...
         'zlim' , lim(3,:) , ...
         'dataaspectratio' , [1 1 1] , ...
         'nextplot' , 'add' , ...
         'view'     , [20 30] );
   patch('vertices',v,'faces',f,'edgecolor','none')

   lgt = light('style','infinite','color','w');
  
   lightangle(lgt,60,60)

   clear v f

end

%***************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [x,y,z] = hepirot(x,y,z,v);

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

Nin = nargin;

if Nin < 3
   error('X, Y and Z required.');
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

if Nin < 4
   v = [];
end

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

v = v - 360 * floor( v / 360 );   % [ 0 .. 360 )

if all( v == 0 )
   return
end

v = v * pi/180;

sv = sin(v);
cv = cos(v);

ex = [ 1  0  0 ];
ey = [ 0  1  0 ];
ez = [ 0  0  1 ];

T = cat( 1 , ex , ey , ez );

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

%------------------------------------------
% Transformation

xyz = cat(2,x(:),y(:),z(:));

xyz = xyz * T;

%------------------------------------------
% Reshape back

x = reshape(xyz(:,1),sz);
y = reshape(xyz(:,2),sz);
z = reshape(xyz(:,3),sz);

%***************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function p = comb(n,k)

if any( k == [ 0  1  n ] )

   if     k == 0
      p = zeros(1,0);
   elseif k == 1
      p = ( 1 : n )';
   else
      p = ( 1 : k );
   end   

   return

end

%------------------------------------------------------

m = nok(n,k);

p = zeros(m,k);

z = 0;

for ii = 1 : n-1
 
     m = nok(n-ii,k-1);
    jj = z + (1:m);

    p(jj,1) = ii;
    p(jj,2:k) = comb(n-ii,k-1) + ii;

     z = z + m;

end

%******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function n = nok(n,k)

n = prod(n-k+1:n) / prod(1:k);
