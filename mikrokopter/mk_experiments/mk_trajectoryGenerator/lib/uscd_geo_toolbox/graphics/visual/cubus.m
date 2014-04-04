function [v,f] = cubus(len,pos,hpr)

% CUBUS  Returns Faces and Vertices of a Cubus
%
% [V,F] = CUBUS(LEN,POS,HPR)
%
%--------------------------------------------------------
%
% LEN = [ LX LY LZ ]  Length in 3 Dimensions
%                     ZERO at Center, Range: -LEN/2 .. LEN/2
%
% POS = [ PX PY PZ ]  Position of Center
%
% POS = [ PX PY PZ ] + DEV * i   Deviation of Center in CUBUS
%
% DEV = [ DX DY DZ ] normalized to half LEN: POS + DEV*LEN/2
%
% HPR = [ Heading  Pitch  Roll ]  Pitch: positive down
%
%--------------------------------------------------------
%
%  V    XYZ-Coordinates of Vertices
%
%  F    VerticeIndex for Faces
%
% CUBUS without an Output plots the Cube,
%       default: regular Cube
%
%--------------------------------------------------------
%
% see also: TUBUS, CORPUS, GLIDER, HEPIRO, PATCH
%

Nin = nargin;

Nout = nargout;

if Nin < 1, len = []; end
if Nin < 2, pos = []; end
if Nin < 3, hpr = []; end

%--------------------------------------------------------------
% Basic Check of Inputs

for v = { len pos hpr }
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

if isempty(len)
   len = [ 1  1  1 ];
elseif prod(size(len)) < 3
   len = len([(1:end) end end]);
   len = len(1:3);
   len = len(:)';
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

n = ~( len == 0 );  % Non Zero Length
s = sum(n);
n = find(n);

if ( s == 0 )
   return
end

len = len/2;

% Take care on the Order! no alternating Signs

op = [  1  -1  -1 
        1   1  -1
       -1   1  -1 
       -1  -1  -1
       -1  -1   1 
       -1   1   1 
        1   1   1
        1  -1   1  ];

%   f = comb(size(op,1),4);
%   k = cat(3,op(f,1),op(f,2),op(f,3));
%   k = reshape(k,size(f,1),size(f,2),prod(size(k))/prod(size(f)));
%   k = find( any( all( ( diff(k,1,2) == 0 ) , 2 ) , 3 ) );
%   f = f(k,:)

f = [  1     2     3     4
       1     2     7     8
       1     4     5     8
       2     3     6     7
       3     4     5     6
       5     6     7     8  ];


m = 2^s;  % 2/4/8  Line/Plain/Cube

v = zeros( m , 3 );

o =  ones( m , 1 );

v(:,n) = op(1:m,1:s) .* len(o,n);

if ( s < 3 )   % Line / Patch 
   f = f( 1 , 1 : m );
end

%--------------------------------------------------------------
% Position, Deviation and Rotation

dev = imag(pos) .* ( len + ( len == 0 ) );
pos = real(pos);

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
   patch('vertices',v,'faces',f)

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
