function [z,nx,ny,nz] = surfarea(x,y,z)

% SURFAREA  Surface Area
%
% C = SURFAREA( X , Y , Z )
%
% returns the corresponding area for each point C
%  of the 3-D gridded surface with the components (X,Y,Z).
%
% The area C is the magnitude of the cross-product of
%  the surface normals in each point of the grid (X,Y,Z).
%
%
% C = SURFAREA(Z) with Z = [ M by N ] 
%  assumes X = ( 1 : N ) and Y = ( 1 : M ) 
%
%
% [C,NX,NY,NZ] = SURFAREA( X , Y , Z )
%
% returns the surface normal components (NX,NY,NZ)
%  for the surface Z  by SURFNORM.
% The normal is normalized to length 1.
%
%
% see also: SURFNORM, POLYAREA
%

%***************************************************
% Inquire Inputs

Nin = nargin;

if ~any( Nin == [ 1  3 ] )
    error('Required Inputs are Z  or  X, Y and Z.');
end

%----------------------------------------------------
% Check of Z

if Nin == 1
   z = x;
end

s = size(z);
m = s(1);
n = s(2);

if ~( ( prod(s) == m*n ) & all( s >= 3 ) )
    error('Z must be 2-D and at least 3-by-3.')
end


if Nin == 1
   x = ( 1 : n );
   y = ( 1 : m );
end

%----------------------------------------------------
% Check Size of X and Y

sx = size(x); px = prod(sx); vx = ( px == max(sx) );
sy = size(y); py = prod(sy); vy = ( py == max(sy) );

if ~isequal(sx,s)
    if vx & ( px == n )
       x = x(:)';
       x = x(ones(m,1),:);
      sx = s;
    end
end

if ~isequal(sy,s)
    if vy & ( py == m )
       y = y(:);
       y = y(:,ones(1,n));
      sy = s;
    end
end

if ~isequal(sx,sy,s)
    error('MatrixDimensions of XYZ must be agree.');
end

%***************************************************
% Surface Normals like in SURFNORM

% Expand x,y,z so interpolation is valid at the boundaries.

x = [3*x(1,:)-3*x(2,:)+x(3,:);x;3*x(m,:)-3*x(m-1,:)+x(m-2,:)];
x = [3*x(:,1)-3*x(:,2)+x(:,3),x,3*x(:,n)-3*x(:,n-1)+x(:,n-2)];

y = [3*y(1,:)-3*y(2,:)+y(3,:);y;3*y(m,:)-3*y(m-1,:)+y(m-2,:)];
y = [3*y(:,1)-3*y(:,2)+y(:,3),y,3*y(:,n)-3*y(:,n-1)+y(:,n-2)];

z = [3*z(1,:)-3*z(2,:)+z(3,:);z;3*z(m,:)-3*z(m-1,:)+z(m-2,:)];
z = [3*z(:,1)-3*z(:,2)+z(:,3),z,3*z(:,n)-3*z(:,n-1)+z(:,n-2)];


rows = ( 2 : m+1 ); 
cols = ( 2 : n+1 );

stencil1 =  [  1    0   -1 ] / 2;
stencil2 =  [ -1 ;  0 ;  1 ] / 2;

xx = filter2(stencil1,x); xx = xx(rows,cols);
yy = filter2(stencil1,y); yy = yy(rows,cols);
zz = filter2(stencil1,z); zz = zz(rows,cols);

x = filter2(stencil2,x); x = x(rows,cols);
y = filter2(stencil2,y); y = y(rows,cols);
z = filter2(stencil2,z); z = z(rows,cols);

% Perform cross product to get normals

nx = -( yy.*z - zz.*y );
ny = -( zz.*x - xx.*z );
nz = -( xx.*y - yy.*x );

z  = sqrt( nx.*nx + ny.*ny + nz.*nz );

if nargout <= 1
   return
end

nx = nx ./ z;
ny = ny ./ z;
nz = nz ./ z;
