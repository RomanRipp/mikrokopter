function  w = wekman(x,y,tx,ty,rho);

% WEKMAN  Calculates vertical Ekmanvelocity from Windstress
%
%  We = WEKMAN( Lon , Lat , Tx , Ty , [Rho] )
%
%  Lon = [ Ny by Nx ] or Nx-Element-Vector [deg]
%  Lat = [ Ny by Nx ] or Ny-Element-Vector [deg]
%
%  Tx  = [ Ny by Nx by Nt ]  X-Windstress [N/m²]
%  Ty  = [ Ny by Nx by Nt ]  Y-Windstress [N/m²]
%
% Note: [Tx,Ty] in same Direction as Windspeed !!!
%
%  Rho = [ Ny by Nx by Nt ] or [ Ny by Nx ] or Scalar [kg/m³]
%
% default: Rho = 1027
%
% Formula:  We = curl_z(T/f) / rho    [m/s]
%
%           curl_z(T/f) = d(Ty/f)/dx - d(Tx/f)/dy
%
%           f = 2 * (2*pi/86164) * sin(Lat)   CoriolisParameter
%

Nin = nargin;

if Nin < 4
   error('Not enough InputArguments');
end

if Nin < 5
   rho = [];
end

if isempty(rho)
   rho = 1027;
end

%********************************************************

r1 = ( prod(size(rho)) == 1 );

if r1
   [msg,x,y,tx,ty] = chk3d(x,y,tx,ty);
else
   [msg,x,y,tx,ty,rho] = chk3d(x,y,tx,ty,rho);
end

if ~isempty(msg)
    error(msg)
end

sx = size(tx);
sy = size(ty);
sr = size(rho);  r2 = ~r1 & ( ( prod(sr) == sr(1)*sr(2) ) );

ok = ( r1 | isequal(sx,sr) | ( r2 & isequal(sx([1 2]),sr) ) );

if ~( ok & isequal(size(tx),size(ty)) ) 
    error('Size of Tx and Ty (and Rho) must be agree');
end

%********************************************************

f = 2 * (2*pi/86164) * sin(y*pi/180); % CoriolisParameter

x = x * 60 * 1852 .* cos(y*pi/180);   % [m]
y = y * 60 * 1852;                    % [m]

o3 = ones(1,sx(3));

w = grad(ty./f(:,:,o3),x(:,:,o3),1,2) - ...
    grad(tx./f(:,:,o3),y(:,:,o3),1,1);

if r2
   w = w ./ rho(:,:,o3);
else
   w = w ./ rho;
end


%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function y = grad(y,x,ord,dim)

% GRAD   Gradient in 1. and 2. Order
%
%  G = GRAD( Y , X , Order , DIM )
%
%  Order == 1  ==>  G = dY/dX
%  Order == 2  ==>  G = d²Y/dX²
%
%  Use numeric estimation:
%
%  dY/dX = ( Y(n+1) - Y(n-1) ) / ( X(n+1) - X(n-1) )
%
%  d²Y/dX² = 2 * (  ( Y(n+1) - Y(n) ) / ( X(n+1) - X(n) )
%                 - ( Y(n) - Y(n-1) ) / ( X(n) - X(n-1) ) )
%              / ( X(n+1) - X(n-1) )
%
%
%   See also DIFF, GRADIENT
%

Nin = nargin;

if Nin < 1
   y = [];
   return
end

if isempty(y)
   return
end

if Nin < 2
   x = [];
end

if Nin < 3
  ord = [];
end

if Nin < 4 
  dim = [];
end

%-------------------------------------------------------------------------

si = size(y);

acc = 1e3*eps; % Accuracy

%-------------------------------------------------------------------------

if isempty(ord)
   ord = 1;
elseif ~any( ord(1) == [ 1  2 ] )
   error('Order must be any of 1 | 2 .')
end

ord = ord(1);


if isempty(dim)
  dim = min(find(~(si==1))); 
  if isempty(dim), dim = 1; end
end

dim = dim(1);

si = cat( 2 , si , ones(1,dim-size(si,2)) );
n  = si(dim);

if n < 2
   error('Length of Y in DIM must be >= 2.')
end

%-------------------------------------------------------------------------
% Permute and reshape so that DIM becomes the 1. dimension of a 2-D array

perm = cat( 2 , ( dim : max(length(size(x)),dim) ) , ( 1 : dim-1 ) );

y = reshape( permute(y,perm) , n , prod(si)/n );

sy = size(y);

if isempty(x)
   x = ones(sy); x = cumsum(x,1);
else
   sx = size(x);
   if ( sx(dim) == n ) & ( prod(sx) == n )
      x = x(:) * ones(1,sy(2));
   elseif ~isequal(sx,si)
      error('Size of Y and X must be agree.')
   else 
      x = reshape( permute(x,perm) , n , prod(si)/n );
   end
end 

%-------------------------------------------------------------------------

% Expand linear at Start and End

y = cat( 1 , 2*y(1,:)-y(2,:) , y , 2*y(n,:)-y(n-1,:) );
x = cat( 1 , 2*x(1,:)-x(2,:) , x , 2*x(n,:)-x(n-1,:) );

ind = ( 1 : n ) + 1;

switch ord

   case 1

     dy = ( y(ind+1,:) - y(ind-1,:) );
     dx = ( x(ind+1,:) - x(ind-1,:) );

     
     y = dy ./ ( dx + ( dx == 0 ) );
 
     is_inf = find( ( dx == 0 ) |  ( y > 1/acc ) );
 
     if ~isempty(is_inf)

       sy = sign( dy(is_inf) );
       sx = sign( dx(is_inf) + ( dx(is_inf) == 0 ) );

       y(is_inf) = sx*sy * inf;

     end

   case 2

     dx1 = ( x(ind+1,:) - x(ind+0,:) );
     dx2 = ( x(ind+0,:) - x(ind-1,:) );

     dx  = ( x(ind+1,:) - x(ind-1,:) );

     dy1 = ( y(ind+1,:) - y(ind+0,:) );
     dy2 = ( y(ind+0,:) - y(ind-1,:) );

     d1 = dy1 ./ ( dx1 + ( dx1 == 0 ) );
     d2 = dy2 ./ ( dx2 + ( dx2 == 0 ) );

     i1 = find( ( dx1 == 0 ) );
     i2 = find( ( dx2 == 0 ) );
 
     d1(i1) = sign(dy1(i1)) * inf;
     d2(i2) = sign(dy2(i2)) * inf;

     dy = ( d1 - d2 );

     y  = 2 * dy ./ ( dx + ( dx == 0 ) );

     is_inf = find( ( dx == 0 ) |  ( y > 1/acc ) );
 
     if ~isempty(is_inf)

       sy = sign( dy(is_inf) );
       sx = sign( dx(is_inf) + ( dx(is_inf) == 0 ) );

       y(is_inf) = sx*sy * inf;

       if ~isempty(is_inf)     

          is_inf = is_inf( find( ~( dx(is_inf) == 0 ) ) );
     
          if ~isempty(is_inf)     

             dy(is_inf) = dy1(is_inf) - dy2(is_inf);

             y(is_inf) = 4 * dy(is_inf) ./ ( dx(is_inf).^2 );
 
             jj = find( ( y(is_inf) > 1/acc ) );
 
             if ~isempty(jj)

               sy = sign( dy(is_inf(jj)) );
               sx = sign( dx(is_inf(jj)) + ( dx(is_inf(jj)) == 0 ) );

               y(is_inf(jj)) = sx*sy * inf;

             end

          end

       end

     end

end

% Permute and reshape back

y = ipermute(reshape(y,si(perm)),perm);


%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,x,y,varargout] = chk3d(x,y,varargin)

% CHK3D Check arguments to 3-D data routines.
%
%   [MSG,X,Y,Z1,Z2,...] = CHK3D(X,Y,Z1,Z2,...)
%

msg = '';

Nin  = nargin - 2;
Nout = nargout - 3;

if Nin < 0
   error('Inputs X ynd Y are missing.');
end


nn = 0;

if Nout > 0

   varargout = cell(1,Nout);

   if Nin > 0
      nn = min(Nin,Nout);
      varargout(1:nn) = varargin(1:nn);
   end

end


sx = size(x);  px = prod(sx);  vx = ( px == max(sx) );
sy = size(y);  py = prod(sy);  vy = ( py == max(sy) );

%****************************************************************
% X & Y only

if Nin == 0

   if isequal(sx,sy)
      return
   end

   if     vx & vy
      x = ones(py,1) * x(:)';
      y = y(:) * ones(1,px);
   elseif vx & ( px == sy(2) )
      x = ones(sy(1),1) * x(:)';
   elseif vy & ( py == sx(1) )
      y = y(:) * ones(1,sx(2));
   else
      msg = 'Size of X and Y must be agree.';
   end
   
   return

end

%****************************************************************
% Get Size of Matrices

sz = NaN * zeros(Nin,2);

for ii = 1 : Nin
 
    s = size(varargin{ii});

    if size(s,2) <= 3
       sz(ii,[1 2]) = s([1 2]);
    end

end

if any(isnan(sz(:,1)))
   msg = 'Following Inputs must be Matrices with max. 3 Dimensions.';
   return
end

ds = sz - sz(ones(Nin,1),:);

if ~all( ds(:) == 0 )
    pz = prod(sz,2);
    vz = ( pz == max(sz,[],2) );
    if ~all( vz & ( pz == pz(1) ) )
        msg = 'Following Inputs must have the same Size in X and Y.';
        return
    end
    for ii = 1 : nn
          varargout{ii} = varargout{ii}(:);
          if ( sz(1,2) == pz(1) )
             varargout{ii} = varargout{ii}';
          end
    end 
end
 
%****************************************************************
% Check Size of Matrices with X and Y
  
sz = sz(1,[1 2]);  pz = prod(sz);  vz = ( pz == max(sz) );

if isequal(sx,sy,sz)
   return
end

%---------------------------------------------------------------
% Check for Vectors

if vz

   %---------------------------------------
   % X and Y not Vectors
   if ~( vx & vy ) 
      msg = 'When Z is a vector, X and Y must also be vectors.';
      return
   end

   %---------------------------------------
   % Vectors of same Lenght
   if isequal(px,py,pz)
      x = x(:);
      y = y(:);
      for ii = 1 : nn
          varargout{ii} = varargout{ii}(:);
      end
      if ( sz(2) == pz )
         x = x';
         y = y';
         for ii = 1 : nn
             varargout{ii} = varargout{ii}';
         end
      end
      return
   end

   %---------------------------------------
   % X and Y build Z
   if isequal([px py],sz)
      if     ~( px == pz )
          x = x * ones(sz);
      elseif ~( py == pz )
          y = y * ones(sz);
      end
      return
   end

   msg = 'Vectors X and Y must match Size of Vector Z.';
   return

end

%---------------------------------------------------------------
% Check MatriceSize

if     vx & vy
   x = ones(py,1) * x(:)';
   y = y(:) * ones(1,px);
elseif vx & ( px == sy(2) )
   x = ones(sy(1),1) * x(:)';
elseif vy & ( py == sx(1) )
   y = y(:) * ones(1,sx(2));
elseif ~isequal(sx,sy)
   msg = 'Size of X and Y must be agree.';
   return
end

if ~isequal(size(x),size(y),sz)
    msg = 'Size of X and Y must match Size of Z.';
    return
end
