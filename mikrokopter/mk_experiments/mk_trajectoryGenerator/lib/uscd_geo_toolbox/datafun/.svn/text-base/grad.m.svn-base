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
   x = 1;
end

x1 = ( prod(size(x)) == 1 );

if ~x1
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

if ~x1
    x = cat( 1 , 2*x(1,:)-x(2,:) , x , 2*x(n,:)-x(n-1,:) );
end

ind = ( 1 : n ) + 1;

switch ord

   case 1

     dy = ( y(ind+1,:) - y(ind-1,:) );

     if x1
        dx = 2 * x;
     else
        dx = ( x(ind+1,:) - x(ind-1,:) );
     end

     y = dy ./ ( dx + ( dx == 0 ) );
 
     is_inf = find( ( dx == 0 ) |  ( y > 1/acc ) );
 
     if ~isempty(is_inf)

         sy = sign( dy(is_inf) );
         sx = sign( dx(is_inf) + ( dx(is_inf) == 0 ) );

         y(is_inf) = sx.*sy * inf;

     end

   case 2

     if x1
         dx1 = x;
         dx2 = x;
         dx  = 2*x;
     else
         dx1 = ( x(ind+1,:) - x(ind+0,:) );
         dx2 = ( x(ind+0,:) - x(ind-1,:) );
         dx  = ( x(ind+1,:) - x(ind-1,:) );
     end

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

         y(is_inf) = sx.*sy * inf;

         is_inf = is_inf( find( ~( dx(is_inf) == 0 ) ) );
     
         if ~isempty(is_inf)     

             dy(is_inf) = dy1(is_inf) - dy2(is_inf);

              y(is_inf) = 4 * dy(is_inf) ./ ( dx(is_inf).^2 );
 
              jj = find( ( y(is_inf) > 1/acc ) );
 
             if ~isempty(jj)

                 sy = sign( dy(is_inf(jj)) );
                 sx = sign( dx(is_inf(jj)) + ( dx(is_inf(jj)) == 0 ) );

                 y(is_inf(jj)) = sx.*sy * inf;

             end

         end

     end

end

% Permute and reshape back

y = ipermute(reshape(y,si(perm)),perm);
