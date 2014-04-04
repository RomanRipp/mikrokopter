function zi = interp21(x,y,z,xi,yi,mode);

% INTERP21  interpolates a Single Point into a 2-D function
%
%  ZI = INTERP21( X , Y , Z , XI , YI );
%
%   X = [  1 by Nx ]   GridVector for X 
%   Y = [ Ny by  1 ]   GridVector for Y
%   Z = [ Ny by Nx ]   MatriceFunction
%
%  XI = [ 1 by 1 ]  Target X
%  YI = [ 1 by 1 ]  Target Y
% 
%  ZI = [ 1 by 1 ]  
%
% If XI or YI are out of range, ZI will set to NaN.
%
% The Interpolation Method is an invers weight.
%
%  ZI = INTERP21( X , Y , Z , XI , YI , 'geo' );
%  
%  Interpretes X as Longitude, Y as Latitude to determine
%   invers weighted Distance
%


zi = [];

Msg = '';

nl = char(10);


if nargin < 5
 error('Number of Inputs must be 5.' )
end


if nargin < 6
  mode = 'none';
end


x = x(:);
y = y(:);



if ~isequal( size(z) , [ size(y,1)  size(x,1) ] )
  Msg = 'Size of Z must be  [ length(Y)  by  length(X) ].';
end

if ~all( ( diff(x) < 0 )  |  ( diff(x) > 0 ) )  |  ...
   ~all( ( diff(y) < 0 )  |  ( diff(y) > 0 ) ) 
  Msg = [ Msg nl(1:(end*(~isempty(Msg))))   ...
           ' X and Y must be monotonic. ' ];
end

if ( prod(size(xi)) ~= 1  )  |  ...
   ( prod(size(yi)) ~= 1  )
  Msg = [ Msg nl(1:(end*(~isempty(Msg))))   ...
          'XI, YI must define a single Target.' ];
end

mode = mode(:)';

if ~ischar(mode) 
  Msg = [ Msg nl(1:(end*(~isempty(Msg))))   ...
            ' Additional Inputs  could be String ''geo''.' ];  
end


if ~isempty(Msg)
  error([ ' INTERP21: ' Msg ])
end

zi = NaN;

if strcmp(mode,'geo')

  y  = deg2merc( y  , 1 );
  yi = deg2merc( yi , 1 );

end


 fx = 1 - 2 * ( x(1) > x(end) );  % Monotonic
 fy = 1 - 2 * ( y(1) > y(end) );

 ix1 = sum( fx*x <= fx*xi );
 iy1 = sum( fy*y <= fy*yi );

 ix2 = size(x,1) - sum( fx*x >= fx*xi ) + 1;
 iy2 = size(y,1) - sum( fy*y >= fy*yi ) + 1;


if ( ix1 >= 1 ) & ( iy1 >= 1 ) & ...
   ( ix2 <= size(x,1)  ) &  ( iy2 <= size(y,1)  )      
  
  if ( ix1 == ix2 ) & ( iy1 == iy2 )

    zi = z( iy1 , ix1 );

  else
       
    ix = [ ix1 ; ix2 ; ix2 ; ix1 ];
    iy = [ iy1 ; iy1 ; iy2 ; iy2 ];

    x = x(ix);
    y = y(iy);
    z = z( (ix-1)*size(z,1) + iy );


    % Normalized to Grid !!!

    gx = ( x(1) - x(2) );
    gy = ( y(2) - y(3) );

    dx = abs( ( x - xi ) / ( gx + ( gx == 0 ) ) );  
    dy = abs( ( y - yi ) / ( gy + ( gy == 0 ) ) );

    d  = sqrt( dx.^2 + dy.^2 );
    d  = ( sqrt(2) - d ) .* ( 1 - dx ) .* ( 1 - dy );

    zi = sum( z.*d ) / sum(d);

  end
          
end
