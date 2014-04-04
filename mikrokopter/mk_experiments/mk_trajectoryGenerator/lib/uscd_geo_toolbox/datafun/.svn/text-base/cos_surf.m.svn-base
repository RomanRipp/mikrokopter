function  z = cos_surf(x0,y0,z0,x,y,e);

% COS_SURF  2-D Interpolation via Cosine
%
%  Z = COS_SURF(X0,Y0,Z0,X,Y,CosineExponent);
%
%   the Vectors X0 and Y0 defines the TargetsPoints ,
%    the Values at this Points are defined by Matrice Z0.
%
%   X and Y are Matrices or Vectors by of same size, 
%     on this Coordinates the Data of Z0 will interpolated.   
%
%  CosineExponent is the  Exponent for Cosine, [ CosExpX  CosExpY ]
%
%

Nin = nargin;

if Nin < 5
   error('Not enough InputArguments.');
end

if nargin < 6
   e = [] ;
end

msg = cell(0,1);

%---------------------------------------------------------------------

if isempty(e)
   e = [ 1  1 ];
else
   p = prod(size(e));
   if p > 2
      msg = cat(1,msg,{'CosineExponent must have max. 2 Elements.'});
   elseif p == 1
      e = e([1 1]);
   end
end

%---------------------------------------------------------------------

x0 = x0(:)';
y0 = y0(:);

sx = size(x0,2);
sy = size(y0,1);

msg = cell(0,1);

if ~isequal([ sy  sx ] , size(z0) ) 
    msg = cat( 1 , msg , {'Length of Vectors X0, Y0 must correpond with ' ...
                          'Size of Z0.'} );
end

%---------------------------------------------------------------------

sx = size(x);  px = prod(sx);  vx = ( px == max(sx) );
sy = size(y);  py = prod(sy);  vy = ( py == max(sy) );

if     vx & vy
      x = ones(py,1) * x(:)';
      y = y(:) * ones(1,px);
elseif vx & ( px == sy(2) )
      x = ones(sy(1),1) * x(:)';
elseif vy & ( py == sx(1) )
      y = y(:) * ones(1,sx(2));
elseif ~isequal(sx,sy)
      msg = cat(1,msg,{'Size of X and Y must be agree.'});
end

%---------------------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Invalid Inputs.\n%s',msg));
end


%---------------------------------------------------
% Sort Targets

  [x0,s_x] = sort(x0);
  [y0,s_y] = sort(y0);

   z0      = z0(s_y,s_x);


  s1 = size(x,1);
  s2 = size(x,2);


  x = x(:);
  y = y(:);


% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%  Length > 10000  ==>  Use parts of Length == 10000  

z = zeros(s1*s2,1);

m = s1*s2 + ( 1e4 - s1*s2 ) * ( s1*s2 > 1*1e4 );

n = ceil(s1*s2/m);

for ii = 1 : n

 ende = m*ii + (s1*s2-m*ii)*( s1*s2 < m*ii );

 ind = ( m*(ii-1)+1 :  ende );

 z(ind) = get_surf(x0(:),y0(:),z0,x(ind),y(ind),e);

end

  z = reshape(z,s1,s2);



%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  z = get_surf(x0,y0,z0,x,y,e)


  sx = size(x0,1);
  sy = size(y0,1);


%---------------------------------------------------
% Set Borders

  x = x + ( x0( 1) - x ) .* ( x < x0( 1) );
  x = x + ( x0(sx) - x ) .* ( x > x0(sx) );

  y = y + ( y0( 1) - y ) .* ( y < y0( 1) );
  y = y + ( y0(sy) - y ) .* ( y > y0(sy) );


%--------------------------------------------------
% Define Cosine

 % ZeroValues
   zx = ( z0(:,2:sx)+z0(:,1:sx-1) )/2;
   zy = ( z0(2:sy,:)+z0(1:sy-1,:) )/2;

   zx = zx( : , [ (1:sx-1) sx-1 ]); 
   zy = zy([ (1:sy-1) sy-1 ] , : ); 

 % Amplitudes
   ax = -( z0(:,2:sx)-z0(:,1:sx-1) )/2;  
   ay = -( z0(2:sy,:)-z0(1:sy-1,:) )/2; 

   ax = ax( : , [ (1:sx-1) sx-1 ]); 
   ay = ay([ (1:sy-1) sy-1 ] , : ); 

 % Periods
   lx = 2 * ( x0(2:sx)-x0(1:sx-1) ); 
   ly = 2 * ( y0(2:sy)-y0(1:sy-1) );

   lx = lx([ (1:sx-1) sx-1 ]); lx = lx(:);
   ly = ly([ (1:sy-1) sy-1 ]); ly = ly(:);
 
 % Phase

   px = x0([ (1:sx-1) sx-1 ]); px = px(:);
   py = y0([ (1:sy-1) sy-1 ]); py = py(:);
 


%---------------------------------------------------
% Interpolate x and y to Targets

  xi = interp1( x0 , ( 1 : sx ) , x(:) );  % Index to lx
  yi = interp1( y0 , ( 1 : sy ) , y(:) );  % Index to ly


  % Surrounding PointIndize of x0, y0
  xi = [ floor(xi)  ceil(xi) ];
  yi = [ floor(yi)  ceil(yi) ];

  i2 = ones(1,2);

  dx = x(:,i2)-x0(xi);
  dy = y(:,i2)-y0(yi);

  dx = dx  + (sum(dx,2)==0)*i2;
  dx = dx ./ (sum(abs(dx),2)*i2);

  dy = dy  + (sum(dy,2)==0)*i2;
  dy = dy ./ (sum(abs(dy),2)*i2);


  % Weigth Zero and Amplitude with Inverse of y , x

%  zx = sum( ( zx( yi + (xi(:,i2)-1) * sy ) .* (1-abs(dy)) ) , 2 );
%  ax = sum( ( ax( yi + (xi(:,i2)-1) * sy ) .* (1-abs(dy)) ) , 2 );

%  zy = sum( ( zy( yi(:,i2) + (xi-1) * sy ) .* (1-abs(dx)) ) , 2 );
%  ay = sum( ( ay( yi(:,i2) + (xi-1) * sy ) .* (1-abs(dx)) ) , 2 );

   cx = cos( 2*pi*( (x-px(xi(:,1)))./lx(xi(:,1))));
   cy = cos( 2*pi*( (y-py(yi(:,1)))./ly(yi(:,1))));

  
  z = sum( ( zx( yi + (xi(:,i2)-1) * sy ) .* (1-abs(dy)) ) , 2 )  + ...
      sum( ( ax( yi + (xi(:,i2)-1) * sy ) .* (1-abs(dy)) ) , 2 ) .* ...   
       sign(cx).*(abs(cx).^e(1))                                 + ...
      sum( ( zy( yi(:,i2) + (xi-1) * sy ) .* (1-abs(dx)) ) , 2 )  + ...
      sum( ( ay( yi(:,i2) + (xi-1) * sy ) .* (1-abs(dx)) ) , 2 ) .* ...
       sign(cy).*(abs(cy).^e(2));

%       (1-2*(cx<0)).*(((1-2*(cx<0)).*cx).^e)                                 + ...
%       (1-2*(cy<0)).*(((1-2*(cy<0)).*cy).^e);

  % Mean Value of x  y 

  z = z/2;