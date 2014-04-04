function z = plateau(x,y,r,c,m,p);

% PLATEAU  Create Plateau-Surface around Center
%
% Z = PLATEAU( X , Y , Radius , Center , Mode , Potenz )
%
% Creates a PLATEAU around Center defined by an (elliptic)
%  Area (Radius). 
%
% The Values of Z inside the Region of Radius around Center ranges
%  from ONE to ZERO.
%
% Radius  = [ RadiusX  RadiusY  [Rotation+Offset*i] ] 
%
%           Rotation of Radius anticlockwise, default: 0
%
%           The Distance of [X,Y] to Center is scaled 
%            to Radius, that the Value at the Elliptic Border is ONE.
%
%           Offset set the RadiusScaled Distance around Center to ZERO,
%             Z returns 1 in that Region.
%           Offset = [ 0 .. 1 ), default: 0
%
% Center  = [ CenterX  CenterY ], default: [ mean(X)  mean(Y) ]
%
% Mode    = Mode of window: 'triangle' | {'cosine'} |  'gauss'  |  'exp' 
%
%           gauss: exp( -(pi*x)^2 / 2 );  exp: exp(-x);
%
%          In case of Mode "gauss" or "exp" the Values 
%             nears asymtotic to ZERO!
%
% Mode    = positive Numeric for Power of Window: ( -1 .. 1 )
%             Mode == 0    ==>  Impuls
%             Mode == 1    ==>  Triangle
%             Mode == Inf  ==>  Linear (constant)
%
% Potenz  = Potenz for RadiusScaled Distance 
%           Potenz = [ 0 .. 1 ], default:  1
%
%-------------------------------------------------------------------------
%
% For Demonstration for Shape of Plateau see PLATDEMO, type: >> platdemo 
%
% see also: WINDOW, WININT, INELLIP
%
%-------------------------------------------------------------------------
% Examples:
%
%  x = linspace(-1,1,201); y = x';
%  z = plateau(x,y,[0.8 0.8 0.2*i],[],'cos');
%
%  figure('colormap',flipud(bone(128)));
%
%  axes( 'xlim' , x([1 end])  , ...
%        'ylim' , y([1 end])' , ...
%        'clim' , [-0.5 1.5 ] , ...
%        'view' , [ 20  30 ] , ...
%'dataaspectratio' , [ 1  1  1 ], ...
%       'nextplot' , 'add' );
%
%  surface( 'xdata' , x , ...
%           'ydata' , y , ...
%           'zdata' , z , ...
%           'cdata' , z , ...
%    'cdatamapping' , 'scaled' , ...
%    'edgecolor'    , 'none'   , ...
%    'facecolor'    , 'interp' , ...
%    'facelighting' , 'gouraud' );
%
%  lightang( gca  , [ 90  60 ] ,  ...
%          'color', [1 1 1]   , ... 
%          'style','infinite' , ...
%          'visible' , 'on' );
%
%--------------------------------------------------------------
%% Overlay Image with Shadow from lower left Corner
%
%  c = [ 10 10 80]/255;         % A Blue Shadow
%  c = permute(c,[1 3 2]);
%
%  load earth  %  X  map
%
%  s = size(X);
%  n = size(map,1);
%  X  = map(cat(3,X,X+n,X+2*n));  % RGB
%
%  z = plateau( (1:s(2)) , (1:s(1))' , 2/3*min(s) , s , 'exp' );
%
%  X = ( 1 - z(:,:,ones(1,3)) ) .* X + ...
%            z(:,:,ones(1,3))   .* c(ones(1,s(1)),ones(1,s(2)),:);
%
%  figure, image(X), axis image
%

Nin = nargin;

if Nin < 3
   error('Not enough InputArguments.');
end

if Nin < 4
   c = [];
end

if Nin < 5
   m = 'c';
end

if Nin < 6
   p = 1;
end

msg = cell(0,1);

%---------------------------------------------------------------------
% X and Y

sx = size(x);
sy = size(y);

if ~isequal(sx,sy)

    px = prod(sx);  vx = ( px == max(sx) );
    py = prod(sy);  vy = ( py == max(sy) );

    if     vx & vy
      x = ones(py,1) * x(:)';
      y = y(:) * ones(1,px);
    elseif vx & ( px == sy(2) )
      x = ones(sy(1),1) * x(:)';
    elseif vy & ( py == sx(1) )
      y = y(:) * ones(1,sx(2));
    else   ~isequal(sx,sy)
      msg = cat(1,msg,{'Size of X and Y must be agree.'});
    end

end

%---------------------------------------------------------------------
% Radius

pr = prod(size(r));

if ~( ( 0 < pr ) & ( pr <= 3 ) & isnumeric(r) )

    msg = cat(1,msg,{'Radius must be a numeric with 1 .. 3 Elements.'});

else

   w = 0;

   if     pr == 1
      r = r([1 1]);
   elseif pr == 3
      w = r(3);
      r = r(1:2);
   end

   o = imag(w);   % Offset
   w = real(w);

   if ~all( ( r >= 0 ) & isfinite(r) )
      msg = cat(1,msg,{'Radius must be positive.'});
   end
 
   if ~( ( 0 <= o ) & ( o < 1 ) )
      msg = cat(1,msg,{'Offset must be a positive numeric in Interval [ 0  1 ).'});
   end
       
end

%---------------------------------------------------------------------
% Center

pc = prod(size(c));

if     pc == 0 
   if isempty(msg) & ~isempty(x)
      c = [ mean(x(:)) mean(y(:)) ];
   end
elseif ( pc == 1 ) &  isnumeric(c)
   c = c([1 1]);
elseif ( pc >  2 ) | ~isnumeric(c)
    msg = cat(1,msg,{'Center be a numeric with max. 2 Elements.'});
end

%---------------------------------------------------------------------
% Mode

e = 1;

ok = ( ischar(m) & ( prod(size(m)) == size(m,2) ) & ~isempty(m) );

if ok
   m = lower(m(1));
else
   ok = ( isnumeric(m)  &  ( prod(size(m)) == 1 ) );
   if ok
      ok = ( ~isnan(m) & ( m >= 0 ) );
   end
   if ok
      e = m;
      m = 'p';
   end
end

if ~ok
    msg = cat( 1 , msg , {'Mode must be a String or positive Numeric.'} );
end

%---------------------------------------------------------------------
% Potenz

ok = ( isnumeric(p)  &  ( prod(size(p)) == 1 ) );
if ok
   ok = ( ~isnan(p) & ( p >= 0 ) );
end
    
if ~ok
    msg = cat( 1 , msg , {'Potenz must be a  positive Numeric.'} );
end

%---------------------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Invalid Inputs.\n%s',msg));
end

if isempty(x) | any( r == 0 ) | ( p == 0 )
   z = zeros(sx);
   return
end

%*********************************************************************
% Rotation

x = x - c(1);
y = y - c(2);

if ~( mod(w,180) == 0 )

    if mod(w,90) == 0

       r = r([2 1]);

    elseif ~( r(1) == r(2) )

       w = w * pi/180;

       s = sin(w);
       c = cos(w);
 
       z = x;
       x =  c * z + s * y;
       y = -s * z + c * y;

    end

end

%*********************************************************************
% Distance, weight by Radius and Potenz

z = sqrt( (x/r(1)).^2 + (y/r(2)).^2 );

if ~( p == 1 )

     z = z .^ p;

end

%*********************************************************************

nn = find( isnan(z) );

if ~( o == 0 )

    z = ( z - o ) / ( 1 - o );
    z = max( z , 0 );

end

if ~any( m == 'ge' )
    z = min( z , 1 );
end

%-----------------------------------------------------

switch lower(m(1))

  %-------------------------------------------------
  % Cosine

  case 'c'

    z = ( 1 + cos( pi * z ) ) / 2;

  %-------------------------------------------------
  % Gauss

  case 'g'

    z = exp( (-1) * (pi*z).^2 / 2 );

  %-------------------------------------------------
  % Gauss

  case 'e'

    z = exp( -z );

  %-------------------------------------------------
  % Triangle | Potenz

  otherwise

    if isinf(e)

       z = ( abs(z) < 1 );

    else
 
       if e == 0
          z = double( z == 0 );
       elseif e == 1
          z = 1 - abs(z);
       else    
          z = 1 - abs(z) .^ e;
       end

    end

end

if ~isempty(nn)
    z(nn) = NaN;
end
