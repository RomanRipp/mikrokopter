function  [y,varargout] = cos_line(xy,x,varargin);

% COS_LINE  1-D Interpolation via Cosine
%
%  Y     = COS_SURF(XY,X,Orientation,Option);
%  [Y,X] = COS_SURF(XY,N,Orientation,Option);
%
%  XY == 2-Column|Row Matrice, defines the TargetPoints in X0 and Y0
%
%   X == Vector or Matrice to interpolate a CosineLine trough the TargetPoints,
%         corresponding with first (second in Orientation'X') Row|Column of XY
%
%   N == Single Number with Value N uses an N-Point-Line in Intervall:
%  
%   X = [ min(X0) , Y0(min(X0)) ] .. [ max(X0) , Y0(max(X0)) ];  Orientation == 'X'
%   X = [ min(Y0) , X0(min(Y0)) ] .. [ max(Y0) , X0(max(Y0)) ];  Orientation == 'Y'
%
%
%  Orientation == 'X' | 'Y'; Orientation for V, default: 'X'
%
%  Option == Bearing + i*CosineExponent
%
%   Values for Bearing:
%
%     Single Number      == Angle [deg] to X-Axes
%     2-Element Vector   == [ dx  dy ],       where Angle = atan2( dy , dx )
%     [ 2 x 2 ] Matrice  == [ x0 y0 ; x1 y1 ] where Angle = atan2( y1-y0 , x1-x0 ) 
%
%   NaN or empty Option  == auto:
%
%    Angle = atan2( Y0(max(X0))-Y0(min(X0)) , max(X0)-min(X0) );  Orientation == 'X'
%    Angle = atan2( X0(max(Y0))-X0(min(Y0)) , max(Y0)-min(Y0) );  Orientation == 'Y'
%
%
%  Note: Angels with oddinary multiplys of 90deg (pi/2) for Mode == 'X' or
%                    even multiplys of 90deg for Mode == 'Y' are not allowed!
%
%
%  Multiple Outputs of Vector V and Angle:
%
%   [Y,X,Angle,Orientation] = COS_SURF( XY , X , ... );
%   [Y,X,Angle,Orientation] = COS_SURF( XY , N , ... );
%
%
% defaults:         Bearing = auto
%            CosineExponent = 1     ( if NaN, ZERO or empty Option )
%
%


Nin  = nargin;
Nout = nargout-1;

varargout = cell(1,Nout);

if Nin < 2
   error('Not enough InputArguments.')
end

%*******************************************************
% Check for Action
if ischar(xy)
   [y,out] = cosline(xy,x,varargin{:});
   n = min(Nout,size(out,2));
   varargout(1:n) = out(1:n); 
   return
end

%*******************************************************
% Get Option and Mode
 
opt  = [];
mode = '';

for ii = 1 : Nin-2
   v = varargin{ii};
   if ischar(v);
      mode = v;
   elseif isnumeric(v)
      opt = v;
   end
end

mode = mode( 1 : ~isempty(mode) );

mode = strcmp( lower(mode) , 'y' );

if Nout >= 3
   t = 'XY';
   varargout{3} = t(1+mode);
end

%*******************************************************
% Check XY-Line

s  = size(xy);
ok = ( ( prod(s) == s(1)*s(2) )  &  ( prod(s) >= 4 ) );
if ok
   ok = any( s == 2 );
   if ok  &  ( s(2) == 2 )  &  ~( s(1) == 2 )
      xy = permute(xy,[2 1]);
   end
end

if ~ok
   error('XY must define a Line, 2-Rows or 2-Columns with min. 2 Points.');
end

%---------------------------------------------------
% Flip with Mode == 1 (Y)

if mode
   xy = xy([2 1],:);
end

s = size(xy,2);

%*******************************************************
% Check X

if Nout > 0
   varargout{1} = x;
end

if isempty(x)
   y = zeros(size(x));
   return
end

si = size(x);

auto = ( prod(si) == 1 );

if auto
   n  = ceil(abs(x));
   si = [ 1  n ];
else
   if Nout >= 1
      varargout{1} = x;
   end
   x = x(:)';
end

%*******************************************************
% Check Exponent

if isempty(opt)

   e = 1;

else

   e = cat( 2 , abs(imag(opt(1))) , 1 );

   e = e( 1 + ( isnan(e(1)) | ( e(1) == 0 ) ) );

end

%*******************************************************
% Check Bearing

i1 = find( xy(1,:) == min(xy(1,:)) );
i2 = find( xy(1,:) == max(xy(1,:)) );


lim = cat( 2 , min(xy,[],2) , max(xy,[],2) );
 dl = diff(lim,1,2);

if isempty(opt)

   phi = atan2( dl(2) , dl(1) );

else
  
  phi = real(opt);
 
  isn = find(isnan(phi));

  sp = size(phi);

  if prod(sp) == 1
     phi = cat( 2 , phi , 0 );
     phi = phi( 1 + isnan(phi(1)) ) * pi/180;
  else
     if isequal(sp,[2 2])
        phi(isn) = lim(isn);
        phi = diff(phi,1,2);
        sp  = size(phi);
     elseif ( prod(sp) == 2 )
        phi(isn) = dl(isn);
     else
        error('Invalid Bearing.')
     end
     phi = atan2(phi(2),phi(1));
  end

end

if Nout >= 2
   varargout{2} = phi * 180/pi;
end

phi = (1-2*mode) * phi + pi/2 * mode;
 
%*******************************************************
% Check Bearing for 90deg-Parts

p = phi / ( pi/2 );

acc = 1e2*eps;

mp = mod(p,2);

if abs( mp -  round(mp) ) < acc

   mp = round(p);

   if mp == 1   
      t = 'XY';
      p = ( (1-mode) + [ 0  2 ] ) * 90;
      c = char(176);
      t = sprintf('%.0f%s or %.0f%s in Orientation %s',p(1),c,p(2),c,t(1+mode));
      error(['Singularity at Bearing ' t '.' ]);
   end

end

transform = ~( mp == 0 );

%*******************************************************
if transform

  %-------------------------------------
  % Basis CoordinateSystem

  ex0 = [ 1 ; 0 ];
  ey0 = [ 0 ; 1 ];

  %-------------------------------------
  % New CoordinateAxes in Basis System

  ex1 = [  cos(phi) ; sin(phi) ];
  ey1 = [ -sin(phi) ; cos(phi) ];

  %--------------------------------------------
  % TransformationMatrice

  E0 = [  ex0(:)  ey0(:)  ];
  E1 = [  ex1(:)  ey1(:)  ];

  T = E0  * inv( E1 );

  % Reverse: inv(T) = T .* [ 1  -1 ;  -1  1 ];

  %***************************************************
  if auto

     xy = T * xy;

     [h,s_i] = sort(xy(1,:),2);

     xy = xy(:,s_i);

     n = si(2);
     x = cat( 2 , xy(1,1)+(0:n-2)*(xy(1,s)-xy(1,1))/(n-1) , xy(1,s) );

  %***************************************************
  else

    %---------------------------------------------------
    % Sort new in Orientation

     [h,s_i] = sort( ( T(1,1)*xy(1,:) + T(1,2)*xy(2,:) ) , 2 );

     xy = xy(:,s_i);
 
    %---------------------------------------------------
    % Expand xy

     m  = tan(phi);

     x01 = [ min(x)  max(x) ];
     x01 = x01( [1 2] + [1 -1] * ( xy(1,1) > xy(1,s) ) );
  
     y01 = xy(2,[1 s]) + m * ( x01 - xy(1,[1 s]) ); 

     xy = cat( 2 , [ x01(1) ; y01(1) ] , xy , [ x01(2) ; y01(2) ] );

      s = s + 2;

     x0 = xy(1,:);

     if ~( all(diff(x0,1,2)<=0) | all(diff(x0,1,2)>=0) )
        t = 'XY';
        c = char(176);
        t = sprintf('Orientation %s at Angle %.0f%s',t(mode+1),phi*180/pi,c);
        error(['TargetPoints XY must be monotonic in ' t '.'  char(10)  ...
               'Use COS_LINE( XY , N , ... ) instead.' ])
     end

    %---------------------------------------------------
    % Transform

     xy = T * xy;

    %---------------------------------------------------
    % Check Expansion

     i0 = 1 + ( xy(1,2-1) >= xy(1,2) );
     i1 = s - ( xy(1,s-1) >= xy(1,s) );

     xy(2,[1 s]) = xy(2,[1 s]+[1 -1]);  % Y-Value at Border

     xy = xy( : , i0 : i1 );
     x0 = x0( i0 : i1 );

      s = i1 - i0 + 1;

    %---------------------------------------------------
    % Projection of X

      [x0,s_i] = sort(x0);

       ok = ones(1,s);
       ok( find( diff(x0,1,2) == 0 ) + 1 ) = 0;
       ok = find(ok);

      x = interp1(x0(ok),xy(1,s_i(ok)),x);

  end
  % auto

%*******************************************************
else

  [h,s_i] = sort(xy(1,:),2);

  xy = xy(:,s_i);

  if auto
     n = si(2);
     x = cat( 2 , xy(1,1)+(0:n-2)*(xy(1,s)-xy(1,1))/(n-1) , xy(1,s) );
     if Nout >= 1
        varargout{1} = x;
     end
   else
     x = min( max( x , xy(1,1) ) , xy(1,s) );
   end

end
%*******************************************************

%---------------------------------------------------
% Get Index x in xy(1,:)

  ok = ones(1,s);
  ok( find( diff(xy(1,:),1,2) == 0 ) + 1 ) = 0;
  ok = find(ok);

  if prod(size(ok)) == 1
     t = 'XY';
     c = char(176);
     t = sprintf('Orientation %s at Angle %.0f%s',t(mode+1),phi*180/pi,c);
     error(['Minimum 2 unique TargetPoints XY for ' t ' required.'])
  end

  ii = interp1( xy(1,ok) , ok , x );

  ii = floor(ii);

%--------------------------------------------------
% Define Cosine

  i0  = ( 1 : (s-1) );
  i1  = cat( 2 , i0 , i0(end) );
 
 % ZeroValue
   y =  ( xy(2,i0+1) + xy(2,i0) ) / 2;
   y = y(i1); 

 % Amplitude
   a = -( xy(2,i0+1) - xy(2,i0) ) / 2;  
   a = a(i1); 

 % Period
   l =  ( xy(1,i0+1) - xy(1,i0) ) * 2; 
   l = l(i1);
 
 % Phase
   p = xy(1,i1);

%--------------------------------------------------
% Combine Cosine with Zero and Amplitude

   c = cos( 2*pi * ( ( x - p(ii) ) ./ l(ii) ) );
  
   y = y(ii) + a(ii) .* sign(c).*(abs(c).^e);


%*******************************************************
if transform
   
  T = inv(T);    % T .* [ 1  -1 ; -1  1 ]

  if auto & ( Nout >= 1 )

     y = T * [ x ; y ];

     varargout{1} = y(1,:);
   
     y = y(2,:); 

  else

     y = x*T(2,1) + y*T(2,2);

  end

end

%*******************************************************
% Reshape back

if ~isequal( si , size(y) )
   y = reshape(y,si);
end
