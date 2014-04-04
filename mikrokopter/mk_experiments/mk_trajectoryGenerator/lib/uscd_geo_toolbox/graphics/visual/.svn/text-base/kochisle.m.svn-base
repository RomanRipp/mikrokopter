 function [x,y] = kochisle(xi,yi,xg,yg,N)

% KOCHISLE  Generates a KochIslandFractal (SnowFlakeCurve)
%
% [XF,YF] = KOCHISLE(X,Y,XG,YG,N)
%
% The Elements of the Curve [ X , Y ]
% will be N times replaced with the Elements
% of the Generator [ XG , YG ].
%
% All Vectors are handled as ColumnVectors.
% Note, that the Sense of the Curves
%  is important.
%
% [XF,YF] = KOCHISLE(X,Y,XG,YG)
%   use the default: N = 6
%
% [XF,YF] = KOCHISLE(N)
%   use the defaults from the following Example
%
%------------------------------------------------
%
% Example:
%
%  % Closed regular Triangle
%  X  = [ 0 ;    0.5    ; 1 ; 0 ];
%  Y  = [ 0 ; sqrt(3)/2 ; 0 ; 0 ];
%  
%  % Generator _/\_
%  XG = [ 0 ; 1/3 ;   1/2     ; 2/3 ; 1 ];
%  YG = [ 0 ;  0  ; sqrt(3)/6 ;  0  ; 0 ]; 
%
%  % SnowFlakeCurve
%  [XF,YF] = kochisle(X,Y,XG,YG,4);
%  plot(XF,YF)
%
%  % try also 
%  [XF,YF] = kochisle(flipud(X),flipud(Y),XG,YG,4);
%
%

Nin = nargin;

if ~any( Nin == [ 0  1  4  5 ] )
    error('Invalid Number of InputArguments.');
end

%************************************************

if Nin <= 1

   if Nin == 0
      N = [];
   else
      N  = xi;
   end

   % Closed regular Triangle
   xi = [ 0 ;    0.5    ; 1 ; 0 ];
   yi = [ 0 ; sqrt(3)/2 ; 0 ; 0 ];
   
   % Generator _/\_
   xg = [ 0 ; 1/3 ;   1/2     ; 2/3 ; 1 ];
   yg = [ 0 ;  0  ; sqrt(3)/6 ;  0  ; 0 ]; 

elseif Nin == 4

   N  = [];

end
   
if isempty(N)
   N = 6;
end

%************************************************

msg = cell(0,1);

%************************************************
% make ColumnVector of Initiator
%------------------------------------------------
if size(xi,1) == 1
 xi = xi';
else
 xi = xi(:,1);
end

if size(yi,1) == 1
 yi = yi';
else
 yi = yi(:,1);
end

if length(xi) ~= length(yi)
   msg = cat(1,msg,{'X and Y must have the same Length.'})
end

%************************************************
% Make ColumnVector of Generator
%------------------------------------------------
if size(xg,1) == 1
 xg = xg';
else
 xg = xg(:,1);
end

if size(yg,1) == 1
 yg = yg';
else
 yg = yg(:,1);
end

if length(xg) ~= length(yg)
   msg = cat(1,msg,{'XG and YG must have the same Length.'})
end


%************************************************

if ~( ( mod(N,1) == 0 ) & ( N >= 0 ) )
   msg = cat(1,msg,{'N must be a positive Integer.'})
end

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Invalid Inputs.\n%s',msg));
end

%************************************************

if N == 0
   x = xi;
   y = yi;
   return
end

%************************************************
% Extract Data from Generator
%------------------------------------------------

ng = length(xg);

% Generator symetrical to [ 0  0 ]
xg = xg - (xg(1)+xg(ng))/2;
yg = yg - (yg(1)+yg(ng))/2;

% NormalVector of GeneratorElements
pg = atan2(diff(xg),-diff(yg));

% LengthScale of Generator
dg0 = sqrt( (xg(ng)-xg(1)).^2 + (yg(ng)-yg(1)).^2 ); 

% NormalVector of Generator
pg0 = atan2( (xg(ng)-xg(1)) , -(yg(ng)-yg(1)) );


oo = ones(1,ng);


% Generation at first only for one BasisVector
% that makes the Code longer, but it saves Memory
% in the order of length(xi)-1

% Symetrical Basis
xb=[-0.5 0.5];
yb=[ 0   0 ];

x=xb;
y=yb;


%************************************************
% Main Loop
%------------------------------------------------

for ii = 1 : N

   n = length(x);

   % NormalVector of Initiators
   p = atan2( diff(x) , -diff(y) );

   % Length of Initiators
   d = sqrt( diff(x).^2 + diff(y).^2 );

   % RotationMatrix
   M = reshape( [ cos(p-pg0) -sin(p-pg0) sin(p-pg0) cos(p-pg0) ]' , ...
                 2 , 2*(n-1) )';
   clear p

   xy = M * [ xg  yg ]';

   clear M 

   % Scale and Move
   x = xy(1:2:2*(n-1),:).*[d/dg0*oo] + [x(1:n-1)+diff(x)/2] * oo;
   y = xy(2:2:2*(n-1),:).*[d/dg0*oo] + [y(1:n-1)+diff(y)/2] * oo;

   clear xy d 

    xe = x(n-1,ng);
    ye = y(n-1,ng);
    x(:,ng) = [];
    y(:,ng) = [];

   % Make ColumnVector
   x = reshape( x' , (ng-1)*(n-1) , 1 );
   y = reshape( y' , (ng-1)*(n-1) , 1 );

    x( (ng-1)*(n-1)+1 ) = xe;
    y( (ng-1)*(n-1)+1 ) = ye;

end
% ii = 1 : N


%****************************************************
% Now bring the Vector [x y] "over" the Elements of [xi,yi]
%----------------------------------------------------

% LengthScale of Basis
db0 = sqrt( diff(xb).^2 + diff(yb).^2 ); 

% NormalVector of Basis
pb0 = atan2( diff(xb) , -diff(yb) );

ni = length(xi);
nb = length(x);

oo = ones(1,nb);

   % NormalVector of Initiators
   pi = atan2( diff(xi) , -diff(yi) );

   % Length of Initiators
   di = sqrt( diff(xi).^2 + diff(yi).^2 );

   % RotationMatrix
   M = reshape( [ cos(pi-pb0) -sin(pi-pb0) sin(pi-pb0) cos(pi-pb0) ]' , ...
                 2 , 2*(ni-1) )';
   clear pi

   xy = M * [ x  y ]';

   clear M 

   % Scale and Move
   x = xy(1:2:2*(ni-1),:).*[di/db0*oo] + [xi(1:ni-1)+diff(xi)/2] * oo;
   y = xy(2:2:2*(ni-1),:).*[di/db0*oo] + [yi(1:ni-1)+diff(yi)/2] * oo;

   clear xy di 

    xe = x(ni-1,nb);
    ye = y(ni-1,nb);
    x(:,nb) = [];
    y(:,nb) = [];

   % Make ColumnVector
   x = reshape( x' , (nb-1)*(ni-1) , 1 );
   y = reshape( y' , (nb-1)*(ni-1) , 1 );

    x( (nb-1)*(ni-1)+1 ) = xe;
    y( (nb-1)*(ni-1)+1 ) = ye;

if nargout == 0
   figure, axis equal, hold on, box on, plot(x,y)
   clear x y
end