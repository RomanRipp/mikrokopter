function [upp,axesize,x,y] = ppunit(axe);

% PPUNIT  Returns AxesUnitPerPixel and AxesPixelSize
%
%  [ UnitPerPixel, PixelSize, XC , YC ] = PPUNIT( AxesHandle )
%
%  UnitPerPixel = [ XUnitsPerPixel  YUnitsPerPixel  ZUnitsPerPixel ] ;
%  PixelSize    = [ PixelLeft  PixelBottom  PixelWidth  PixelHight ];
%
%  XC / YC      = PixelCoordinates for Corners of AxesBox,
%                  (relative to [ PixelLeft  PixelBottom ]):
%
%                 [ (X0,Y0,Z0) (X1,Y0,Z0) (X1,Y1,Z0) (X0,Y1,Z0) ...
%                   (X0,Y0,Z1) (X1,Y0,Z1) (X1,Y1,Z1) (X0,Y1,Z1)     ]
%
%           with: XLim = [ X0  X1 ]; YLim = [ Y0  Y1 ]; ZLim = [ Z0  Z1 ]
%                 XDir = 'normal';   YDir = 'normal';   ZDir = 'normal';
%                 Projection = 'orthographic'
%
% PPUNIT use the AxesProperties Position, DataAspectratio and View;
%         assuming a orthographic Projection of Axes.
%
%--------------------------------------------------------------------------
%
%  Units --> Points:   
% 
%    [points] = [unit] / UnitsPerPixel / ScreenPixelPerInch * PointsPerInch
%
%   PointsPerInch = 72
%
%--------------------------------------------------------------------------
% Example: 
%
% figpos = [200 100 400 300];  % Pixel: [ Left Bottom Width Hight ]
%
% figure('units','pixels','position',figpos,'nextplot','add')
% set(gca,'color','none','box','on','view',[20 30], ...
%         'xlim',[10 20],'ylim',[20 30],'zlim',[30 40],'nextplot','add')
%
% xlabel('X-Axis'),  zlabel('Z-Axis'),  zlabel('Z-Axis')
%
% [upp,axepos,xc,yc] = ppunit(gca);
%
% axepos = ( axepos - 1 ) ./ figpos([3 4 3 4]); % normalized
% xc = ( xc - 1 ) / figpos(3); xc = xc + axepos(1);
% yc = ( yc - 1 ) / figpos(4); yc = yc + axepos(2);
%
%% Dummy Axes for Surrounding Box and Corners
% 
% axe = axes( 'units' , 'normalized' , 'position' , [0 0 1 1] , ...
%             'xlim'  , [0 1]  , 'ylim'  , [0 1] , 'clim' , [1 4] , ...
%             'view'  , [0 90]  , 'color' , 'none' , 'nextplot' , 'add' );
%
%% Surrounding Box
%
% axebox = axepos + [ 0  0 axepos([1 2]) ];
%
% plot( axebox([1 3 3 1 1]) , axebox([2 2 4 4 2]) , 'r--' );
%
%% Lower and Upper Corners
%
% patch( 'xdata' , xc(1:4) , 'ydata' , yc(1:4) , 'marker' , 'v' );
% patch( 'xdata' , xc(5:8) , 'ydata' , yc(5:8) , 'marker' , '^' );
%
% set( findobj(axe,'type','patch') , 'cdata' , (1:4) , ...
%      'edgecolor' , 'interp' , 'facecolor' , 'none' , ...
%      'linestyle' , 'none' , 'linewidth' , 1 , 'markersize' , 10 , ...
%      'markeredgecolor' , 'k' , 'markerfacecolor' , 'flat' );
%    
%    

if nargin == 0
   fig = get(0,'currentfigure');
   axe = get(fig,'currentaxes');
end

if isempty(axe)

   upp     = zeros(0,2);
   axesize = zeros(0,4);

   x = [];
   y = [];

   return

end


ok = ( isnumeric(axe) & ( prod(size(axe)) == 1 ) );
if ~ok
    ok = ishandle(axe);
    if ok
       ok = strcmp( get(axe,'type') , 'axes' );
    end
end

if ~ok
   error('Input must be an AxesHandle.');
end


  
axeuni  = get(axe,'units');
          set(axe,'units','pixels')
axesize = get(axe,'position');
          set(axe,'units',axeuni);


%********************************************
if ~isequal( get(axe,'view') , [ 0  90 ] );

  [upp,axesize,x,y] = ppunit3(axe,axesize);

  return

end
%********************************************


dx = diff(get(axe,'xlim'));
dy = diff(get(axe,'ylim'));

mode   = get(axe,'dataaspectratiomode');

if strcmp(mode,'manual')

   aspect = get(axe,'dataaspectratio');

   %  w/dx*ax == h/dy*ay 
   %
   % w = h * dx/dy * ay/ax;
   % h = w * dy/dx * ax/ay; 
   %

   pos = zeros(2,2);

   pos(1,2) = axesize(4);
   pos(1,1) = axesize(4) * dx/dy * aspect(2)/aspect(1);
   pos(2,1) = axesize(3);
   pos(2,2) = axesize(3) * dy/dx * aspect(1)/aspect(2);

   ii = find( sum( ( pos <= axesize([1 1],[3 4]) ) , 2 ) == 2 );

   pos = pos(ii(1),:);
 
   axesize([1 2]) = axesize([1 2])+axesize([3 4])/2-pos/2;
   axesize([3 4]) = pos;

end

 upp = [ dx  dy  0 ] ./ axesize([3 4 4]) ; % UnitPerPixel

%-----------------------------------
% Index for Cube

  ix = [ 1 2 2 1 1 2 2 1 ];
  iy = [ 1 1 2 2 1 1 2 2 ];

  x = [ 0  axesize(3) ];
  y = [ 0  axesize(4) ];
  
  x = x(ix);
  y = y(iy);

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [upp,axesize,x,y] = ppunit3(axe,axesize);

% PPUNIT  Returns AxesUnitPerPixel and AxesPixelSize for 3D-Axes
%


  xl = get(axe,'xlim');
  yl = get(axe,'ylim');
  zl = get(axe,'zlim');

   v = get(axe,'view');

  dx = diff(xl);
  dy = diff(yl);
  dz = diff(zl);


%-----------------------------------
% Index for Cube

  ix = [ 1 2 2 1 1 2 2 1 ];
  iy = [ 1 1 2 2 1 1 2 2 ];
  iz = [ 1 1 1 1 2 2 2 2 ];


  % Index for Axes in ix, iy, iz
  %  i01(1,:) --> i01(2,:)
  %       x  y  z

  i01 = [ 1  1  1 
          2  4  5 ];

%-----------------------------------

mode = get(axe,'dataaspectratiomode');

if strcmp(mode,'manual')

 x = xl(ix);
 y = yl(iy);
 z = zl(iz);

 aspect = get(axe,'dataaspectratio');

else

 x = ix-1;  % 0 | 1
 y = iy-1;
 z = iz-1;

 aspect = [ 1  1  1 ];

end

[x,y] = xy_proj( x , y , z , v , aspect );

% Normalize
x = x - min(x);
y = y - min(y);

xr = max(x);  % Range
yr = max(y);  % Range


if strcmp(mode,'manual')

   %  w/dx*ax == h/dy*ay 
   %
   % w = h * dx/dy * ay/ax;
   % h = w * dy/dx * ax/ay; 
   %

   pos = zeros(2,2);

   pos(1,2) = axesize(4);
   pos(1,1) = axesize(4) * xr/yr;
   pos(2,1) = axesize(3);
   pos(2,2) = axesize(3) * yr/xr;

   ii = find( sum( ( pos <= axesize([1 1],[3 4]) ) , 2 ) == 2 );

   pos = pos(ii(1),:);

   axesize([1 2]) = axesize([1 2])+axesize([3 4])/2-pos/2;
   axesize([3 4]) = pos;

end

 % Cube in Pixel

 x = x/xr * axesize(3);
 y = y/yr * axesize(4);

 % AxesLength in Pixel

 d01 = sqrt( ( x(i01(1,:)) - x(i01(2,:)) ) .^2 + ...
             ( y(i01(1,:)) - y(i01(2,:)) ) .^2   );
       

 upp = ~( d01 < 1 ) .* [ dx  dy  dz ] ./ ( d01 + ( d01 < 1 ) );

 
%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [x,y,z] = xy_proj(x,y,z,v,aspect)

% function  [x,y,z] = xy_proj(X,Y,Z,View,DataAspectRatio)
%
% View = [ Azimut Elevation ]
%

si = num2cell(size(x));

az = v(1)*pi/180; 
el = v(2)*pi/180;

aspect = aspect(:)' / min(aspect(:));

x = x(:)';
y = y(:)';
z = z(:)';

T = [ cos(az)           sin(az)           0
     -sin(el)*sin(az)   sin(el)*cos(az)  cos(el)
      cos(el)*sin(az)  -cos(el)*cos(az)  sin(el) ];

T = T ./ aspect(ones(3,1),:);

xyz = T * [ x ; y ; z ];

x = reshape(xyz(1,:),si{:});
y = reshape(xyz(2,:),si{:});
z = reshape(xyz(3,:),si{:});
 
