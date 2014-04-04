function  [hh,dd] = stripes(h,varargin);

% STRIPES  draw stripes in patches
%
% H = STRIPES( PatchHandle , Angle , LineDistance ,  LineWidth , LineStyle )
%
% PatchHandle   Handle for PatchObject, one Polygon per each Column
%                of XData and YData
%
% Angle         Angle of Stripes to X-Axes
% LineDistance  LineDistance in Points
% LineWidth     LineWidth in Points
% LineStyle     LineStyle
%
% H             Handle of LineObject for Stripes
%                one Line per Polygon, NaN if the Patch
%                is not a closed Polygon or, however, 
%                no Stripe was drawn.
%
% use a negative Value of LineDistance for Units of X-Axis
%
%   ... , -LineDistance[X_Units] , ...
%
% use a NonZero imaginary part to LineDistance, to specify
%  the Number of StripeLines to draw:
%
%   ... ,  LineNumber + i , ...
%
%
% STRIPES( PatchHandle , Angle , LineDistance , Target , ... )
% 
% Creates the Stripes, that the Distances of the Lines to the
%   Target: [ X  Y ] is integer-multiple of LineDistance.
%
%---------------------------------------------------------------------
%
% XYL = STRIPES( XY , Angle , LineDistance , [Target] )
%
%  Returns the Koordinates of Stripes for a Polygon,
%   specified by XY ( 2 Rows or 2 Columns )
%
% XYL = [ 2 by N by 2 ];  N Stripes: [ P12  by  N  by  XY ]
%
%---------------------------------------------------------------------
%
% [ ... , Distance ] = STRIPES( ... ) returns the used Distance
%
%     between the Stripes (scaled if Distance in points),
%     that can be useful for using  LineNumber + i
%  
%---------------------------------------------------------------------
%
% example:
%
%  figure
%
%   h  = patch(rand(10,1),rand(10,1),'r','linewidth',2);
%   hh = stripes(h,30,10,2);
%
% % or, if   kochisle.m   is available
%
%   % Closed regular Triangle
%   X  = [ 0 ;    0.5    ; 1 ; 0 ];
%   Y  = [ 0 ; sqrt(3)/2 ; 0 ; 0 ];
%   
%   % Generator _/\_
%   XG = [ 0 ; 1/3 ;   1/2     ; 2/3 ; 1 ];
%   YG = [ 0 ;  0  ; sqrt(3)/6 ;  0  ; 0 ]; 
%  
%   % SnowFlakeCurve
%   [XF,YF] = kochisle(X,Y,XG,YG,3);
%
%   figure
%   axes('dataaspectratio',[1 1 1], ...
%        'xlim' , [ -0.5  1.5 ] , 'nextplot' , 'add' )
%
%   h  = patch(XF,YF,'r','linewidth',1);
%
%   h1 = stripes(h, 30,10,1);
%   h2 = stripes(h, 90,10,1);
%   h3 = stripes(h,150,10,1);
%
%%-----------------------------------------------------
%
%   xv = rand(6,1); yv = rand(6,1);
%   xy = stripes([xv yv],40,20+i);
%   ab = stripes([xv yv],00,0.1);
%
%   figure, hold on
%
%   patch(xv,yv,'w');
%   patch( 'xdata' , xy(:,:,1) , ...
%          'ydata' , xy(:,:,2) , ...
%          'facecolor' , 'none' , ...
%          'edgecolor' , 'r' );
%   patch( 'xdata' , ab(:,:,1) , ...
%          'ydata' , ab(:,:,2) , ...
%          'facecolor' , 'none' , ...
%          'edgecolor' , 'g' );
%  

Nin = nargin;

if Nin == 0
   error('Not enough InputArguments.');
end

[msg,ish,isn,x,y,p,ld,lw,ls,T,scl,axe] = checkin(h,varargin{:});

if ~isempty(msg)
    error(msg);
end

%********************************************

hh = [];
dd = [];

if ~isequal(size(x),size(y)) | isempty(x)
    return
end

%********************************************

M = size(x,2);
N = size(x,1);

if M <= 2
   return
end

if ish

   HoldState = get(axe,'nextplot');
               set(axe,'nextplot','add');

   edgcol = get(h,'edgecolor');
   hh     = NaN*ones(N,1);
   dd     = NaN*ones(N,1);

else

   hh = zeros(2,0,2);

   dd = [];

end

% Index to extract [Start;End] from 1. Line of Row

i01 = [ ( 1 : M ) ; [ ( 2 : M ) 1 ] ];
i01 = 1 + 2 * ( i01 - 1 );

% [ 1 3 5 ... 2*M-1 ; ...
%   3 5 7 ...   1        ]
 
for ii = 1 : N

   [xx,yy,dd(ii)] = get_stripe(x(ii,:),y(ii,:),scl,T,p,ld,isn,ish);

   if ~isempty(xx)

       if ish

           hh(ii) = patch( 'parent'    , axe    , ...
                           'xdata'     , xx     , ...
                           'ydata'     , yy     , ...
                           'linestyle' , ls     , ...
                           'linewidth' , lw     , ...
                           'marker'    , 'none' , ...
                           'edgecolor' , edgcol , ...
                           'facecolor' , 'none'        );

       else

           hh = cat(3,xx,yy);

       end

   end


end

if ish
   set(axe,'nextplot',HoldState);
end


%***************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,ish,isn,x,y,p,ld,lw,ls,T,scl,axe] = checkin(h,varargin);

% varargin = { ang ld [p] lw ls }

Nin = nargin - 1;

msg = cell(0,1);

ish = 0;
isn = 0;

x   = [];
y   = [];
p   = [];
ld  = [];
lw  = [];
ls  = '';

T   = zeros(2,2);

scl = ones(1,2);
axe = [];


%--------------------------------------------------------------------------------
% Handle | XY

if ~isnumeric(h)
    msg = cat(1,msg,{'First Input must be numeric.'});
else
    ish = ( prod(size(h)) == 1 );
    if ish
       ish = ishandle(h);
       if ish
          ish = strcmp(get(h,'type'),'patch');
       end
       if ~ish
           msg = cat(1,msg,{'Single first Input must be a Handle of an PatchObject.'});
       end
    elseif ~isempty(h)
       si = size(h);
       if ~( any( si(1:2) == 2 ) & ( si(1)*si(2) == prod(si) ) )
           msg = cat(1,msg,{'MatriceInput must have 2 Rows or 2 Columns for X and Y.'});
       end
    end
end

%--------------------------------------------------------------------------------
% Angle

if Nin < 1
   ang = [];
else
   ang = varargin{1};
end

if isempty(ang)
   ang = 0;
else
   ok = ( isnumeric(ang) & ( prod(size(ang)) == 1 ) );
   if ok
      ok = isfinite(ang);
   end
   if ~ok
       msg = cat(1,msg,{'Angle must be a single finite Numeric.'});
   end
end

%--------------------------------------------------------------------------------
% LineDistance

if Nin < 2
   ld = [];
else
   ld = varargin{2};
end

if isempty(ld)
   isn = 1;
   ld  = 9;
else
   ok = ( isnumeric(ld) & ( prod(size(ld)) == 1 ) );
   if ok
      isn = ~( imag(ld) == 0 );
       ld = real(ld);
       ok = ( ~( ld == 0 ) & isfinite(ld) );
   end
   if ~ok
       msg = cat(1,msg,{'LineDistance must be a single finite nonzero Numeric.'});
   elseif isn 
       ld = ceil(abs(ld));    % !!! Number of Stripes
   end
end

isd = ( ld < 0 );
ld  = abs(ld);

%--------------------------------------------------------------------------------
% Position

if Nin >= 3
   v = varargin{3};
   if isnumeric(v) & ( prod(size(v)) == 2 )
      p = v(:);
      varargin = varargin([ 1 2 (4:Nin) ]);
      Nin      = Nin - 1;
   end
end

%--------------------------------------------------------------------------------
if ish
%--------------------------------------------------------------------------------
 
   % LineWidth
   
   if Nin < 3
      lw = [];
   else
      lw = varargin{3};
   end

   if isempty(lw) 
      lw = get(h,'LineWidth');
   else
      ok = ( isnumeric(lw) & ( prod(size(lw)) == 1 ) );
      if ok
         ok = ( isfinite(lw) & ( lw > 0 ) );
      end
      if ~ok
          msg = cat(1,msg,{'LineWidth must be a single finite positive Numeric.'});
      end
   end

   % LineStyle

   if Nin < 4
      ls = '';
   else
      ls = varargin{4};
   end

   if isempty(ls)
      ls = get(h,'LineStyle');
   else
      ok = ( ischar(ls) & ( prod(size(ls)) == size(ls,2) ) );
      if ~ok
          msg = cat(1,msg,{'LineStyle must be a nonempty String.'});
      end
   end

%--------------------------------------------------------------------------------
end
%--------------------------------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n','Invalid Inputs.',msg{:});
else
    msg = '';
end

if ~isempty(msg) | isempty(h)
    return
end

%********************************************
% Get X and Y, One Patch per Row

if ish

   x = get(h,'xdata');
   y = get(h,'ydata');

   x = permute(x,[2 1]);
   y = permute(y,[2 1]);

   a = get(h,'parent');
   if xor( strcmp( get(a,'xdir') , 'reverse' ) , ...
           strcmp( get(a,'ydir') , 'reverse' )       )
      ang = -ang;
   end

else

   if ~( si(1) == 2 ) & ( si(2) == 2 )
      h = permute(h,[2 1]);
   end

   x = h(1,:);
   y = h(2,:);

end

%********************************************
% LineDirection

ang = ang - 180 * floor( ang / 180 );  % --> [   0 .. 180 )
ang = ang - 180 * ( ang > 90 );        % --> ( -90 ..  90 ]

ang = ang * pi/180;

%-------------------------------------
% Basis CoordinateSystem

ex0 = [ 1 ; 0 ];
ey0 = [ 0 ; 1 ];

%-------------------------------------
% New CoordinateAxes in Basis System

ex1 = [  cos(ang) ;  sin(ang) ];
ey1 = [ -sin(ang) ;  cos(ang) ];

%--------------------------------------------
% TransformationMatrice

E0 = [  ex0(:)  ey0(:) ];
E1 = [  ex1(:)  ey1(:) ];

T = E0  * inv( E1 );

%********************************************
% DataUnits in Points

if ~ish
    return
end

axe = get(h,'parent');

[upp,axesize] = ppunit(axe);         % Data per Pixel

if isd

   scl = [ 1   upp(2)/upp(1)  ];        % Use Distance true to X-Axis

else

   spp = get(0,'screenpixelsperinch');  % Pixel per Inch
   ppi = 72;                            % Point per Inch

   ppp = spp/ppi;                       % Pixel per Point

   scl = spp/ppi * upp;                 % Data per Point

end

%***************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [upp,axesize,figpos] = ppunit(axe);

% PPUNIT  Returns AxesUnitPerPixel  and AxesPixelSize
%
%  [ UnitPerPixel, AxePixelSize, FigurePixelSize ] = PPUNIT( AxesHandle )
%
%  UnitPerPixel = [ XUnitsPerPixel  YUnitsPerPixel ] ;
%  PixelSize    = [ PixelLeft  PixelBottom  PixelWidth  PixelHight ];
%
 

fig = get(axe,'parent');

 
fig_uni = get(fig,'units');
set(fig,'units','pixels')
figpos = get(fig,'position');
set(fig,'units',fig_uni);

axe_uni = get(axe,'units');
set(axe,'units','normalized')
axepos = get(axe,'position');
set(axe,'units',axe_uni);

dx = diff(get(axe,'xlim'));
dy = diff(get(axe,'ylim'));


axesize = axepos.*figpos([3 4 3 4]);

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

   ii = find( sum( ( pos <= [1;1]*axesize([3 4]) ) , 2 ) == 2 );

   pos = pos(ii(1),:);
 
   axesize([1 2]) = axesize([1 2])+axesize([3 4])/2-pos/2;
   axesize([3 4]) = pos;

end

 upp = [ dx  dy ] ./ axesize([3 4]) ; % UnitPerPixel


%***************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [x,y,ld] = get_stripe(x,y,scl,T,p,ld,isn,ish)


if any( isnan(x) | isnan(y) ) 
  x = [];
  y = [];
  return
end

n = size(x,2);

% Index to extract [Start;End] from 1. Line of Row
% [ 1 3 5 ... 2*M-1 ; ...
%   3 5 7 ...   1        ]
 
ii = [ ( 1 : n ) ; [ ( 2 : n ) 1 ] ];
ii = 1 + 2 * ( ii - 1 );

%------------------------------------------------
% Rotate, weighted with Scale before !!!
        
xy = T * [ x/scl(1) ; y/scl(2) ];

if ~isempty(p)
    p = T * (p./scl(:));
end

%------------------------------------------------
% Determine the Lines
%  [ Start  End ]

xx = xy(ii+0)';   % from 1. Row  
yy = xy(ii+1)';   % from 2. Row
 
jj = find( yy(:,1) == yy(:,2) );  % horizontal Lines

xx(jj,:) = [];
yy(jj,:) = [];

yl = [ min(xy(2,:))  max(xy(2,:)) ] + [ 1  -1 ] * (1e2*eps);

if isn
    nl = ld;
    ld = ( yl(2) - yl(1) ) / ( nl + 1 );
else
    if ish
       ld = min( ld , ( yl(2) - yl(1) )/2 );  % Minimum 1 Stripe
    end
    nl = ceil( ( yl(2) - yl(1) + 1e3*eps ) / ld ) - 1;  % Number of YLevels
end

if nl == 0
   x = [];
   y = [];
   return
end

% Y-Value for Stripes
if isempty(p)
   yl = yl(1) + ld * ( 1 : nl ); 
else
   y1 = yl(2);
   y0 = p(2) + ld * ceil( ( yl(1) - p(2) ) / ld );
   yl = y0 + ld * ( ( 0 : nl ) + ( abs(y0-yl(1)) < 1e2*eps ) );
   nl = min( nl , sum( yl < y1 ) );
   if nl == 0
      x = [];
      y = [];
      return
   end
   yl = yl(1:nl);
end

s1 = size(xx,1);    % New Number of Lines <= M
i1 = ones(1,s1);
i2 = ones(1,nl);

% Calculate X-Value of each Line at yl


yl = yl(i1,:);

m = (xx(:,2)-xx(:,1)) ./ (yy(:,2)-yy(:,1));
n =  xx(:,1) - m .* yy(:,1);

xl = m(:,i2) .* yl + n(:,i2);   
 
% Sort Y01-Values of PatchLines

ii = ( 1 : s1 )';
ii = ii(:,[1 1]);

jj      = ( yy(:,1) < yy(:,2) );
ii(:,1) = ii(:,1) + ~jj * s1;
ii(:,2) = ii(:,2) +  jj * s1 ;

xx = xx(ii);
yy = yy(ii);


% Find correct Lines, matching yl

ii = ( ( yy(:,1*i2) <= yl )  &  ( yl < yy(:,2*i2) ) );


% Sort X-Values

[ xl , jj ] = sort(xl,1);

jj = jj + ones(s1,1) * ( 0 : nl-1 ) * s1;

yl = yl(jj);
ii = ii(jj);

kk = cumsum( sum(ii(:,1:(end-1)),1) , 2 );  % see below

ii = find(ii(:));
nn = size(ii,1);

% Build the PointPairs

jj = ( 1 : (nn-1) );
ii = [ ii(jj)  ii(jj+1) ];

xl = xl(ii);
yl = yl(ii);

if nn == 2  % Just 1 Pair
   xl = xl(:)';
   yl = yl(:)';
end

% After each nl-Pairs the yl-Level switch ==> NaN-Line

kk( find( kk > size(xl,1) ) ) = [];

xl(kk,:) = [];
yl(kk,:) = [];

if isempty(xl)
   x = []; 
   y = [];
   return
end

% Look, if MeanPoint of Lines inside Polygon

ii = inside( (xl(:,1)+xl(:,2))/2 , (yl(:,1)+yl(:,2))/2 , xy(1,:) , xy(2,:) );

ii = find(ii);

xl = xl(ii,:)';
yl = yl(ii,:)';

xy1 = T' * [ xl(1,:) ; yl(1,:) ];
xy2 = T' * [ xl(2,:) ; yl(2,:) ];
    
x = scl(1) * [ xy1(1,:) ; xy2(1,:) ];
y = scl(2) * [ xy1(2,:) ; xy2(2,:) ];


%***************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function in = inside(x,y,xv,yv)

% INSIDE True for points inside polygonal region.
%
% IN = INSIDE(X, Y, XV, YV) returns a matrix IN the size of
%   X and Y.  IN(p,q) = 1 if the point (X(p,q), Y(p,q)) is
%   strictly inside the polygonal region whose vertices are
%   specified by the vectors XV and YV; 
%   %%  IN(p,q) is 0.5 if the point is on the polygon; 
%   otherwise IN(p,q) = 0.
%
%   Example:
%
% xv = rand(6,1); yv = rand(6,1);
% xv = [xv ; xv(1)]; yv = [yv ; yv(1)];
% x = rand(1000,1); y = rand(1000,1);
% in = inside(x,y,xv,yv);
% plot(xv,yv,x(in),y(in),'.r',x(~in),y(~in),'.b')

%   Steven L. Eddins, July 1994
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 1.5 $  $Date: 1997/11/21 23:40:39 $
%
%   $Revision: 1.6 $  $Date: 2004/03/21 22:15:00 $  Ch.Begler
%     110% Time   in compare with INPOLYGON
%      50% Memory in compare with INPOLYGON


% PolygonVertices: RowVectors

n  = size(xv,2);

% Close the Polygon

if ~( ( xv(1) == xv(n) ) & ( yv(1) == yv(n) ) )
  xv = [ xv , xv(1) ];
  yv = [ yv , yv(1) ];
  n  = n + 1;
end

% Matrice: ColumnVectors

m = size(x,1);

xv = xv(ones(1,m),:) - x(:,ones(n,1));
yv = yv(ones(1,m),:) - y(:,ones(n,1));

%--------------------------------------------------------------
% Compute the Signum of the cross product and the dot product
% of adjacent vertices.

ii = ( 1 : n-1 );

% SignCrossProduct
scp = sign(  xv(:,ii  ) .* yv(:,ii+1) ...
           - xv(:,ii+1) .* yv(:,ii  )     );

% DotProduct
% dp = xv(:,ii) .* xv(:,ii+1) + yv(:,ii) .* yv(:,ii+1);

%--------------------------------------------------------------
% Compute the quadrant number for the vertices relative
% to the test points.

xv = ( xv > 0 );
yv = ( yv > 0 );

% Compute the vertex quadrant changes for each test point.

dq = diff( ( ~xv & yv ) + 2 * ( ~xv & ~yv ) + 3 * ( xv & ~yv ) , 1 , 2 );

%--------------------------------------------------------------
% Fix up the quadrant differences.  Replace 3 by -1 and -3 by 1.
% Any quadrant difference with an absolute value of 2 should have
% the same sign as the cross product.

% dq = dq - 4 * sign(dq) .* ( abs(dq) == 3 );
% dq = dq - 4 * ( dq == 3 ) + 4 * ( dq == -3 );

dq(find(dq== 3)) = -1;
dq(find(dq==-3)) =  1;

dq = dq + (  2 * scp - dq ) .* ( abs(dq) == 2 );

dq(find(isnan(dq))) = 0;

%--------------------------------------------------------------
% Find the points on the polygon.  If the cross product is 0 and
% the dot product is nonpositive anywhere, then the corresponding
% point must be on the contour.

in = ( sum(dq,2) ~= 0 );

% Indicate points on the polygon with a value of 0.5.
% in = in  - 0.5 * any( ( scp == 0 ) & ( dp <= 0 ) , 2 );

