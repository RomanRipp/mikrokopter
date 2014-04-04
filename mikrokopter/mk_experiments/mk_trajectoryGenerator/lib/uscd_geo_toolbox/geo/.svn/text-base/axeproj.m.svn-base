function [s,hb,hg] = axeproj(axe,fcn,c,p,n);

% AXEPROJ   Projection of Axes using Mapping Toolbox
%
% [MapStruct,HB,HG] = AXEPROJ( AXE , FCN , Origin , [Parallels] )
%
% FCN    = Function to call of MappingToolbox/MAPPROJ
%          specials: 'azma'  --> 'eqaazim'
%                    'azmd'  --> 'eqdazim'
%                    'cona'  --> 'eqaconic'
%                    'cond'  --> 'eqdconic'
%                    'cyla'  --> 'eqacylin'
%                    'cyld'  --> 'eqdcylin'
%                    'lamb'  --> 'lambert'
%                    'merc'  --> 'mercator'
%
% Origin = [ OriginLon OriginLat [Rotation] ]
%
% MapStruct = 
%    geoid         : [ 6378137 0.08181919131097 ]  % WGS84
%    aspect        : 'normal'
%    origin        : Origin([2 1 [3]])
%    mapparallels  : Parallels
%    flatlimit     : YLim
%    flonlimit     : XLim
%    angleunits    : 'degrees'
%    scalefactor   : 1
%    falseeasting  : 0
%    falsenorthing : 0
%
% Note:  AXEPROJ requires the NON-MappingToolbox-Function MERCATOR.
%        To avoid confusions with the Function MERCATOR of the
%         MappingToolbox, this Function is added as Subfunction MERCAT
%         at the end of this File.
%
% see also: MPROJ, MERCATOR (required, NON-Mappingtoolbox)
%


Nin = nargin;

app = mfilename;
fsp = ( app == filesep );
if any(fsp)
    app = app(max(find(fsp))+1:end);
end
app = upper(app);

dp = 5;  % PixelResolution for Grid and BorderPatches

spt = struct( 'trimmed' , {[]} , ...
              'clipped' , {[]}        );

s = struct( ...
'projection'    , { '' }        , ...
'geoid'         , { [ 6378137 0.08181919131097 ] } , ...
'aspect'        , { 'normal' }   , ...
'origin'        , { zeros(1,2) } , ...
'mapparallels'  , { zeros(1,1) } , ...
'flatlimit'     , { zeros(0,2) } , ...
'flonlimit'     , { zeros(0,2) } , ...
'angleunits'    , { 'degrees'  } , ...
'scalefactor '  , { 1 } , ...
'falseeasting'  , { 0 } , ...
'falsenorthing' , { 0 } , ...
'savedpoints'   , { spt }      );

%*****************************************************
% Check Inputs

if Nin < 2
   error('Not enough InputArguments.');
end

if Nin < 3
   n = [];
end

if Nin < 4
   c = [];
end

if Nin < 5
   p = [];
end

%-----------------------------------------------------

[msg,s,n] = checkin(s,axe,fcn,c,p,n,Nin);

if ~isempty(msg)
    error(msg);
end

fcn = s.projection;

%*****************************************************
% MapBorder, Limits, Ticks

if isempty(n)
   n = {};
else
   n = {n};
end

hb = border(axe,n{:});

xl = get(axe,'xlim');
yl = get(axe,'ylim');

xt = get(axe,'xtick');
yt = get(axe,'ytick');

xtl = get(axe,'xticklabel');
ytl = get(axe,'yticklabel');

[upp,si] = ppunit(axe);

upp = upp*dp;

si  = ceil(si([3 4])/dp);   % [ Width Height ] Pixel

%*********************************************************
% Set Missing ProjectionProperties

if isempty(s.flonlimit)
   s.flonlimit = xl;
end

if isempty(s.flatlimit)
   s.flatlimit = yl;
end

if isempty(s.origin)
   s.origin = [ mean(yl) mean(xl) 0 ];
end

if isempty(s.mapparallels)
   s.mapparallels = s.origin(1);
end

%*********************************************************
% Boundary

xb = cat(2,xl,mean(xl),s.origin(2));
yb = cat(2,yl,mean(yl),s.origin(1));

[xb,yb] = meshgrid(xb,yb);

[xb,yb] = feval(fcn,s,yb,xb,'none','forward',spt);

xb = [ min(xb(:)) max(xb(:)) ];
yb = [ min(yb(:)) max(yb(:)) ];


%*********************************************************
% Transform XY-Data of Children of Axes

%---------------------------------------------
% Set LabelUnits to normalized !!!

lb = { 'xlabel' 'ylabel' 'zlabel' 'title' };

nl = size(lb,2);
hl = zeros(1,nl);
ul =  cell(1,nl);

for ii = 1 : nl
    hl(ii) = get(axe,lb{ii}); 
    ul{ii} = get(hl(ii),'units');
end

set(hl,'units','normalized');

%---------------------------------------------
% Project AxesChildren

set_data(axe,s);

%---------------------------------------------
% Set AxesLimits, ApplicationData

set( axe , 'xlim' , xb , ...
           'ylim' , yb , ...
        'visible' , 'off'   );

setappdata(axe,app,s);

%---------------------------------------------
% Set LabelUnits back

for ii = 1 : nl
    set(hl(ii),'units',ul{ii});
end

%---------------------------------------------
% Grid

hg = set_grid(axe,s,si,xl,yl,xt,yt);

%*********************************************************
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  hb = [];

  return

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%*********************************************************
% Nice Border

nx = size(xt,2);
ny = size(yt,2);

[xt,yp] = feval(fcn,s,yb(:)*ones(1,nx),xt([1 1],:),'none','forward',spt);
[xp,yt] = feval(fcn,s,yt([1 1],:),xb(:)*ones(1,ny),'none','forward',spt);

bw = get_int(axe,xb,yb);

keyboard

figure,axis([-20 60 30 70]),grid on, axeproj(gca,'cona')

for h = hb(:)'

    x = get(h,'xdata');
    y = get(h,'ydata');

    dx = diff(x,1,1);
    dy = diff(y,1,1);

    n = ceil( max( abs(dx)/upp(1) , abs(dy)/upp(2) ) );

    [x,y] = feval(fcn,s,y,x,'none','forward',spt);

    set( h , 'xdata' , x , ...
             'ydata' , y , ...
             'tag'   , 'BORDER'  );

end


%******************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,s,n] = checkin(s,axe,fcn,c,p,n,Nin);


msg = cell(0,1);

%-----------------------------------------------------
% Axes

ok = ( isnumeric(axe) & ( prod(size(axe)) == 1 ) );
if ok
   ok = ishandle(axe);
   if ok
      ok = strcmp(get(axe,'type'),'axes');
   end
end

if ~ok
    msg = cat(1,msg,{'First Input must be a single AxesHandle.'});
end

%-----------------------------------------------------
% FCN | MapStruct

si = size(fcn);

is_fcn = (   ischar(fcn) & ~isempty(fcn) & ( prod(si) == si(2) ) );
is_str = ( isstruct(fcn) & ( prod(si) == 1 ) );

if is_fcn

   fcn = lower(fcn);

   if si(2) == 4 
      switch fcn
       case 'lamb'
         fcn = 'lambert';
       case 'merc'
         fcn = 'mercator'
       otherwise
         suf = '';
         switch fcn(1:3)
          case 'azm'
                suf = 'azim';       
          case 'con'
                suf = 'conic';
          case 'cyl'
                suf = 'cylin';
        end
        if ~isempty(suf)
            fcn = sprintf('eq%s%s',fcn(4),suf);
        end
      end
   end

   if exist(fcn,'file') == 2
      s.projection = fcn;
   else
      msg = cat(1,msg,{sprintf('Function "%s" not found.',fcn)});
   end

elseif is_str

    [m,s] = structcmp( s , fcn , 1 )

    if ~isempty(m)
        msg = cat(1,msg,{sprintf('Invalid MapStructure.\n%s',m)});
    end

else

    msg = cat(1,msg,{'First Input must be a String or a MapStructure.'});

end

%-----------------------------------------------------
% Origin

isn = ( ( Nin == 3 ) & ( prod(size(c)) == 1 ) );
if isn
   isn = ( ( real(c) == 0 ) & ~( imag(c) == 0 ) );
end
if isn
   n = c;
   c = [];
end

if ~isempty(c)
    sc = prod(size(c));
    if ~( isnumeric(c) & ( sc >= 2 ) )
       msg = cat(1,msg,{'Origin be numeric with min. 2 Elements: [ Lon Lat ].'});
    else
       c = c( 1 : min(3,sc) );
       c = c(:)';
       if sc == 2
          c = [ c  0 ];
       end
       c([1 2]) = c([2 1]);     % [ Lat  Lon ]
       s.origin = c;
    end
elseif is_str
    c  = s.origin;
    sc = prod(size(c));
    c  = c( 1 : min(3,sc) );
    c  = c(:)';
    if sc == 2
       c = [ c  0 ];
    end
    s.origin = c;
elseif is_fcn
    s.origin = [];
end

%-----------------------------------------------------
% Parallel

isn = ( ( Nin == 4 ) & ( prod(size(p)) == 1 ) );
if isn
   isn = ( ( real(p) == 0 ) & ~( imag(p) == 0 ) );
end
if isn
   n = p;
   p = [];
end

if ~isempty(p)
    sp = prod(size(p));
    if ~( isnumeric(p) & ( sp <= 2 ) )
       msg = cat(1,msg,{'Parallels be numeric with max. 2 Elements.'});
    else
       s.mapparallels = p;
    end
elseif ~is_str & ~isempty(s.origin)
    s.mapparallels = s.origin(1);
else
    s.mapparallels = p;
end

%-----------------------------------------------------
% Intervall

if isn
   n = imag(n);
elseif ~isempty(n)
   sn = prod(size(n));
   if ~( isnumeric(p) & ( sn == 1 ) )
       msg = cat(1,msg,{'Intervall be a single numeric.'});
   elseif ( n == 0 )
       msg = cat(1,msg,{'Intervall must be NOT ZERO.'});
   elseif ( real(n) == 0 )
       n = imag(n);
   else
       n = real(n);
   end
end   

%-----------------------------------------------------

if isempty(msg)
   msg = '';
else
   msg = sprintf('%s\n',msg{:});
end

%******************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function set_data(axe,s)

% SET_DATA  Project XY-Data of AxesChildren


fcn = s.projection;
spt = s.savedpoints;

%---------------------------------------------
% ObjectInitialisation

ini = { 'image'       'surface'
        'surface'     'surface'
        'light'       'light'
        'text'        'text'
        'line'        'line'
        'patch'       'patch'
        'rectangle'   'patch'   };

%---------------------------------------------
% Get Children, transform positions

shh = get(0,'showhiddenhandles');
      set(0,'showhiddenhandles','on');
chl = get(axe,'children');
      set(0,'showhiddenhandles',shh);

if any(strcmp(get(chl,'type'),'image'));
   warning('Projection may not be correct for Data of Type IMAGE.');
end

for ch = chl(:)'

    typ = get(ch,'type');

    obj = strcmp(ini(:,1),typ);
    if any(obj)
       obj = find(obj);
       obj = ini{obj,2};
    else
       obj = 'none';
    end

    if strcmp(typ,'text')

       if strcmp( get(ch,'units') , 'data' )

          pos = get(ch,'position');
 
          [pos(1),pos(2)] = feval(fcn,s,pos(2),pos(1),obj,'forward',spt);

          set(ch,'position',pos);

       end

    elseif 1 %%% ~any( ch == hb )

         x = get(ch,'xdata');
         y = get(ch,'ydata');

         if ( prod(size(x)) == 1 ) & ( prod(size(y)) == 1 )
            obj = 'point';
         end

         [x,y] = feval(fcn,s,y,x,obj,'forward',spt);

         set(ch,'xdata',x,'ydata',y);

    end

end  

%******************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function hg = set_grid(axe,s,si,xl,yl,xt,yt)

% SET_GRID  Draw a Grid

gg = { get(axe,'xgrid')  get(axe,'ygrid') };
gg = strcmp(gg,'on');

if ~any(gg)
    hg = [];
    return
end


   set( axe , 'xgrid' , 'off' , ...
              'ygrid' , 'off'        );

   xg = xt( 1 : (end*gg(1)) );
   yg = yt( 1 : (end*gg(2)) );

   xn = gg(2) * si(1);
   yn = gg(1) * si(2);
 
   xx = zeros(xn,1);
   yy = zeros(yn,1);
   if gg(2)
      xx = linspace(xl(1),xl(2),xn)';
   end
   if gg(1)
      yy = linspace(yl(1),yl(2),yn)';
   end

   nx = size(xg,2);
   ny = size(yg,2);

   xg = cat(1,xg(ones(1,yn),:),NaN*ones(1,nx));
   yx = cat(1,yy(:,ones(1,nx)),NaN*ones(1,nx));

   yg = cat(1,yg(ones(1,xn),:),NaN*ones(1,ny));
   xy = cat(1,xx(:,ones(1,ny)),NaN*ones(1,ny));

   xg = cat(1,xg(:),xy(:));
   yg = cat(1,yx(:),yg(:));

   [xg,yg] = feval(s.projection,s,yg,xg,'none','forward',s.savedpoints);

   zl = get(axe,'zlim');
   zl = zl( 1 + strcmp(get(axe,'zdir'),'normal') );

   hg = line( 'parent'    , axe , ...
              'xdata'     , xg  , ...
              'ydata'     , yg  , ... 
              'zdata'     , zl + 0*xg , ...
              'linestyle' , get(axe,'gridlinestyle') , ...
              'linewidth' , get(axe,'linewidth')     , ...
              'color'     , 'k' , ...
              'tag'       , 'GRID'          );

%******************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function hb = border(axe,int)

% BORDER  Set Geographic AxesBorder
%
% HH = BORDER( Axe , MaxInt )
%
%  MaxInt < 0  ==> round Limits to even parts of Intervall
%
% HH = [ Corner WhiteX WhiteY BlackX BlackY ]  PatchHandles
%

if nargin < 2
   int = 5;
end

is_fix = ( int < 0 );

int = abs(int);

%*****************************************************************

xl = get(axe,'xlim');
yl = get(axe,'ylim');

%*****************************************************************
% Get Ticks

[xt,yt,tick_int,mode] = get_tick(axe,xl,yl,int);

%*****************************************************************
% Get Border

[bw,intv] = get_int(axe,xl,yl,xt,yt);

%*****************************************************************
% Check for AutoLimit

if is_fix & ~isempty(intv) & ~isempty(tick_int)

    acc = 1e-10;

    op = [ -1  1 ];

    xl = intv * op .* ceil( op.*(xl-2*op*acc) / intv );
    yl = intv * op .* ceil( op.*(yl-2*op*acc) / intv );

    if isequal(mode{1},'auto')
       if isempty(xt)
          xt = xl;
       else
          xc = xt([1 end]) + op*tick_int;
          ok = ( abs( xc - xl ) > acc );
          xt = cat( 2 , xc(1) , xt , xc(2) );
          xt = xt( 1+ok(1) : end-ok(2) );
       end
    end

    if isequal(mode{2},'auto')
       if isempty(yt)
          ytic = yl;
       else
          yc = yt([1 end]) + op*tick_int;
          ok = ( abs( yc - ylim ) > acc );
          yt = cat( 2 , yc(1) , yt , yc(2) );
          yt = yt( 1+ok(1) : end-ok(2) );
       end
    end

    bw = get_int(axe,xl,yl,xt,yt);

end

%*****************************************************************
% TickLabels

[xtl,ytl] = get_label(xt,yt);

%********************************************************
% Set Axes

 
set(axe, 'dataaspectratio' , [ 1  1  1 ] , ...
          'tickdir'     , 'out'     , ...
          'box'         , 'on'      , ...
          'xlim'        , xl        , ...
          'xtick'       , xt        , ...
          'xticklabel'  , xtl       , ...
          'ytick'       , yt        , ...
          'ylim'        , yl        , ...
          'yticklabel'  , ytl       , ...
          'nextplot'    , 'add'                   )  


%********************************************************
% BorderPatches

hb    = zeros(1,5); % Corner WhiteX WhiteY BlackX BlackY

lb = { 'C' 'XW' 'YW' 'XK' 'YK' };

hb(1) = patch( 'parent' , axe      , ...
               'tag'    , [ 'BORDER_' lb{1} ]  , ...
            'facecolor' , [ 1  1  1 ]          , ...
            'edgecolor' , [ 0  0  0 ]          , ...
            'marker'    , 'none'               , ...
            'linestyle' , get(axe,'linestyle') , ...
            'linewidth' , get(axe,'linewidth') , ...
             'clipping' , 'off'                , ...
            'erasemode' , 'background'               );

for ii = 2 : 5
    hb(ii) = copyobj(hb(1),axe);
    set( hb(ii) , 'tag' , [ 'Border_' lb{ii} ] , ...
            'facecolor' , lower(lb{ii}(2))           );
end

set_border(axe,hb,xl,yl,xt,yt,bw,intv);


%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [bw,intv] = get_int(axe,xl,yl,xt,yt);

bw   = [];
intv = [];

%********************************************************
% Get BorderWidth

  [mdiff,kk]= max([ diff(xl) diff(yl) ]);

   bw = 0.010*mdiff; % BorderWidth
 
   fig = get(axe,'parent');

   [upp,axepos,figpos] = ppunit(axe);

   papuni = get(fig,'paperunits');
   set(fig,'paperunits','inches');
   pappos = get(fig,'paperposition');
   set(fig,'paperunits',papuni);

   % Minimum 1.8 mm  == 10 Points (150 dpi)
   bmin1 = 1.8/25.4 / (pappos(4)*axepos(4)/figpos(4)) * diff(yl); 
   
   % Minimum 5 Pixels
   bmin2 = 5 / axepos(4) * diff(yl);        % Minimum 5 Pixels
   
   bmin = max( bmin1 , bmin2 );

   bw = bw + (bmin-bw)*( bw < bmin );

if nargout < 2
   return
end

%********************************************************
% Get Intervall

   jj = find( ( xt < xl(1) ) | ( xt > xl(2) ) );
   xt(jj) = [];

   jj = find( ( yt < yl(1) ) | ( yt > yl(2) ) );
   yt(jj) = [];

   dxt = diff(xt);
   dyt = diff(yt);

   dok_x = isempty(dxt);
   if ~dok_x
     dok_x = ( ( max(dxt) - min(dxt) ) < 1e-3*mean(dxt) );
   end

   dok_y = isempty(dyt);
   if ~dok_y
     dok_y = ( ( max(dyt) - min(dyt) ) < 1e-3*mean(dyt) );
   end


if  ~( dok_x  &  dok_y )
     return
end


   mdxt = [];
   if ~isempty(dxt)
      mdxt = mean(dxt);
   end

   mdyt = [];
   if ~isempty(dyt)
      mdyt = mean(dyt);
   end

   intv = [ mdxt ; mdyt ];

   intv = sort([ 4*intv ; 2*intv ; intv ; intv/2 ]);

     jj = find( ( intv - 3*bw )  >= 0 );

  
    if isempty(jj)
       intv = [];
    else
       intv = intv(jj(1));
    end

%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function set_border(axe,hb,xl,yl,xt,yt,bw,intv);

%*******************************************************************
% Inner, White Border

    i12 = [ 1  1  2  2 ];
    i21 = [ 1  2  2  1 ];

    y_for_x_0 = yl(i12) + [ -0.5  0    0    0.5 ]*bw;
    x_for_y_0 = xl(i12) + [ -0.5  0    0    0.5 ]*bw;
    y_for_x   = yl(i12) + [ -1   -0.5  0.5  1   ]*bw;
    x_for_y   = xl(i12) + [ -1   -0.5  0.5  1   ]*bw;

    xx_0 = xl([1 2]);
    yy_0 = yl([1 2]);

    i34 = i12(:);
    i34 = i34(:,[1 1]);
    i34(:,2) = i34(:,2) + 2;  % [ 1 1 2 2 ; 3 3 4 4 ]'
    
    y_x_0 = cat( 2 , y_for_x_0(i34) , y_for_x(i34) );
    x_y_0 = cat( 2 , x_for_y_0(i34) , x_for_y(i34) );

    xx_0 = xx_0(i21); xx_0 = xx_0(:);
    xx_0 = xx_0 * ones(1,4);

    yy_0 = yy_0(i21); yy_0 = yy_0(:);
    yy_0 = yy_0 * ones(1,4);


    %--------------------------------------------------------------- 
    % 4 Corners - White

 
       y_for_x_1 = yl(i12) + [-1 0  0  1]*bw;
       x_for_y_1 = xl(i12) + [-1 0  0  1]*bw;

       xx_1 = x_for_y_1( i34([1 3 3 1],i21) );

       yy_1 = y_for_x_1( i34(:,i12) );

       if any(size(xx_1)==1)
       % Matlab4
        xx_1 = reshape(xx_1(:),4,4);
        yy_1 = reshape(yy_1(:),4,4);
       end

%*******************************************************************
% Outher, Black/White Border

if ~isempty(intv)

    nx = prod(size(xt));
    ny = prod(size(yt));

    xx = cat( 2 , xl(1),  ...
     fliplr( xt(1)  : -intv : xl(1)  ) , ...
           ( xt(1)  :  intv : xt(nx) ) , ...
           ( xt(nx) :  intv : xl(2)  ) , xl(2) );

    xx(find(abs(diff(xx))<1e-10))=[];

    
    yy = cat( 2 , yl(1),  ...
     fliplr( yt(1)  : -intv : yl(1)  ) , ...
           ( yt(1)  :  intv : yt(ny) ) , ...
           ( yt(ny) :  intv : yl(2)  ) , yl(2) );

    yy(find(abs(diff(yy))<1e-10))=[];

    nx = size(xx,2);
    ny = size(yy,2);

    xdif = abs(diff(xx([1 2 nx-1 nx])));
    ydif = abs(diff(yy([1 2 ny-1 ny])));

    % first and last segment to small
    xsm = ( xdif([1 3]) < bw );
    ysm = ( ydif([1 3]) < bw );


    ox = ones(1,nx-1);

    y_x = cat( 2 , y_for_x(i34(:,1*ox)) , y_for_x(i34(:,2*ox)) );

    xx = xx([(1:nx-1);(2:nx);(2:nx);(1:nx-1)]);
    xx = permute(xx,([ 1  2 ] + [ 1  -1 ]*(size(xx,2)~=(nx-1))) );
    xx = [ xx  xx ]; 
  
    oy = ones(1,ny-1);

    x_y = cat( 2 , x_for_y(i34(:,1*oy)) , x_for_y(i34(:,2*oy)) );

    yy = yy([(1:ny-1);(2:ny);(2:ny);(1:ny-1)]);
    yy = permute(yy,([ 1  2 ] + [ 1  -1 ]*(size(yy,2)~=(ny-1))) );
       
    yy = [ yy  yy ]; 

    %--------------------------------------------------------------- 
    % Black Patches

    ix = cat( 2 , ( (xsm(1)+1) : 2 : nx-1 ) , ( (xsm(1)+1+nx-1) : 2 : 2*(nx-1) ) );
    iy = cat( 2 , ( (ysm(1)+1) : 2 : ny-1 ) , ( (ysm(1)+1+ny-1) : 2 : 2*(ny-1) ) );

else

    nn = NaN*ones(2,1);

    xx = nn; ix = 1;
    yy = nn; iy = 1;

    x_y = nn;
    y_x = nn;

end

%*******************************************************************
% Set Patches: Corner WhiteX WhiteY BlackX BlackY

  % 4 Corners

  set( hb(1) , 'xdata' , xx_1 , ...
               'ydata' , yy_1 , ...
               'zdata' , [] , ...
               'cdata' , []         );

  % White Stripes

  set( hb(2) , 'xdata' ,  xx_0 , ...
               'ydata' , y_x_0 , ...
               'zdata' , [] , ...
               'cdata' , []         );

  set( hb(3) , 'xdata' , x_y_0 , ...
               'ydata' ,  yy_0 , ...
               'zdata' , [] , ...
               'cdata' , []         );

  % Black Patches

  set( hb(4) , 'xdata' ,  xx(:,ix) , ...
               'ydata' , y_x(:,ix) , ...
               'zdata' , [] , ...
               'cdata' , []         );

  set( hb(5) , 'xdata' , x_y(:,iy) , ...
               'ydata' ,  yy(:,iy) , ...
               'zdata' , [] , ...
               'cdata' , []         );


   %------------------------------------------------
   % Set AxesTickLenght

    tl = get( axe , 'ticklength' );

    tl(1) = 1.5*min(bw./[diff(xl) diff(yl)]);

    set(axe,'ticklength',tl)
    
   %------------------------------------------------
   % Set TitleYPosition

    tt  = get(axe,'title');
    uni = get(tt,'units');
          set(tt,'units','normalized');

    pp    = get(tt,'position');
    pp(2) = 1 + 2*tl(1);
 
    set(tt,'position',pp);

    set(tt,'units',uni)


%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [xt,yt,tick_int,mode] = get_tick(axe,xl,yl,int)


xt = get(axe,'xtick');
yt = get(axe,'ytick');

mode = { get(axe,'xtickmode')  get(axe,'ytickmode') };

tick_int = [];

%*****************************************************************
% Prepare Ticks

if ~any(strcmp(mode,'auto'))
    return
end


if int == 0

   if strcmp(mode{1},'auto')
      xt = [];
   end

   if strcmp(mode{2},'auto')
      yt = [];
   end

   return

end

      
  tick_int = [ 90 ( [ 60 30 20 10  5  2  1 ]/1         ) , ...    % deg
                  ( [    30 20 10  5  2  1 ]/60        ) , ...    % min
                  ( [    30 15 12  6       ]/3600      ) , ...    % min/10
                  ( [    30 15 12  6       ]/3600/10   ) , ...    % min/100
                  ( [    30    12  6       ]/3600/100  )       ]; % min/1000

  lims = cat( 1 , xl , yl );

  dl = diff(lims,1,2);

  nt = size(tick_int,2);

  ii = abs( dl(:,ones(1,nt)) / int - tick_int([1 1],:) );

  [ hilf , ii ] = min(ii,[],2);

  [ hilf , kk ] = max(tick_int(ii));

   ok = find( ( rem( tick_int , tick_int(ii(kk)) ) == 0 )   &  ...
                 (  tick_int >= tick_int(ii(kk))  ) ); 

  [hilf,jj] = min( abs( dl(3-kk) / int - tick_int(ok)) );

  ii(3-kk) = ok(jj);

  tick_int = tick_int(ii(1));

  tic = cell(1,2);

  for ii = 1 : 2

    tic{ii} = tick_int * ( ceil(lims(ii,1)/tick_int) : ...
                          floor(lims(ii,2)/tick_int)       );

    if isempty(tic{ii}) 
       tic{ii} = tick_int * (floor(lims(ii,1)/tick_int) : ...
                              ceil(lims(ii,2)/tick_int)       );
       if size(tic{ii},2) == 2
          tic{ii} = lims(ii,:);
       end     
    end

  end

  if strcmp(mode{1},'auto')
     xt = tic{1};
  end

  if strcmp(mode{2},'auto')
     yt = tic{2};
  end


%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [xtl,ytl] = get_label(xt,yt)

xtl = '';
ytl = '';

%*****************************************************************
% Prepare TickLabels

if any( abs( xt-round(xt) ) > 1e-6 )  |  ...
   any( abs( yt-round(yt) ) > 1e-6 )
   not_flag = 'sexges';  % sexages
else
   not_flag = 'deg';  % degree
end

%----------------------------------------------------
% X

if ~isempty(xt)

   %  -->  [ -180 .. 180 )
   xt0 = xt - 360 * floor( (xt+180) / 360 );

   xt0(end) = xt0(end) + ( 180 - xt0(end) ) * ...
                       ( xt0(end) == -180 );

   xtl = geolabel(xt0,'lon',not_flag);

end
 
%----------------------------------------------------
% Y

if ~isempty(yt)

   ytl = geolabel(yt,'lat',not_flag);

end


%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = geolabel(tick,ini,not);

% GEOLABEL  Transform Geographical Coordinates into String
%
% Label = GEOLABEL( X , INI , NotationFlag )
%
%  X             Vector for Coordinates
%  INI           'lon'      |  'lat'
%  NotationFlag  'sexages'  |  'deg'
%
 
if isempty(tick)
  str = '';
  return
end

lab = [ 'WE' ; 'SN' ];

lab = lab( 1 + strcmp(ini,'lat') , : );

if nargin < 3
  if any( abs( tick-round(tick) ) > 1e-6 )
    not = 'deg';
  else
    not = 'sexages';
  end
end


   grad =   fix(tick);
   minu =   fix(  60*(tick-grad));
   secu = round(3600*(tick-grad-minu/60));

    minu_secu =  fix(secu/60);
   secu = secu - 60*minu_secu;
   minu = minu + minu_secu;

    grad_minu =  fix(minu/60);
   minu = minu - 60*grad_minu;
   grad = grad + grad_minu;

   dec_secu = [];
   sec_form = '';

  if any(diff(1e3*grad+minu)==0)
   dec_secu = round(1e2*( 60*(tick-grad) - minu)) ;
   sec_form = '.%2.2d';
  end
  if ~isempty(dec_secu)
   if any(diff(1e3*grad+minu+1e-3*dec_secu)==0)
    dec_secu = round(1e3*( 60*(tick-grad) - minu)) ;
    sec_form = '.%3.3d';
   end
  end

  is_dec = ~isempty(dec_secu);


  nt = prod(size(tick));

  str = cell(nt,1);
  str(:) = { '' };

  for ii = 1 : nt

     if strcmp(not,'deg')

       str{ii} = [ sprintf('%.0f',abs(grad(ii))) char(176)  ];

     else

       str{ii} = [ sprintf('%.0f',abs(grad(ii))) char(176) , ...
                   sprintf('%2.2d',abs(minu(ii)))           , ...
                   sprintf(sec_form,abs(dec_secu(ii*(1:is_dec)))) char(39)  ];

     end

  end


  for ii = [ 1  nt(1:(nt>1)) ]

    i01 = ( ( tick(ii) > 0 )  |  ( ( ii == 1 )  &  ( tick(ii) == 0 ) ) );
   
    str{ii} = [ str{ii}  lab( 1 + i01 ) ];

  end


