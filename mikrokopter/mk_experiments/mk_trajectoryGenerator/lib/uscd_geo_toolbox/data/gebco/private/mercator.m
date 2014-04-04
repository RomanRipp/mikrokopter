function [hb,ext,ext0,int] = mercator(varargin)

% MERCATOR  MercatorProjection for Axes or Data
%
%---------------------------------------------------------------
% Axes Transformation
%
% MERCATOR( Axe , EXT_New , EXT_Old , NTickInt )
%  
%   Axe       AxesHandle for Projection,        default: GCA
%
%   EXT_New   Scale for new Projection,         default: 1
%   EXT_Old   Scale for old Projection,         default: 0 (no Projection)
%       
%   NTickInt  Number of TickIntervalls on Axis, default: 5
%
%   NOTE: If NTickInt has an NonZero imaginary Part, the AxesLimits
%          will rounded to even parts of the selected TickIntervall.
%
%   EXT_New, EXT_Old, NTickInt must be single finite numerics,
%   EXT_New, EXT_Old with Values between 0 and 1.
%
%
% [ BorderHandle , EXT_New , EXT_Old , NTickInt ] = MERCATOR( Axe , ... ) 
%
%   returns the Handle of the Mercator-Border-Patch
%    which is a Patch for White and one for Black Colors: [ 1 by 2 ]
%
%   The Patches have the tags: 'MERCATOR_BORDER'
%
%   The Black Patches starts at the lower-left Corner.
%   Use a single Character-Argument to start with a White Patch:
%
%    MERCATOR( Axe , ... , C , ... )
%
%   where C is:  'x'  Start White at X-Axis
%                'y'  Start White at Y-Axis
%                'z'  Start White at X-Axis and Y-Axis
%
%  NOTE: A nonzero imaginary Part of EXT_New defines
%         the BorderWidth in Pixel.
%        If the imaginary Part is negative or less then 5
%         No Black Patches are drawn.
%
%---------------------------------------------------------------
% Axes Information about MERCATOR-Settings
%
% [ EXT , NTickInt , AxesHandle , BorderHandle ] = MERCATOR( Axe , 'Info' )
%
% Note: In all cases the first Input "Axe" can be an AxesHandle 
%        or a Handle of an AxesChildren.
%
%---------------------------------------------------------------
% Data Transformation, give a Single-Element-CellArray
%
% 1) YMercator = MERCATOR( {YData} , EXT_New , EXT_Old )
%
%     Transforms the Values Y for MercatorProjection.
%
% 2) YMercator = MERCATOR( Axe , {YData} , EXT_Old )
%
%     Transforms the Values Y for MercatorProjection of Axe, 
%      use EXT_New from Axe.
%
% 3) YMercator = MERCATOR( Axe , EXT_New , {YData} )
%
%     Transforms the Values Y from MercatorProjection of Axe into new, 
%      use EXT_Old from Axe
%
%---------------------------------------------------------------
% Example for usage with AXE_ZOOM:
%
%   figure('toolbar','none')
%   axis([-180 180 -90 90])
%   grid on
%
%   mercator(gca,1,0,10);
%
%   x = linspace(-120,120,100);
%   y = linspace( -60, 60,100);
%
%   plot( x , mercator(gca,{y}) , 'k-' );
%
%   axe_zoom('new',gca,'mercator'); axe_zoom('on',gca);
%
%
%------------------------------------------------------------
%
% MERCATOR changes the  !!!  ydata  !!! of Axes Children
%          
%           and following AxesProperties:
%                            ylim       
%                            ytick
%                            yticklabels
%                            xtick
%                            xticklabels
%                            'dataaspectratio' , [1 1 1]
%
% MANUAL AxesTicks will not be changed, 
%  if the Input "NTickInt" is NOT given !!!
%
% WARNING:  MERCATOR didn't work correctly for TYPE IMAGE
%            Scale YData of the Image to an equal spaced Vector
%            in the new Projection BEFORE.
%
% Example for using IMAGE: 
%
%    xc = ( 1 : 2 : 65 );
%    yc = ( 1 : 2 : 65 )';
%    cc = rand(size(yc,1),size(xc,2));
%
%    ym = deg2merc(yc([1 end]),1);             % Mercator
%    ym = linspace(ym(1),ym(2),size(yc,1))';   % Equal spaced
%    ym = merc2deg(ym,1);                      % Back
%    ym([1 end]) = yc([1 end]);
%
%    cm = interp1(yc,cc,ym);
%
%    % Now use:  image(xc,ym,cm)
%
%**************************************************
% 
% The Calculation for the Projection of the Y-Data:
%
%%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%
% function  y = deg2merc(y,ext)
%
% y = 180/pi * log(abs(tan(pi/4 + y*ext*pi/180/2))) / ext;
%
%
%%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%
% function  y = merc2deg(y,ext)
%
% y = 2*180/pi * ( atan(exp(y*ext*pi/180)) - pi/4 ) / ext; 
%
%%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%
%


Nin  = nargin;
Nout = nargout;

app = mfilename;
fsp = ( app == filesep );
if any(fsp)
    app = app(max(find(fsp))+1:end);
end
app = upper(app);


tag = sprintf('%s_BORDER',app);

 hb = [];

%*****************************************************************
% Check Inputs

zpp = 'AXE_ZOOM';   % To Check for AXE_ZOOM

[msg,axe,ext,ext0,int,pxw,yy,iy,md,inf] = checkin(zpp,varargin{:});

if ~isempty(msg)
    error(msg);
end

if isempty(axe) & isempty(yy)
   return
end

%*****************************************************************

mode = cell(0,2);
if ~isempty(int)
    mode = { 'auto' 'auto' };
end

%*****************************************************************
% Check for YData to transform

% YMercator = MERCATOR( {YData} , Ext_New , Ext_Old )
% YMercator = MERCATOR( Axe , {YData} , Ext_Old )
% YMercator = MERCATOR( Axe , ExtNew , {YData} )

if ~isempty(yy)

   if ( iy == 2 ) & isempty(ext0)
      ext0 = 0;
   end

   if ~isempty(axe)
      [hb,ext,ext0] = check_axe(axe,app,ext,ext0,int,pxw,mode,md);
   end

   if isempty(ext) & ~isempty(axe)

      hb = yy;

   else

      [ext,ext0] = defaults(ext,ext0);

      hb = deg2merc( merc2deg(yy,ext0) , ext );

   end

   return

end
 
%*****************************************************************
% Check for existing Mercator-Patches, if AxesInput

[hb,ext,ext0,int,pxw,mode,md] = check_axe(axe,app,ext,ext0,int,pxw,mode,md);

if inf

   %% [hb,ext,ext0,int] = mercator(varargin)
   %% [ext,int,axe,hb]  = mercator(axe,'info')

   ext0 = int;
    int = hb;
     hb = ext;
    ext = ext0;
   ext0 = axe;

   return

end

[ext,ext0,int,pxw,md] = defaults(ext,ext0,int,pxw,md);

%*****************************************************************

xlim =          get(axe,'xlim');
ylim = merc2deg(get(axe,'ylim'),ext0);

%*****************************************************************
% Fit YLimits
%
% smallest TickInt: 6/3600/100 = 1/60000, see below

sc = 100*3600/6;

yl  = ylim * sc;
ryl = round(yl);

ok = find( abs(yl-ryl) < 1e-3 );

ylim(ok) = ryl(ok)/sc;

if ylim(1) == ylim(2)
   ylim = ylim(1) + [ -1  1 ] * sc/2;
end

%*****************************************************************
% Maximum YAxeIntervalls for MercatorScale

   ymax = [ 84.0   1.0
            89.9   0.6 
            90.0   0.0  ];

   cc =  polyfit(ymax(:,2),ymax(:,1),2);

if ~( ( ext == 0 ) & ( ext0 == 0 ) )

   ylim = limit(cc,ext,ext0,sc,ylim);

end

  ylimM = deg2merc( ylim , ext );


aspect = {};
if ~( ext == 0 ) | strcmp( get(axe,'dataaspectratiomode') , 'auto' )
   aspect = { 'dataaspectratio' , [ 1  1  1 ] };
end

 
set(axe,  aspect{:}     , ...
          'xlim'        , xlim      , ...
          'ylim'        , ylimM     , ...
          'nextplot'    , 'add'                   )

%*****************************************************************
% Get Ticks

[xtic,ytic,tick_int] = get_tick(axe,xlim,ylim,mode,int,ext0);

yticM = deg2merc( ytic , ext );

%*****************************************************************
% Get Border

[bw,upp,intv] = get_int(axe,pxw,xlim,ylimM,xtic,ytic,yticM);

%*****************************************************************
% Check for AutoLimit

acc = 1e-10;
op  = [ -1  1 ];

if ~isequal(imag(int),0) & ~isempty(intv) & ~isempty(tick_int)

    xlim = intv * op .* ceil( op.*(xlim-2*op*acc) / intv );
    ylim = intv * op .* ceil( op.*(ylim-2*op*acc) / intv );

    if ~( ( ext == 0 ) & ( ext0 == 0 ) )
        ylim = limit(cc,ext,ext0,sc,ylim);
    end

    if isequal(mode{1},'auto')
       if isempty(xtic)
          xtic = xlim;
       else
          xt = xtic([1 end]) + op*tick_int;
          ok = ( abs( xt - xlim ) > acc );
          xtic = cat( 2 , xt(1) , xtic , xt(2) );
          xtic = xtic( 1+ok(1) : end-ok(2) );
       end
    end

    if isequal(mode{2},'auto')
       if isempty(ytic)
          ytic = ylim;
       else
          yt = ytic([1 end]) + op*tick_int;
          ok = ( abs( yt - ylim ) > acc );
          ytic = cat( 2 , yt(1) , ytic , yt(2) );
          ytic = ytic( 1+ok(1) : end-ok(2) );
       end
    end

    ylimM = deg2merc( ylim , ext );
    yticM = deg2merc( ytic , ext );

    set(axe,'ylim',ylimM);

    [bw,upp] = get_int(axe,pxw,xlim,ylimM,xtic,ytic,yticM);

end

%*****************************************************************
% TickLabels

[xticl,yticl] = get_label(axe,xtic,ytic);

%*****************************************************************
% Set Title, XYZ_Labe to Normalized

hlb = get( axe , { 'title'  'xlabel'  'ylabel'  'zlabel' } );
hlb = cat(1,hlb{:});
uni = textunit(hlb,'normalized');
 
%*****************************************************************
% Set YData

if ~( ext == ext0 )
    setydata(axe,ext,ext0,cc,sc,tag,zpp);
end

%********************************************************
% Check Ticks at Start and End with AxisLimits

if ~isempty(xtic)
    ii = [ 1  prod(size(xtic)) ];
    jj = ( ( abs( xtic(ii) - xlim ) < acc ) &  ( op.*xtic(ii) > op.*xlim ) );
    if any(jj)
       jj = find(jj);
       xtic(ii(jj)) = xlim(jj);
    end
end

if ~isempty(yticM)
    ii = [ 1  prod(size(yticM)) ];
    jj = ( ( abs( yticM(ii) - ylimM ) < acc ) &  ( op.*yticM(ii) > op.*ylimM ) );
    if any(jj)
       jj = find(jj);
       yticM(ii(jj)) = ylimM(jj);
    end
end

%********************************************************
% Set Axes

set(axe,  'box'         , 'on'      , ...
          'tickdir'     , 'out'     , ...   % !!!!
          'xlim'        , xlim      , ...
          'xtick'       , xtic      , ...
          'xticklabel'  , xticl     , ...
          'ytick'       , yticM     , ...
          'ylim'        , ylimM     , ...
          'yticklabel'  , yticl     , ...
          'nextplot'    , 'add'                   )  

%*****************************************************************
% Set Title, XYZ_Labe back

textunit(hlb,uni);

%********************************************************
% Set Border

if isempty(hb)
   hb = [ patch('parent',axe,'tag',tag,'visible','on') ...
          patch('parent',axe,'tag',tag,'visible','on')     ];
end

border(axe,hb,xlim,ylim,ylimM,xtic,ytic,ext,bw,upp,intv,md);

ud = struct( 'hb'   , { hb   } , ...
             'ext'  , { ext  } , ...
             'int'  , { int  } , ...
             'mode' , { mode } , ...
             'bw'   , { bw   } , ...
             'md'   , { md   } , ...
            'pxw'   , { pxw   } , ...
             'intv' , { intv }       );

setappdata(axe,app,ud);

if Nout == 0
   clear hb
end 

%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function varargout = defaults(varargin);

def = [ 1  0  5  0  0 ];   % Defaults for [ ext  ext0  axeint  pxw  MD ]

nn = min(size(def,2),min(nargin,nargout));

if nn == 0
   return
end

varargout = cell(1,nargout);

for ii = 1 : nn
    if isempty(varargin{ii})
       varargout{ii} = def(ii);
    else
       varargout{ii} = varargin{ii};
    end
end

%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,axe,ext,ext0,int,pxw,yy,iy,md,inf] = checkin(zpp,varargin)

Nin = nargin - 1;

msg    = '';
axe    = [];
ext    = [];
ext0   = [];
int    = [];
yy     = [];
iy     = NaN;

pxw    = [];     % PixelWidth

md     = '';    % Mode

inf    = 0;

is_zoom = NaN;

%*****************************************************************
% Check AxesInput

if Nin == 0
   fig = get(0,'currentfigure');
   if ~isempty(fig)
       axe = get(fig,'currentaxes');
   end
   return
elseif isempty(varargin{1})
   return
end

%*****************************************************************
% Get other Inputs from VARARGIN

val = cell(1,4);  % { axe ext ext0 axeint }

ok = ones(1,5);

n = min(Nin,5);

im = NaN;  % Index for "md"

for ii = 1 : n

    v = varargin{ii};

    ok(ii) = ( prod(size(v)) == 1 );

    jj = ii - ~isnan(im);

    if ok(ii)

       is_cell = iscell(v);
       if is_cell
          v = v{1};
       end

       ok(ii) = isnumeric(v);

       if ok(ii)
          if  is_cell
              ok(ii) = ( isempty(yy) & ( ii <= 3 ) );
              if ok(ii)
                 yy = v;
                 iy = ii;
              end
          else
              ok(ii) = ( isfinite(v) & ( v >= 0 ) );
              if ok(ii)
                 if jj == 1
                    ok(ii) = ishandle(v);
                    if ok(ii)
                       ok(ii) = strcmp(get(v,'type'),'axes');
                       if ~ok(ii)
                           v = get(v,'parent');
                           ok(ii) = strcmp(get(v,'type'),'axes');
                       end
                       if ok(ii)
                          is_zoom = isappdata(v,zpp);
                       end
                    end
                 elseif any( jj == [ 2  3 ] )  % ext | ext0
                    ok(ii) = ( abs(real(v)-0.5) <= 0.5 );
                    if ok(ii) & ( jj == 2 )
                       pxw = imag(v);           % PixelWidth
                    end
                    v = real(v);
                 end
                 if ok(ii)
                    val{jj} = v;
                 end
              end
          end
       else
          ok(ii) = ( isnan(im) & ischar(v) );
          if ok(ii)
             md = lower(v);
             im = ii;
          end
       end

    elseif ~isnan(is_zoom) & ( ii == 2 ) & isequal(size(v),[1 4])

       % AXE_ZOOM | Info

       if ( is_zoom == 1 ) & isnumeric(v)

           % varargin = { axe  lim }

           axe = val{1};

           lim = cat(2,get(axe,'xlim'),get(axe,'ylim'));
           ok(ii) = all( abs(lim-v) <= 1e-10 );

           if ok(ii)
              break
           end

       elseif ischar(v)

           inf = strcmp(lower(v),'info');

       end

       ok(ii) = ( ok(ii) | inf );

       if ok(ii)
          break
       end

    end

end

if ~all(ok(1:n))
    msg = 'Invalid Inputs.';
end

axe  = val{1};
ext  = val{2};
ext0 = val{3};
int  = val{4};

%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [hb,ext,ext0,int,pxw,mode,md] = ...
             check_axe(axe,app,ext,ext0,int,pxw,mode,md);

%*****************************************************************
% Check for existing Mercator-Properties

hb = [];

if isappdata(axe,app);

   ud = getappdata(axe,app);

   hb = ud.hb;

   if ~isempty(hb)
       ok = ishandle(hb);
       if any(~ok)
          delete(hb(find(ok)));
          hb = [];
       end
   end

   if isempty(ext)
      ext = ud.ext;
   end
   if isempty(ext0)
      ext0 = ud.ext;
   end
   if isempty(int)
      int = ud.int;
   end

   if isempty(pxw)
      pxw = ud.pxw;
   end

   if isempty(mode)
      mode = ud.mode;
   end

   if isempty(md)
      md = ud.md;
   end

end

if isempty(mode)
   mode = { get(axe,'xtickmode')  get(axe,'ytickmode') };
end

%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function setydata(axe,ext,ext0,c,sc,tag,zpp);

%***********************************************************************
% Set Y-Data

 shh = get(0,'showhiddenhandles');
       set(0,'showhiddenhandles','on');

 children = get(axe,'children');

       set(0,'showhiddenhandles','off');


 if any(strcmp(get(children,'type'),'image'));
    ww = warnstat;
    if ~strcmp(ww,'off')
        warning('on');
        warning('Projection may not be correct for Data of Type IMAGE.');
        warning(ww);
    end
 end

 for child = children(:)'

   if strcmp(get(child,'type'),'text')

      if strcmp( get(child,'units') , 'data' )

         pos = get(child,'position');
 
         pos(2) = merc2deg(pos(2),ext0);

         pos(2) = deg2merc(pos(2),ext);

         set(child,'position',pos )

      end

   elseif strcmp(get(child,'type'),'light')

          pos = get(child,'position');
 
          pos(2) = merc2deg(pos(2),ext0);

          pos(2) = deg2merc(pos(2),ext);

          set(child,'position',pos )

   elseif ~strcmp(get(child,'tag'),tag)

      %****************************************************************
      % Special for Image (much deficit in Quality !!!)
      %****************************************************************
      if 000 %%% strcmp( get(child,'type') , 'image' ) %%%!!!!!!!!!!!!!
      %****************************************************************

         yy = get(child,'ydata');
         yy = yy(:);

         n  = size(yy,1);

         ok = ( n > 1 );
         if ok
            cc = get(child,'cdata');
            ok = ( ( size(cc,1) == n ) & ~isempty(cc) );
            if ok
               ok = ~all(isnan(cc(:)));
            end
         end

         if ok
          
            y  = yy([1 n]);
            y  = merc2deg(y,ext0);
            ym = deg2merc(y,ext);

            ym = ym(1) + ( ym(2) - ym(1) ) * ( 0 : (n-1) )' / (n-1);

            ym = merc2deg(ym,ext);
            ym([1 n]) = y;
            ym = deg2merc(ym,ext0);
            ym([1 n]) = yy([1 n]);

            ok = ~all( abs(yy-ym) < 1e-10 );

         end

         if ok
 
            is_uint8 = strcmp(class(cc),'uint8');
            if is_uint8
               cc = double(cc);
            end
            for ii = 1 : size(cc,3)
                cc(:,:,ii) = interp1(yy,cc(:,:,ii),ym);
            end
            if is_uint8
               cc = uint8(round(cc));
            end

            ym = merc2deg(ym,ext0);
            ym = deg2merc(ym,ext);

            set(child,'ydata',ym,'cdata',cc);

         else

            y = merc2deg(yy,ext0);

            set( child , 'ydata' , deg2merc(y,ext) );
            

         end

      %****************************************************************
      else
      %****************************************************************

         y = merc2deg(get(child,'ydata'),ext0);

         set( child , 'ydata' ,  deg2merc(y,ext) )

      %****************************************************************
      end
      %****************************************************************

   end
   % Text | Light | ...

 end
 % child

 %******************************************
 % Check Camera

 for p = { 'CameraTarget' 'CameraPosition' }

     if ~strcmp( get(axe,[p{1} 'Mode']) , 'auto' )

         pos = get( axe , p{1} );
 
         pos(2) = merc2deg(pos(2),ext0);

         pos(2) = deg2merc(pos(2),ext);

         set( axe , p{1} , pos );

     end

 end

 %******************************************
 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 % Special for AXE_ZOOM
 
   if isappdata(axe,zpp)
 
      ud = getappdata(axe,zpp);

      hy = ud.History(:,[3 4]);

      hy = merc2deg(hy,ext0);

      hy = limit(c,ext,ext0,sc,hy);

      ud.History(:,[3 4]) = deg2merc(hy,ext);
            
      setappdata( axe , zpp , ud );

    end

 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 %******************************************
 
%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function y = limit(c,ext,ext0,sc,y)


   op  = [ -1   1 ];  % Operator

   yl  = op * polyval(c,ext);
   yl0 = op * polyval(c,ext0);

    n  = ones(size(y,1),1);

   % Set Old Values near yl0 to yl !!!

   ok = ( abs( n*yl0 - y ) <  1e-3/sc );
    y = y.*(~ok) + (n*yl).*ok;

   % Set Values to valid Range

   ok = ( ( n*(yl.*op) - y.*(n*op) ) > 0 );
    y = y.*ok + (n*yl).*(~ok);

%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [bw,upp,intv] = get_int(axe,pxw,xlim,ylimM,xtic,ytic,yticM);

bw   = [];
intv = [];

pm   = 5;   % Minimum 5 Pixel for Black/White

%********************************************************
% Get BorderWidth

   [upp,axepos,figpos] = ppunit(axe);

   if ~( pxw == 0 )

       bw = abs(pxw) + i * ( pxw < pm );

   else

      bw = 0.010 * max(axepos([3 4]));
      bw = min( bw , 0.005 * figpos(4) );
   
      fig = get(axe,'parent');
   
      papuni = get(fig,'paperunits');
      set(fig,'paperunits','inches');
      pappos = get(fig,'paperposition');
      set(fig,'paperunits',papuni);

      % Minimum 11 Points == 1.86 mm (150 dpi)
      bmin = 11/150 / (pappos(4)/figpos(4));
   
      bw = max( bw , max( pm , bmin ) );

   end

if nargout < 3
   return
end

%********************************************************
% Get Intervall

   jj = find( ( xtic < xlim(1) ) | ( xtic > xlim(2) ) );
   xtic(jj) = [];

   jj = find( ( yticM < ylimM(1) )  | ( yticM > ylimM(2) ) );
   ytic(jj) = [];

   dxt = diff(xtic);
   dyt = diff(ytic);

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

   intv = [ mdxt/upp(1) ; mdyt/upp(2) ];
   upi  = [ 1+0*mdxt    ; 2+0*mdyt    ];

   intv = [ 4*intv ; 2*intv ; intv ; intv/2 ];
   upi  = [   upi  ;   upi  ;  upi ;  upi   ];
 
   [intv,si] = sort(intv);

    upi = upi(si);

     jj = find( ( intv - 3*bw )  >= 0 );

    if isempty(jj)
       intv = [];
    else
         jj = jj(1);
       intv = intv(jj) * upp(upi(jj));
    end

%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function border(axe,hb,xlim,ylim,ylimM,xtic,ytic,ext,bw,upp,intv,md);


nbl = ~( imag(bw) == 0 );  % No Black/White

bw = real(bw);

%*******************************************************************
% Inner, White Border

    i12 = [ 1  1  2  2 ];
    i21 = [ 1  2  2  1 ];

    hlf = ( 1 + nbl ) / 2;  %  Half: 0.5 | 1 
    nlf = ( 1 - hlf );

    y_for_x_0 = ylimM(i12) + [ -hlf  0    0    hlf ] * bw * upp(2);
    x_for_y_0 =  xlim(i12) + [ -hlf  0    0    hlf ] * bw * upp(1);
    y_for_x   = ylimM(i12) + [ -1   -nlf  nlf  1   ] * bw * upp(2);
    x_for_y   =  xlim(i12) + [ -1   -nlf  nlf  1   ] * bw * upp(1);

    xx_0 = xlim([1 2]);
    yy_0 = ylimM([1 2]);

    i34 = i12(:);
    i34 = i34(:,[1 1]);
    i34(:,2) = i34(:,2) + 2;  % [ 1 1 2 2 ; 3 3 4 4 ]'
   
    y_x_0 = cat( 2 , y_for_x_0(i34) , y_for_x(i34) );
    x_y_0 = cat( 2 , x_for_y_0(i34) , x_for_y(i34) );

    xx_0 = xx_0(i21); xx_0 = xx_0(:);
    xx_0 = xx_0 * ones(1,4);

    yy_0 = yy_0(i21); yy_0 = yy_0(:);
    yy_0 = yy_0 * ones(1,4);

    xw = [ xx_0 x_y_0];
    yw = [ y_x_0 yy_0];

    %--------------------------------------------------------------- 
    % 4 Corners - White

 
    y_for_x_1 = ylimM(i12) + [ -1  0  0  1 ] * bw * upp(2);
    x_for_y_1 = xlim(i12)  + [ -1  0  0  1 ] * bw * upp(1);

    xx_1 = x_for_y_1( i34([1 3 3 1],i21) );

    yy_1 = y_for_x_1( i34(:,i12) );

    if any(size(xx_1)==1)
    % Matlab4
      xx_1 = reshape(xx_1(:),4,4);
      yy_1 = reshape(yy_1(:),4,4);
    end

    xw = [ xx_1  xw ];
    yw = [ yy_1  yw ];

%*******************************************************************
% Outher, Black/White Border

if isempty(intv) | nbl

    xk = NaN;
    yk = NaN;

else

    nx = prod(size(xtic));
    ny = prod(size(ytic));

    xx = cat( 2 , xlim(1),  ...
     fliplr( xtic(1)  : -intv : xlim(1)  ) , ...
           ( xtic(1)  :  intv : xtic(nx) ) , ...
           ( xtic(nx) :  intv : xlim(2)  ) , xlim(2) );

    xx(find(abs(diff(xx))<1e-10))=[];

    
    yy = cat( 2 , ylim(1),  ...
     fliplr( ytic(1)  : -intv : ylim(1)  ) , ...
           ( ytic(1)  :  intv : ytic(ny) ) , ...
           ( ytic(ny) :  intv : ylim(2)  ) , ylim(2) );

    yy(find(abs(diff(yy))<1e-10))=[];

    yy = deg2merc(yy,ext);


    nx = size(xx,2);
    ny = size(yy,2);

    xdif = abs(diff(xx([1 2 nx-1 nx])));
    ydif = abs(diff(yy([1 2 ny-1 ny])));

    % first and last segment to small
    xsm = ( xdif([1 3]) < bw );
    ysm = ( ydif([1 3]) < bw );

    %---------------------------------------------------------

    ox = ones(1,nx-1);

    y_x = cat( 2 , y_for_x(i34(:,1*ox)) , y_for_x(i34(:,2*ox)) );

    if ~( size(y_x,1) == 4 )
        y_x = reshape(y_x,4,2*(nx-1));
    end

    xx = xx([(1:nx-1);(2:nx);(2:nx);(1:nx-1)]);
    xx = permute(xx,([ 1  2 ] + [ 1  -1 ]*(size(xx,2)~=(nx-1))) );
    xx = [ xx  xx ]; 
  
    %---------------------------------------------------------

    oy = ones(1,ny-1);

    x_y = cat( 2 , x_for_y(i34(:,1*oy)) , x_for_y(i34(:,2*oy)) );

    if ~( size(x_y,1) == 4 )
        x_y = reshape(x_y,4,2*(ny-1));
    end

    yy = yy([(1:ny-1);(2:ny);(2:ny);(1:ny-1)]);
    yy = permute(yy,([ 1  2 ] + [ 1  -1 ]*(size(yy,2)~=(ny-1))) );
       
    yy = [ yy  yy ]; 

    %--------------------------------------------------------------- 
    % Black Patches

    ix = xsm(1) + 1 + any( md == 'xz' );
    iy = ysm(1) + 1 + any( md == 'yz' );

    ix = ( ix : 2 : nx-1 );
    iy = ( iy : 2 : ny-1 );

    ix = cat( 2 , ix , ix+(nx-1) );
    iy = cat( 2 , iy , iy+(ny-1) );


    xk = [ xx(:,ix)  x_y(:,iy) ];
    yk = [ y_x(:,ix)  yy(:,iy) ];

end

%*******************************************************************
% Set Patches

    % 4 Corners and 2 White Stripes on each Side

  set( hb(1) , 'xdata' , xw , ...
               'ydata' , yw , ...
               'zdata' , [] , ...
               'cdata' , [] , ... 
           'facecolor' , [1 1 1]              , ...
           'edgecolor' , [0 0 0]              , ...
           'linewidth' , get(axe,'linewidth') , ...
            'clipping' , 'off'                , ...
           'erasemode' , 'background'               );

    % Black Patches

  set( hb(2) , 'xdata' , xk , ...
               'ydata' , yk , ...
               'zdata' , [] , ...
               'cdata' , [] , ... 
           'facecolor' , [0 0 0]              , ...
           'edgecolor' , [0 0 0]              , ...
           'linewidth' , get(axe,'linewidth') , ...
            'clipping' , 'off'                , ...
           'erasemode' , 'background'               )


   %------------------------------------------------
   % Set AxesTickLenght

    bw = ( bw ./ [diff(xlim)/upp(1) diff(ylimM)/upp(2)] );  % normalized

    td = get( axe , 'tickdir' );

    out = strcmp( td , 'out' );

    bw = ( 1 + 0.5*out ) * bw;

    if strcmp( td , 'out' )

       set(axe,'ticklength',min(bw) * [ 1  1 ])
    
    end

   %------------------------------------------------
   % Set TitleYPosition

    tt  = get(axe,'title');

    uni = textunit(tt,'normalized');

    pp  = get(tt,'position');

    if ( pp(2) > 1 ) & strcmp(get(tt,'verticalalignment'),'bottom')

       pp(2) = 1 + 1.5*bw(2);
 
       set( tt , 'position' , pp );

       set( tt , 'units' , uni )

    end

%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  y = deg2merc(y,ext)

if ext ~= 0

   y = 180/pi * log(abs(tan(pi/4 + y*ext*pi/180/2))) / ext;

end


%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  y = merc2deg(y,ext)
 
if ext ~= 0

   y = 2*180/pi * ( atan(exp(y*ext*pi/180)) - pi/4 ) / ext; 

end


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

function [xtic,ytic,tick_int] = get_tick(axe,xlim,ylim,mode,int,ext0)

acc = 1e-10;

xtic =          get(axe,'xtick');
ytic = merc2deg(get(axe,'ytick'),ext0);

int = real(int);

tick_int = [];

if ~isempty(xtic)
    xtic = xtic( find( ( xlim(1)-acc < xtic ) & ( xtic < xlim(2)+acc ) ) );
end

if ~isempty(ytic)
    ytic = ytic( find( ( ylim(1)-acc < ytic ) & ( ytic < ylim(2)+acc ) ) );
end


%*****************************************************************
% Prepare Ticks

if any(strcmp(mode,'auto'))

  tick_int = [ 90 ( [ 60 30 20 10  5  2  1 ]/1         ) , ...    % deg
                  ( [    30 20 10  5  2  1 ]/60        ) , ...    % min
                  ( [    30 15 12  6       ]/3600      ) , ...    % min/10
                  ( [    30 15 12  6       ]/3600/10   ) , ...    % min/100
                  ( [    30    12  6       ]/3600/100  )       ]; % min/1000

  lims = cat( 1 , xlim , ylim );

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

  limt = lims/tick_int;
  limr = round(limt);
    jj = ( abs(limt-limr) < acc );
  if any(jj)
     jj = find(jj);
     limt(jj) = limr(jj);
  end

  for ii = 1 : 2

    tic{ii} = tick_int * ( ceil(limt(ii,1)) : floor(limt(ii,2)) );

    if isempty(tic{ii}) 
       tic{ii} = tick_int * (floor(limt(ii,1)) : ceil(limt(ii,2)) );
       if size(tic{ii},2) == 2
          tic{ii} = lims(ii,:);
       end     
    end

  end

  if strcmp(mode{1},'auto')
     xtic = tic{1};
  end

  if strcmp(mode{2},'auto')
     ytic = tic{2};
  end

end

%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [xticl,yticl] = get_label(axe,xtic,ytic)

xticl = '';
yticl = '';

%*****************************************************************
% Prepare TickLabels

if any( abs( xtic-round(xtic) ) > 1e-6 )  |  ...
   any( abs( ytic-round(ytic) ) > 1e-6 )
   flag = 'sexges';  % sexages
else
   flag = 'deg';  % degree
end

%----------------------------------------------------
% X

if ~isempty(xtic)

   xticl = geolabel(axe,xtic,'lon',flag);

end
 
%----------------------------------------------------
% Y

if ~isempty(ytic)

   yticl = geolabel(axe,ytic,'lat',flag);

end


%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = geolabel(axe,tick,ini,flag);

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

isy = strcmp(ini,'lat');

lab = [ 'WE' ; 'SN' ];

lab = lab( 1 + isy , : );

pre = ( isy & strcmp(get(axe,'yaxislocation'),'right') );
pre = ( ~isy | pre );

pre = char( 32 * ones(1,pre) );

app = ' ';

if nargin < 4
  if any( abs( tick-round(tick) ) > 1e-6 )
    flag = 'deg';
  else
    flag = 'sexages';
  end
end

if isy

   while 1
         ii = ( abs(tick) > 90 );
         if ~any(ii)
             break
         end
         ii = find(ii);
         tick(ii) = tick(ii) + ( sign(tick(ii)) * 180 - 2*tick(ii) );
   end

else

   %  -->  [ -180 .. 180 )

   tick = tick - 360 * floor( (tick+180) / 360 );

   tick(end) = tick(end) + ( 180 - tick(end) ) * ...
                           ( tick(end) == -180 );

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

     if strcmp(flag,'deg')

       str{ii} = [ sprintf('%.0f',abs(grad(ii))) char(176)  ];

     else

       if is_dec
          str{ii} = sprintf(sec_form,abs(dec_secu(ii)));
       end

       str{ii} = sprintf( '%.0f%s%2.2d%s%s' , abs(grad(ii)) , char(176) , ...
                           abs(minu(ii)) , str{ii} , char(39) );
 
     end

     if any( ii == [ 1  nt(1:double(nt>1)) ] )
        i01 = ( ( tick(ii) > 0 )  |  ( ( ii == 1 )  &  ( tick(ii) == 0 ) ) );
        str{ii} = [ str{ii}  lab( 1 + double(i01) ) ];
     end
         
     str{ii} = sprintf('%s%s%s',pre,str{ii},app);

  end

%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function org = textunit(h,uni)

% TEXTUNIT   Set Property UNITS of TextObjects
%
% TEXTUNIT preserve the Problem of PositionShift
%  while change from 'data' to normalized Units.
%
% Note: During set to Units 'pixels' a PositionShift
%       may be occur cause of round of PositionValues.
%
% OrgUnit = TEXTUNIT( TextHandle , NewUnit  )
% OrgUnit = TEXTUNIT( TextHandle , NewUnits )
%
% NewUnit(s) can be a CharacterArray or CellArray of Strings.
%
% To get the possible Values of the TextProperty Unit type:
%
%    >> display(set(0,'DefaultTextUnits'))
%
% The Number of Elements of TextHandles and NewUnits 
%  must be equal.
%

Nin  = nargin;
Nout = nargout;

org = cell(0,1);

if Nin < 1
   if Nout == 0, clear org, end
   return
end

%-----------------------------------------------------
% Check Inputs

msg = cell(0,1);

if ~isempty(h)

    h = h(:);
    n = size(h,1);

    ok = isnumeric(h);
    if ok
       ok = all(ishandle(h));
       if ok
          ok = all(strcmp(get(h,'type'),'text'));
       end
    end

    if ~ok
        msg = cat( 1 , msg , {'Invalid TextHandles.'} )
    else
        org = cell(n,1);
        for ii = 1 : n
            org{ii} = get(h(ii),'units');
        end
    end
end

if Nin == 2
   [ok,uni] = chkcstr(uni);
   if ok
      uni = uni(:);
      ok = ~any(strcmp(uni,''));
   end
   if ~ok
       msg = cat( 1 , msg , {'Units must be String(s).'} );
   elseif ~isempty(h)
       def = set(0,'DefaultTextUnits');
       for ii = 1 : size(uni,1)
           if isempty(strmatch(uni{ii},def))
               m = sprintf(' ''%s'' |',def{:});
               msg = cat( 1 , msg , {sprintf('Units must be any of:%s',m(1:(end-1)))} );
               break
           end
       end 
       if ~any( size(uni,1) == [ 1  n ] )
           msg = cat( 1 , msg , {'UnitNumber must be equal to HandleNumber.'} );
       elseif ( size(uni,1) == 1 ) & isempty(msg) 
              uni = uni(ones(1,n));
       end    
   end
end

if ~isempty(msg)
    error(sprintf('%s\n',msg{:}));
end

if ( Nin < 2 ) | isempty(h)
   if n == 1, org = org{1}; end
   if Nout == 0, clear org, end
   return
end

%-----------------------------------------------------
% Set Units

c = 'xyz';
m = size(c,2);

for ii = 1 : n

    
    if strcmp(org{ii},'data')

       % 'data' --> 'norm'

       a = get( h(ii) , 'parent' );

       p = get( h(ii) , 'position' );

       for jj = 1 : m

           l = get( a , [ c(jj)  'lim' ] );

           if strcmp( get( a , [ c(jj)  'dir' ] ) , 'reverse' )
              l = l([2 1]);
           end

           p(jj) = ( p(jj) - l(1) ) / ( l(2) -l(1) );

       end 

       set( h(ii) , 'units' , 'normalized' , ...
                 'position' , p                    )

    end

    set( h(ii) , 'units' , uni{ii} );

end

if n == 1, org = org{1}; end

if Nout == 0, clear org, end

%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,str] = chkcstr(str,opt)


% CHKCSTR  Checks Input for CellString, contains Strings !
%
%  [ok,str] = chkcstr(str,Option)
%
%  Option ~= 0 ==> CharacterArrays not allowed,
%
%   default: Option == 0   ==>  CharacterArrays --> CellString
%
 
if nargin < 2
   opt = 0;
end

if strcmp(class(str),'char') & isequal(opt,0)
   str = cellstr(str);
end

ok = iscellstr(str);
if ~ok
   return
end

try
  s = cat(2,str{:});
catch
  ok = 0;
  return
end
 
ok = ( strcmp(class(s),'char')  &  ( prod(size(s)) == size(s,2) ) );


%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ww = warnstat

% WARNSTAT  Returns global WarningStatus
%
%  WARNSTAT returns the Status of WARNING
%
% Matlab R<13   WARNING
% Matlab R>12   WARNING for Identifier ALL
%

ww = warning;

if isstruct(ww)   % New Matlab R>12 Syntax
   try
      id = strcmp({ww.identifier},'all');
      if any(id)
         id = find(id);
         ww = ww(id(1)).state;
      else
         ww = '';
      end
   catch
      ww = '';
   end
elseif ~chkstr(ww)
   ww = '';
end

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = chkstr(str,opt)


% CHKSTR  Checks Input for String
%
%  ok = chkstr(str,Option)
%
%  Option ~= 0 ==> only true for nonempty Strings
%
%   default:   Option == 0  (true for empty Strings)
%

 
if nargin < 2
   opt = 0;
end

ok = ( strcmp( class(str) , 'char' )      & ...
       ( prod(size(str)) == size(str,2) ) & ...
       ( isequal(opt,0) | ~isempty(str) )         );

