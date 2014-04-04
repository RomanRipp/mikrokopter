function fig = platdemo(varargin)

% PLATDEMO   Run Demonstration for PLATEAU
%
% Displays a 3D-Surface of a Plateau
%
% The Plateau- and ViewParameter can be controled by a GUI.
%
% see also: LIGHTANG
%

fig = [];

Nin  = nargin;
Nout = nargout;

if Nin == 0
   fcn = mfilename;
   fsp = ( fcn == filesep );
   if any(fsp)
      fcn = fcn(max(find(fsp))+1:end);
   end
   fig = plat_new(fcn);
   if Nout == 0
      clear fig
   end
   return
end


fcn = varargin{1};

if ~( ischar(fcn) & ~isempty(fcn) & ...
      ( prod(size(fcn)) == size(fcn,2) ) )
    error('First Input must be a String (SubFunction).');
end


msg = '';

try
   feval(fcn,varargin{2:end});
catch 
   msg = lasterr;
end

if ~isempty(msg)
    error(sprintf('Error call "%s".\n%s\n',fcn,msg));
end

if Nout == 0
   clear fig
end

%***************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function fig = plat_new(fcn)


nc = 100;

ResizeFcn = sprintf('%s(''plat_pos'',gcbo);',fcn);

fig = figure( 'colormap' , rgb(nc)    , ...
                 'color' , 'w'    , ...
               'menubar' , 'none' , ...
               'toolbar' , 'none' , ...
           'numbertitle' , 'off' , ...
                  'name' , 'Demo PLATEAU' , ...
             'ResizeFcn' , ResizeFcn            );

axe = axes( 'parent' , fig       , ...
               'box' , 'off'     , ...
              'xlim' , [ -1  1 ] , ...
              'ylim' , [ -1  1 ] , ...
              'zlim' , [  0  1 ] , ...
              'clim' , [  0  1 ] , ...
             'color' , 'none'    , ...
             'layer' , 'top'     , ...
           'tickdir' , 'out'     , ...
          'userdata' , []        , ...    % YTicks
          'nextplot' , 'add'     , ...
        'projection' , 'orthographic' , ...
              'tag'  , 'Axes'            );

srf = surface( 'parent' , axe      , ...
                'xdata' , []       , ...
                'ydata' , []       , ...
                'zdata' , []       , ...
                'cdata' , []       , ...
            'edgecolor' , 'none'   , ...
            'facecolor' , 'interp' , ...
         'cdatamapping' , 'scaled' , ...
                  'tag' , 'Surface'      );

if 1
set( srf , ...
         'facelighting' , 'gouraud' , ...
         'AmbientStrength'  , 0.4  , ...
         'DiffuseStrength'  , 0.6  , ...
         'SpecularStrength' , 0.2  , ...
         'SpecularExponent' , 10   , ...
  'SpecularColorReflectance', 1          ); 
end

lgt = light( 'parent' , axe , ...
             'style' , 'infinite' , ...
             'color' , [ 1  1  1 ] , ...
             'tag'   , 'Light'            );

%------------------------------------------------------
% ColorBar

nc = size( get(fig,'colormap') , 1 );

bar = axes( 'parent' , fig       , ...
               'box' , 'on'      , ...
              'xlim' , [ 0  1 ]  , ...
              'ylim' , [ 0  1 ]  , ... %  1 nc ] + 0.5* [ -1  1 ] , ...
              'clim' , [ 0  1 ]  , ...
             'xtick' , []        , ...
     'yaxislocation' , 'right'   , ...
             'color' , 'none'    , ...
             'layer' , 'top'     , ...
          'nextplot' , 'add'     , ...
              'tag'  , 'ColorBar'            );

img = image( 'parent' , bar        , ...
             'xdata' , 0.5         , ...
             'ydata' , [ 0  1 ]' + 0.5/nc * [ 1  -1 ]' , ...
             'cdata' , ( 1 : nc )' , ...
      'cdatamapping' , 'direct'    , ...
              'tag'  , 'ColorImage'            );

%------------------------------------------------------
% Controls

ini = plat_ini;

hc = ini(:,[1 3]);
  
for ii = ( 1 : size(ini,1) )

    if iscellstr(ini{ii,2})
       style = 'popupmenu';
    else
       style = 'edit';
    end

    CB = sprintf('%s(''plat_clb'',gcbo);',fcn);

    h = uicontrol( 'parent' ,  fig     , ...
                    'style' ,  style   , ...
                      'min' ,  0       , ...
                      'max' ,  1       , ...
                 'callback' ,  CB      , ...
               'fontweight' , 'normal' , ...
          'foregroundcolor' , 'k'      , ...
          'backgroundcolor' , [0.9 0.95 1.0] , ... 
                      'tag' , ini{ii,1}           );

    if strcmp(style,'edit')

       frm = sprintf('%%.%.0ff',ini{ii,4});  % Format for  Edit <---> Slider

       str = sprintf(frm,ini{ii,3});
       val = eval(str);

       mn = eval(sprintf(frm,ini{ii,2}(1)));
       mx = eval(sprintf(frm,ini{ii,2}(2)));

       s = uicontrol( 'parent' ,  fig          , ...
                       'style' , 'slider'      , ...
                         'min' ,  mn           , ...
                         'max' ,  mx           , ...
                       'value' ,  val          , ...
                    'callback' ,  CB           , ...
                    'userdata' ,  frm          , ...
             'foregroundcolor' , 'k'           , ...
             'backgroundcolor' , [0.8 0.9 1.0] , ... 
                         'tag' ,  ini{ii,1}             );

       set( h , 'string' , str , ...
              'userdata' , str , ...
   'horizontalalignment' , 'center'      );

       h = h + i*s;

    else

       set( h , 'string' , ini{ii,2} , ...
                 'value' , ini{ii,3} , ...
   'horizontalalignment' , 'left' );

    end

    hc{ii,2} = {h};

    %----------------------------------------------
    % TextLabel above

    uicontrol( 'parent' ,  fig       , ...
                'style' , 'text'     , ...
               'string' ,  ini{ii,1} , ...
                  'min' ,  0         , ...
                  'max' ,  1         , ...
             'callback' ,  ''        , ...
           'fontweight' , 'bold'     , ...
      'foregroundcolor' , 'k'        , ...
      'backgroundcolor' ,  get(fig,'color') , ... 
                  'tag' ,  ini{ii,1}           );


end

%------------------------------------------------------
% Initialisation

hc = permute(hc,[2 1]);
hc = struct(hc{:});

plat_clb( imag(hc.Length) );

plat_clb( hc.Azimut );
plat_clb( hc.Elevation );

plat_clb( hc.ColorMap );
plat_clb( hc.ZScale );

plat_pos( fig , 1 );

plat_clb( hc.LightAzim );
plat_clb( hc.LightElev );

%***************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ini = plat_ini

% WIN_INI  UIControl-Initialisation

mode = { 'triangle'
         'cosine'
         'gauss'
         'exp'
         'potenz'   };

mval = 2;

len  = [ 11  301 ];

azm   = [ -180  180 ];
elv   = [    0   90 ];

laz   = azm;
lel   = [ -1  90 ];

rad   = [ 0.01  2.0  ];   
rot   = [ -90   90   ];
off   = [ 0.0   0.99 ];
scl   = [ 0.2   2.0  ];


clm = { 'bone'
        'copper'
        'gray'
        'hot'
        'cool'
        'hsv'
        'jet'
        'rgb'
        'spring'
        'summer'
        'autumn'
        'winter'  };

%         Tag     String/MinMax   Value  NK
ini   = { 'Mode'       mode     mval     0
          'Length'     len      101      0
          'RadiusX'    rad      0.8      2
          'RadiusY'    rad      0.8      2
          'Rotation'   rot      0.0      0
          'Offset'     off      0.0      2
          'ColorMap'   clm      1.0      0
          'ZScale'     scl      0.5      1
          'Azimut'     azm       20      0
          'Elevation'  elv       30      0
          'LightAzim'  laz       60      0
          'LightElev'  lel       40      0    };

%***************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function plat_pos(fig,mode)


figoff = [ 20 30 20 40 ];  % pixels: [ Left Bottom Right Top ]
axeoff = [  6  5  4  3 ];  % characters: [ Left Bottom Right Top ]
baroff = [  2  5  ];       % characters: [ Width Right ] 

axemin = [ 100 100 ];      % pixels: [ Width Height ] (minimum)
  
%--------------------------------------------------------------
% FigurePosition

[fs,ssi,ppi,is_win] = bestfont;

   ssi = ssi([3 4]);

if nargin > 1
   axemin = [ 2/3 1/2 ] .* ssi;
end

figuni = get(fig,'units');
         set(fig,'units','pixels');

figpos = get(fig,'position');

figsiz = ssi - [ sum(figoff([1 3]))  sum(figoff([2 4])) ];  % max

%--------------------------------------------------------------
% ControlExtension

ini = plat_ini;
 ni = size(ini,1);

nv = 2;                      % Number Vertical
nh = ceil(ni/nv);            % Number Horizontal


fh = fs/72 * ppi;

hs = 10+is_win;              %  SliderHeight
bw = 6;                      % Double Width of Border

hc = 2;                      % CharacterHeight of Control

hh = ceil( hc * fh );        % ControlHeight
ww = ceil(  8 * fh );        % ControlWidth
ww = ww + hh;

ht = hh - bw;                % TextHeight

hd = ceil( 2.0 * hh );       % Horizontal Distance
vd = ceil( 0.5 * hh );       % Vertical   Distance

h0 = hs + hh + ht;           % Height of ControlRow

hm = nv * ( h0 + vd );       % Vertical   Extension of ControlField
wm = nh * ( ww + hd ) + hd;  % Horizontal Extension of ControlField

%--------------------------------------------------------------
% AxeExtension

axesum = [ sum(axeoff([1 3]))  sum(axeoff([2 4])) ];
axesum = axesum + [ sum(baroff)  0 ];

axesiz = fh * axesum + axemin;

figorg = figpos;

figpos([3 4]) = max( figpos([3 4]) , axesiz + [ 0  hm ] );
figpos([3 4]) = min( figpos([3 4]) , figsiz );


if axesiz(2)+hm > figpos(4)

   off = ( axemin(2) + nv * ( hs - bw ) );

    nn = ( axesiz(2)+hm - off ) / hh;

    hh = ceil( ( figpos(4) - off ) / nn );

    hd = ceil( 2.0 * hh );  % Horizontal Distance
    vd = ceil( 0.5 * hh );  % Vertical   Distance

    ht = hh - bw;           % TextHeight

    h0 = hs + hh + ht;      % Extension of a Row

    hm = nv * ( h0 + vd );
    wm = nh * ( ww + hd ) + hd;

end
 
%--------------------------------------------------------------
% FigureWidth

if ( wm > figsiz(1) )
   hd = ( figsiz(1) - nh*ww ) / ( nh + 1 );
   hd = max( floor(hd) , ceil(0.5*hh) );
   ww = floor( ( figsiz(1) - (nh+1) * hd ) / nh );
   wm = nh * ( ww + hd ) + hd;
end

if wm < figpos(3)
   hd = ( figpos(3) - nh*ww ) / ( nh + 1 );
   hd = floor(hd);
else
   figpos(3) = wm;
end

%--------------------------------------------------------------
% FigurePosition,  Fit to upper Left

figpos(2)     = sum(figorg([2 4])) - figpos(4);

off = figpos([1 2]) + figpos([3 4]) - ( ssi - figoff([3 4]));

figpos([1 2]) = figpos([1 2]) - max(0,off);

figpos([1 2]) = max( figpos([1 2]) , figoff([1 2]) );

%--------------------------------------------------------------
% AxesPosition

axeoff = axeoff * hh/hc;
baroff = baroff * hh/hc;

axepos        = zeros(1,4);
axepos([1 2]) = axeoff([1 2]) + [ 0  hm ];
axepos([3 4]) = axepos([1 2]) + axeoff([3 4]) + [ sum(baroff) 0 ];
axepos([3 4]) = figpos([3 4]) - axepos([3 4]);

barpos    = axepos;
barpos(1) = axepos(1) + axepos(3) + axeoff(3);
barpos(3) = baroff(1);

axepos = axepos ./ figpos([3 4 3 4]);
barpos = barpos ./ figpos([3 4 3 4]);

%--------------------------------------------------------------

axe = findobj(fig,'type','axes','tag','Axes');
bar = findobj(fig,'type','axes','tag','ColorBar');

if any( abs(figpos-figorg) >= 1 )

   fcn = get(fig,'ResizeFcn');
         set(fig,'ResizeFcn','');

   set( fig , 'units'    , 'pixels'     , ...
              'position' ,  figpos            );

   wygiwys(fig);

   set( fig , 'ResizeFcn' , fcn );

end

set( axe , 'units'    , 'normalized' , ...
           'position' ,  axepos            );

set( bar , 'units'    , 'normalized' , ...
           'position' ,  barpos            );

%--------------------------------------------------------------
% XLabel of Bar

xl = get(bar,'xlabel');

set( xl ,     'units' , 'normalized' , ...
           'position' ,  [ 1  -fh/figpos(4)/axepos(4) 0 ] , ...
'horizontalalignment' , 'center'     , ...
  'verticalalignment' , 'top'            );

%--------------------------------------------------------------
% Title of Axes

tl = get(axe,'title');

set( tl ,     'units' , 'normalized' , ...
           'position' ,  [ 0.5  1+fh/figpos(4)/axepos(4) 0 ] , ...
'horizontalalignment' , 'center'     , ...
  'verticalalignment' , 'bottom'            );

%--------------------------------------------------------------
% Set Controls

for ii = 1 : ni;

    hc = findobj(fig,'type','uicontrol','tag',ini{ii,1});

    if ~isempty(hc)

        x0 = ( ceil(ii/nv) - 1 ) * ( ww + hd ) + hd;

        y0 = mod(ii,nv);
        y0 = y0 + nv * ( y0 == 0 );
        y0 = ( nv - y0 ) * ( h0 + vd ) + vd;

        for h = hc(:)'

            style  = get(h,'style');

            no_sld = ~strcmp(style,'slider');
            is_txt =  strcmp(style,'text');

            y1 = y0 + hs*no_sld + hh*is_txt;
            h1 = hs + ( hh - hs ) * no_sld;
            h1 = h1 + ( ht - h1 ) * is_txt;

            pos = [ x0  y1  ww  h1 ];  

            set( h , 'units' , 'pixels' , ...
                  'position' ,  pos     , ...
                'fontunits'  , 'points' , ...
                'fontsize'   ,  fs             );
 
       end

    end

end

%***************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function plat_clb(h);

% WIN_CLB  Callbacks

ok = ( isnumeric(h) & ( prod(size(h)) == 1 ) );
if ok
   ok = ishandle(h);
   if ok
      ok = strcmp(get(h,'type'),'uicontrol');
   end
end

if ~ok
    warning('UIControlHandle required.');
    return
end

fig = get(h,'parent');
axe = findobj(fig,'type','axes','tag','Axes');

tag = get(h,'tag');
val = get(h,'value');
str = get(h,'string');

msg = '';

%****************************************************************
switch get(h,'style')
%****************************************************************

   %-------------------------------------------------------------
   case 'edit'
   %-------------------------------------------------------------

         hs = findobj(fig,'style','slider','tag',tag);

         if ~isempty(hs)

             frm = get( hs , 'userdata' );

             try
                val = eval(str);
             catch
                val = [];
             end

             if ~( isnumeric(val) & ( prod(size(val)) == 1 ) )
                  msg = 'Invalid Input.';
             else
                  if strcmp(tag,'Length')
                     val = 2 * floor(val/2) + 1;   % Odd Number !!!
                  end
                  str = sprintf(frm,val);
                  val = eval(str);
                  if ( val < get(hs,'min') ) | ( get(hs,'max') < val ) 
                     s1 = sprintf(frm,get(hs,'min'));
                     s2 = sprintf(frm,get(hs,'max'));
                    msg = sprintf('Value must between %s and %s.',s1,s2);
                  end
             end

             if isempty(msg)
                set( hs , 'value' , val );
                set( h  , 'string' , str , ...
                        'userdata' , str      );
             else
                set( h  , 'string' , get(h,'userdata') );
             end

         end

   %-------------------------------------------------------------
   case 'slider'
   %-------------------------------------------------------------

         he = findobj(fig,'style','edit','tag',tag);

         if ~isempty(he)

             if strcmp(tag,'Length')
                val = 2 * floor(val/2) + 1;   % Odd Number !!!
             end

             frm = get( h , 'userdata' );

             str = sprintf(frm,val);
             val = eval(str);

             set( h  , 'value'  , val );
             set( he , 'string' , str , ...
                     'userdata' , str      );

         end

   %-------------------------------------------------------------
   case 'popupmenu'
   %-------------------------------------------------------------

         str = str{val};
 
   %-------------------------------------------------------------
   otherwise
   %-------------------------------------------------------------

         return

%****************************************************************
end   % slider  <-->  edit  /  popupmenu
%****************************************************************

if ~isempty(msg)
    msgbox(msg,'Error','warn');
    return
end

pointer = get(fig,'pointer');
          set(fig,'pointer','watch');

%****************************************************************
switch tag
%****************************************************************

   %-------------------------------------------------------------
   case 'ColorMap'
   %-------------------------------------------------------------

        nc = size(get(fig,'colormap'),1);

        try
           set(fig,'colormap',cmap_lim(str,nc,[0.5 2.5]));
        catch
           msgbox(sprintf('Can''t evaluate ColorMap "%s".\n%s',str,lasterr), ...
                  'Error','warn');
        end
 
   %-------------------------------------------------------------
   case { 'LightAzim' 'LightElev' }
   %-------------------------------------------------------------

        h = findobj(axe,'type','light','tag','Light');

        v = lightang(h);

        v(1+strcmp(tag,'LightElev')) = val;

        vis = { 'on'  'off' };
        vis = vis{ 1 + ( v(2) < 0 ) };

        lightang( h , v , 'visible' , vis );

   %-------------------------------------------------------------
   case { 'ZScale' 'Azimut' 'Elevation' }
   %-------------------------------------------------------------

        set( axe , 'xtickmode' , 'auto' , ...
                   'ytickmode' , 'auto' , ...
                   'xticklabelmode' , 'auto' , ...
                   'yticklabelmode' , 'auto'       );
                 
        if strcmp(tag,'ZScale')
           set( axe , 'dataaspectratio' , [ 1  1  1/val ] );
        else
           v = get(axe,'view');
           v(1+strcmp(tag,'Elevation')) = val;
           set(axe,'view',v);
        end

        for c = 'xy'
            t = get(axe,[c 'tick']);
            l = get(axe,[c 'ticklabel']);
            l = cellstr(l);
            i0 = ( abs(t) < 1e3*eps );
            if any(i0)
               l(find(i0)) = {upper(c)};
               set( axe , [c 'tick'] , t , [c 'ticklabel'] , l );
            end
        end

   %-------------------------------------------------------------
   otherwise
   %-------------------------------------------------------------

        plat_calc(fig)

%****************************************************************
end  % tag
%****************************************************************

set(fig,'pointer',pointer);

%***************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function plat_calc(fig)

axe = findobj(fig,'type','axes'   ,'tag','Axes');
srf = findobj(fig,'type','surface','tag','Surface');

hx = findobj(fig,'style','slider','tag','RadiusX');
hy = findobj(fig,'style','slider','tag','RadiusY');
hr = findobj(fig,'style','slider','tag','Rotation');
ho = findobj(fig,'style','slider','tag','Offset');

hm = findobj(fig,'style','popupmenu','tag','Mode');

hn = findobj(fig,'style','slider','tag','Length');


if isempty(axe) | isempty(srf) | ...
   isempty(hx)  | isempty(hy)  | ...
   isempty(hr)  | isempty(ho)  | ...
   isempty(hm)  | isempty(hn)

   warning('Can''t find required Handles.');

   return

end


mode =      get(hm,'string');
mode = mode{get(hm,'value')};

n = get(hn,'value');

x = get(axe,'xlim');
y = get(axe,'ylim');

x = linspace(x(1),x(2),n);
y = linspace(y(1),y(2),n)';

r = [ get(hx,'value') get(hy,'value') get(hr,'value')+i*get(ho,'value') ];

z = plateau(x,y,r,[mean(x) mean(y)],mode);

set( srf , 'xdata' , x , ...
           'ydata' , y , ...
           'zdata' , z , ...
           'cdata' , z        )

str = sprintf('WINDOW( %s )',upper(mode) );

set( get(axe,'title') , 'string' , str );


%***************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [fs,ssi,ppi,is_win] = bestfont

% BESTFONT  Returns an optimal FontSize for this computer
%
% [ FontSize , ScreenSize , ScreenPixelsPerInch ] = BESTFONT
%

is_win = strcmp( upper(computer) , 'PCWIN' );

uni = get(0,'units');       
      set(0,'units','pixels')
ssi = get(0,'ScreenSize');  
      set(0,'units',uni);
          
ppi = get(0,'ScreenPixelsPerInch');

is_tall = -1 + ( ssi(4) >=  480 ) + ...
               ( ssi(4) >=  600 ) + ...
             0*( ssi(4) >= 1024 );

fs =  8 + 2 * is_tall - 1 * is_win;

%***************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function wygiwys(fig,ornt)

% WYGIWYS  WhatYouGetIsWhatYouSee
%
% Switch the PaperPosition to the same View like on the Screen
% The PaperOrientation will adapted to the Extension of the Figure.
% The Figure will positioned in the Center of the Paper
%
% WYGIWYS( FigureHandle )
%
% WYGIWYS( FigureHandle , Orientation )
%
%   use the defined PaperOrientation: 'portrait' | 'landscape'
%


if nargin < 1
  fig = get(0,'currentfigure');
end

if isempty(fig)
 return
end

ok = ( isnumeric(fig)  &  ( prod(size(fig)) == 1 ) );
if ok
   ok = ishandle(fig);
   if ok
      ok = strcmp( get(fig,'type') , 'figure' );
   end
end

if ~ok
   error('Input must be a FigureHandle.');
end

mode = { 'portrait'  'landscape' };
if nargin == 2
   if ~( ischar(ornt) & ~isempty(ornt) & ...
         ( prod(size(ornt)) == size(ornt,2) ) )
       error('Orientation must be a String.');
   end
   ornt = lower(ornt(1));
   mode = mode{ 1 + strcmp(ornt,'l') };
   set( fig , 'paperorientation' , mode );
end
   
figuni = get(fig,'units');
papuni = get(fig,'paperunits');

set(fig,     'units' , 'pixels' , ...
        'paperunits' , 'inches'       );

figpos = get(fig,'position');

ppi    = get(0,'screenpixelsperinch');

pappos = zeros(1,4);

pappos([3 4]) = figpos([3 4]) / ppi;

if ~ischar(mode)

    set( fig , 'paperorientation' , mode{1} );

   pap_si = get(fig,'papersize');

   sc     = pap_si ./ pappos([3 4]);
   flip   = ( sc(1) < sc(2) );
   flip   = ( flip & ( ( sc(1) < 1 ) | ( pappos(4) < pappos(3) ) ) );
   mode   = mode{ 1 + flip };

   pap_si = pap_si([1 2]+[1 -1]*flip);

   set( fig , 'paperorientation' , mode );

else

   pap_si = get(fig,'papersize');

end

pappos([1 2]) = (pap_si-pappos([3 4])) / 2;

set(fig,'paperposition',pappos);

set(fig,     'units' , figuni   , ...
        'paperunits' , papuni         );

%***************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function c = cmap_lim(fcn,nc,lim)

% CMAP_LIM  returns bright-limitated Colormap
%
% CMAP_LIM( ColorMap    , ColorNumber , Limit )
% CMAP_LIM( ColorMapFcn , ColorNumber , Limit )
%
% Options:
%
%  imag(ColorNumber) ~= 0  ==>  Invert ColorMap
%     
%  real(ColorNumber) >  0  ==>   Bright --> Dark  
%                                  Blue --> Red  
%
%  real(ColorNumber) <  0  ==>     Dark --> Bright  |
%                                   Red --> Blue
%                

if nargin < 1
   fcn = 'hot';
end

if nargin < 2
   nc = 64;
end

if nargin < 3
   lim = [ 0.5  2.5 ];
end

%------------------------------------------------

if real(nc) == 0
  c = zeros(0,3);
  return
end

invert = ~( imag(nc) == 0 );
flip   =  ( real(nc) <  0 );

nc = abs(real(nc));

%------------------------------------------------

ok = ( ischar(fcn) & ~isempty(fcn) & ...
       ( prod(size(fcn)) == size(fcn,2)   ) ); 

if ok

  try  
    c = feval(fcn,2*max(nc,256));
  catch
    ok = 0;
  end

else

  ok = ( isnumeric(fcn) & ~isempty(fcn) & ...
         ( prod(size(fcn)) == size(fcn,1)*size(fcn,2) ) & ...
         ( size(fcn,2) == 3 ) );

  if ok
     ok = ( all(isfinite(fcn)) & all( fcn(:) >= 0 ) & all( fcn(:) <= 1 ) );
  end

  if ok
     c = fcn;
  end

end

if ~ok
   c = hot(2*max(nc,256));
end

%-----------------------------------------------------

c = 1*invert + (1-2*invert) * c ;

%-----------------------------------------------------
   
cs = sum( c , 2 );
ok = find( ( lim(1) <= cs )  &  ( cs <= lim(2) ) );


if prod(size(ok)) > 1
  c = c(ok,:);
end

ncm = size(c,1);

if ncm == 1

  c = c( ones(nc,1) , : );

else

  c = interp1( (1:ncm)' , c , linspace(1,ncm,nc)' );

end

f = 1 - 2*flip;

if f*sum(c(1,:),2) < f*sum(c(nc,:),2)
% Last == Darkest

   c = c(nc:-1:1,:);

else
% Last == Nearest to Red

  h = rgb2hsv(c);
  h = abs(h(:,1)-0.5);

  if f*h(1,1) > f*h(nc,1)

   c = c(nc:-1:1,:);

  end

end

%***************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function c = rgb(n)

% RGB  Red-Yellow-Green-Cyan-Blue Colormap
%
% ColorMap = RGB( N )
%

if nargin < 1
  fig = get(0,'currentfigure');
  if isempty(fig)
    n = size( get(0,'defaultfigurecolormap') , 1 );
  else
    n = size( get(fig,'colormap') , 1 );
  end
end


n0 = 64;

x = linspace( 0 , pi/2 , n0+1 );
y = sin(x) + 1;
y = cat( 2 , ( 2 - y(end:-1:2) ) , y );
y = cat( 2 , y , y(end-1:-1:1) );
y = y/2;

y = y(:);

ny = size(y,1);  %  n0 == 4*n0+1  !!!

n2 = ( ny - 1 ) / 2;

cc = [ 0  0  1   -1
       0  1  1    0
       0  1  0    1
       1  1  0    2
       1  0  0    3   ];

no = n2/2;  % Offset for Red and Blue
% Overlapp: n2
% Red and Blue n2+no !!!
n0 = 4*n2 + 1 + 2*no;
c0 = zeros( n0 , 3 );

ic = ( 1 : ny );

for ii = 1 : size(cc,1);

  jj = cc(ii,4) * n2 + no + ic;

  ok = find( ( 1 <= jj )  &  ( jj <= n0 ) );
 
  for kk = 1 : 3
    c0(jj(ok),kk) = c0(jj(ok),kk) + y(ok) * cc(ii,kk);
  end

end

c0 = c0 + ( 1 - c0 ) .* ( c0 > 1 );

c0 = c0 .^ 0.5;

c = interp1( ( 1 : n0 )' , c0 , linspace( 1 , n0 , n )' );
