function fig = windemo(varargin)

% WINDEMO   Run Demonstration for WINDOW
%
% Displays a 3D-Surface of the Dependence of the Window,
%  specified by Mode and Length, from the Value of Power
%
% The Window- and ViewParameter can be controled by a GUI.
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
   fig = win_new(fcn);
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

function fig = win_new(fcn)


nc = 128;

ResizeFcn = sprintf('%s(''win_pos'',gcbo);',fcn);

fig = figure( 'colormap' , rgb(nc)    , ...
                 'color' , 'w'    , ...
               'menubar' , 'none' , ...
               'toolbar' , 'none' , ...
           'numbertitle' , 'off' , ...
                  'name' , 'Demo WINDOW' , ...
             'ResizeFcn' , ResizeFcn            );

axe = axes( 'parent' , fig       , ...
               'box' , 'on'      , ...
              'xlim' , [ -1  1 ] , ...
              'zlim' , [  0  1 ] , ...
              'clim' , [  0  1 ] , ...
             'color' , 'none'    , ...
             'layer' , 'top'     , ...
           'tickdir' , 'out'     , ...
          'userdata' , []        , ...    % YTicks
          'nextplot' , 'add'     , ...
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

ptc = patch(   'parent' , axe      , ...
                'xdata' , []       , ...
                'ydata' , []       , ...
                'zdata' , []       , ...
            'linestyle' , '-'      , ...
            'linewidth' ,  1       , ...
            'edgecolor' , 'k'      , ...
            'facecolor' , 'none'   , ...
                  'tag' , 'Patch'      ); 

set( get(axe,'ylabel') , 'string' , 'Power' );

%------------------------------------------------------
% ColorBar

nc = size( get(fig,'colormap') , 1 );

bar = axes( 'parent' , fig       , ...
               'box' , 'on'      , ...
              'xlim' , [ 0  1 ]  , ...
              'ylim' , [ 1 nc ] + 0.5* [ -1  1 ] , ...
              'clim' , [ 0  1 ]  , ...
             'xtick' , []        , ...
     'yaxislocation' , 'right'   , ...
             'color' , 'none'    , ...
             'layer' , 'top'     , ...
          'nextplot' , 'add'     , ...
              'tag'  , 'ColorBar'            );

img = image( 'parent' , bar        , ...
             'xdata' , 0.5         , ...
             'ydata' , [ 1   nc ]' , ...
             'cdata' , ( 1 : nc )' , ...
      'cdatamapping' , 'direct'    , ...
              'tag'  , 'ColorImage'            );

%------------------------------------------------------
% Controls

ini = win_ini;
frm = '%.0f';    % Format for  Edit <---> Slider

hc = ini(:,[1 3]);
  
for ii = ( 1 : size(ini,1) )

    if iscellstr(ini{ii,2})
       style = 'popupmenu';
    else
       style = 'edit';
    end

    CB = sprintf('%s(''win_clb'',gcbo);',fcn);

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

win_pos( fig , 1 );

win_clb( hc.Azimut );
win_clb( hc.Elevation );

win_clb( hc.YScale );
win_clb( hc.ZScale );

drawnow

win_clb( imag(hc.Length) );

%***************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ini = win_ini

% WIN_INI  UIControl-Initialisation

mode = { 'linear'
         'triangle'
         'cosine'
         'gauss'
         'binomial'
         'power'   };

mval = 3;

len  = [ 3  501 ];

yscl = set(0,'DefaultAxesYScale');
zscl = set(0,'DefaultAxesYScale');

yval = strmatch('log',lower(yscl));
if isempty(yval)
   yval = 1;
else
   yval = yval(1);
end

zval = strmatch('lin',lower(zscl));
if isempty(zval)
   zval = 1;
else
   zval = zval(1);
end

zdata = { 'absolute'
          'normalized' };

cdata = { 'absolute'
          'logarithmic'
          'normalized'  };

azm   = [ -180  180 ];
elv   = [  -90   90 ];

%         Tag     String/MinMax   Value
ini   = { 'Mode'       mode     mval
          'Length'     len      101
          'ZData'      zdata      2
          'CData'      cdata      1
          'YScale'     yscl     yval
          'ZScale'     zscl     zval
          'Azimut'     azm      150
          'Elevation'  elv       30  };

%***************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function win_pos(fig,mode)


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

ini = win_ini;
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
   hd = max( floor(hd) , ceil(0,5*hh) );
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

%--------------------------------------------------------------
% Check Line for Changed YTicks

if nargin > 1
   return
end

ok = win_line(fig,1);
if ok
   hz = findobj(fig,'style','popupmenu','tag','ZData');
   if ~isempty(hz)
       win_clb(hz);
   end
end

%***************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function win_clb(h);

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
   case { 'YScale'  'ZScale' }
   %-------------------------------------------------------------

        set( axe , tag , str );

        %-------------------------------------------------------------
        % Lines

        if tag(1) == 'Y'
           h = findobj(fig,'type','patch','tag','Patch');
           if ~isempty(get(h,'userdata'))
               ok = win_line(fig,1);      % Check with AxesYTick
               if ok
                  hz = findobj(fig,'style','popupmenu','tag','ZData');
                  if ~isempty(hz)
                      win_clb(hz);
                  end
               end
           end
        end

   %-------------------------------------------------------------
   case { 'Azimut' 'Elevation' }
   %-------------------------------------------------------------

        v = get(axe,'view');

        v(1+strcmp(tag,'Elevation')) = val;

        set(axe,'view',v);

        %-------------------------------------------------------------
        % Lines

           h = findobj(fig,'type','patch','tag','Patch');
           if ~isempty(get(h,'userdata'))
               ok = win_line(fig,1);      % Check with AxesYTick
               if ok
                  hz = findobj(fig,'style','popupmenu','tag','ZData');
                  if ~isempty(hz)
                      win_clb(hz);
                  end
               end
           end

   %-------------------------------------------------------------
   case { 'ZData' 'CData' }
   %-------------------------------------------------------------

        tags = { 'Surface'  'Patch' };

        for tt = tags( 1 : ( 2 - (tag(1)=='C') ) )

            h = findobj(fig,'type',lower(tt{1}),'tag',tt{1});

            w = get(h,'userdata');

            if isempty(w)
               set(fig,'pointer',pointer);
               return
            end

            switch str(1)
               case 'n'     % normalized
                     if strcmp(tt{1},'Surface')
                        w = w ./ ( max(w,[],2) * ones(1,size(w,2)) );
                     else
                        w = w ./ ( ones(size(w,1),1) * max(w,[],1) );
                     end
                   lim = [ 0  1 ];
               case 'l'     % logarithmic
                     w(find(w==0)) = NaN;
                     w = log(w) / log(10);
                   lim = [ min(w(:))  0 ];
               case 'a'
                   lim = [ 0  max(w(:)) ];
               otherwise
                   warning(sprintf('Invalid Option: %s.',str));
                   set(fig,'pointer',pointer);
                   return
            end

            set(h,tag,w);

            if strcmp(tt{1},'Surface')

               set(axe,[tag(1) 'lim'],lim);

               switch tag(1)

                  case 'Z'

                      set( get(axe,'zlabel') , 'string' , str );

                  case 'C'

                      bar = findobj(fig,'type','axes' ,'tag','ColorBar');
                      img = findobj(fig,'type','image','tag','ColorImage');

                      set( get(bar,'xlabel') , 'string' , str );

                      set( bar , 'ylim' , lim );

                      nc = size( get(fig,'colormap') , 1 );
                      dc = ( lim(2) - lim(1) ) / nc;
                      yc = ( lim(1)+dc/2 : dc : lim(2)-dc/2 )';

                      set(img,'ydata',yc);

               end  % tag(1)

            end  % Surface

        end  % Surface , Patch


   %-------------------------------------------------------------
   otherwise
   %-------------------------------------------------------------

        win_calc(fig)

%****************************************************************
end  % tag
%****************************************************************

set(fig,'pointer',pointer);

%***************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function win_calc(fig)


axe = findobj(fig,'type','axes'   ,'tag','Axes');
srf = findobj(fig,'type','surface','tag','Surface');

hz = findobj(fig,'style','popupmenu','tag','ZData');
hc = findobj(fig,'style','popupmenu','tag','CData');

hm = findobj(fig,'style','popupmenu','tag','Mode');

hn = findobj(fig,'style','slider','tag','Length');


if isempty(axe) | isempty(srf) | ...
   isempty(hz)  | isempty(hc)  | ...
   isempty(hm)  | isempty(hn)

   warning('Can''t find required Handles.');

   return

end


mode =      get(hm,'string');
mode = mode(get(hm,'value'));

new = lower(mode{1}(1:3));
old = get(hm,'userdata');
      set(hm,'userdata',new);

is_pot = strcmp(new,'pot');

n = get(hn,'value');

y = get( srf , 'ydata' );

if isempty(y)

   d = 0.01;
   l = round( abs(log(d)) / log(10) );

   y = ( -l : d : l )';
   y = cat( 1 , 0 , 10.^y );

end

m = size(y,1);

w = zeros(m,n);

if is_pot
   mode = cell(0,0);
end

for ii = 1 : m
    [w(ii,:),x] = window(n,mode{:},y(ii));
end

set( srf , 'xdata' , x , ...
           'ydata' , y , ...
           'zdata' , w , ...
           'cdata' , w , ...
        'userdata' , w        )

set( axe , 'xlim' , [ min(x)  max(x) ] , ...
           'ylim' , [ min(y)  max(y) ]        );
      

if is_pot
   str = sprintf('WINDOW( %.0f , [Power] )',n);
else
   str = sprintf('WINDOW( %.0f , ''%s'' , [Power] )',n,new);
end

set( get(axe,'title') , 'string' , str );

win_line(fig,0);

win_clb(hz);
win_clb(hc);

%---------------------------------------------------------------
% Switch Azimut if Switch to "Power"

if xor(is_pot,strcmp(old,'pot'))

   off = 180 * ( 1 - 2*is_pot );

   ha = findobj(fig,'style','slider','tag','Azimut');

   if ~isempty(ha)
       val = get(ha,'value') + off;
       val = val - 360 * floor(val/360);  % [    0 .. 360 )
       val = val - 360 * ( val >= 180 );  % [ -180 .. 180 )
       set(ha,'value',val);
       win_clb(ha);
   end

end

%***************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = win_line(fig,chk)

ok  = 0;

axe = findobj(fig,'type','axes' ,'tag','Axes');
ptc = findobj(fig,'type','patch','tag','Patch');

hm = findobj(fig,'style','popupmenu','tag','Mode');

hn = findobj(fig,'style','slider','tag','Length');


if isempty(axe) | isempty(ptc) | ...
   isempty(hm)  | isempty(hn)

   warning('Can''t find required Handles.');

   return

end

mode =      get(hm,'string');
mode = mode(get(hm,'value'));

n = get(hn,'value');

y = get( axe , 'ytick' );

ok = ~( chk & isequal( y , get(axe,'userdata') ) );

set(axe,'userdata',y);

if ~ok
    return
end

y = y(:);

m = size(y,1);

w = zeros(m,n);

if mode{1}(1) == 'p'
   mode = cell(0,0);
end

for ii = 1 : m
    [w(ii,:),x] = window(n,mode{:},y(ii));
end

w = permute(w,[2 1]);
x = permute(x,[2 1]) * ones(1,m);
y = ones(n,1) * permute(y,[2 1]);

nn = NaN * ones(1,m);

set( ptc , 'xdata' , cat(1,x,nn) , ...
           'ydata' , cat(1,y,nn) , ...
           'zdata' , cat(1,w,nn) , ...
        'userdata' , cat(1,w,nn)        )

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
