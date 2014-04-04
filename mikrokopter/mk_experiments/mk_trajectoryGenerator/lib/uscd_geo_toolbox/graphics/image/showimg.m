function [Msg,fig,txt,inf,c,cmap] = showimg(varargin)

% SHOWIMG  Opens a Figure to show an Image from File or CData
%
% [ Msg , FigureHandle ] = SHOWIMG( FileName , FigureName , FigureTag );
%
%  The Inputs FigureName and FigureTag are optional.
%
% [ ... , InfoText , InfoStruct , CData , ColorMap ] = SHOWIMG( ... )
%
%  returns an InfoText, an InfoStructure as part of the Output from IMFINFO,
%   the CData-Matrice and the ColorMap.
%
%--------------------------------------------------------------------------
%
% Instead of the Input "FileName" you can give a 
%
%   [ M by N by 3 ]-TrueColor- or a [ M by N ]-GrayScale-Image,
%
% or a StructArray with the Fields:
%
%  'cdata'               :   [ M by N ] | [ M by N by 3 ]
%  'colormap'   | 'cmap' :   [ R  G  B ]     (optional)
%  'colorlimit' | 'clim' :   [ cmin cmax ]   (optional)
%
%--------------------------------------------------------------------------
%
% [ Msg , FigureHandle ] = SHOWIMG( 0 , FigureName , FigureTag );
%
% creates a new SHOWIMG-Figure only, without displaying it.
%
% [ Msg , FigureHandle ] = SHOWIMG( FileName , FigureHandle , FigureTag );
% [ Msg , FigureHandle ] = SHOWIMG( FileName , FigureName , FigureHandle );
%
% Shows the Image in the specified Figure, created by SHOWIMG.
%
%
% [ Msg , FigureHandle ] = SHOWIMG( FileName , ... , [ Flip Flop ] , ... );
%
% Flips the Image vertical   if NonZero "Flip" and 
% Flops the Image horizonatl if NonZero "Flop". 
%
% [ Msg , FigureHandle ] = SHOWIMG( FileName , ... , i*Rotation , ... );
%
% Shows the Image with the specified Rotation, 
%  which has to be a multiple of 90 (deg).
%
%--------------------------------------------------------------------------
%
% A Click with the Right-MouseKey shows a ContextMenu.
%   The actual Location and Color of the Pixel is displayed and
%   you can modify the Image, save the Image or load a new Image.
%
% To use a ZOOM-Functionality the Function AXE_ZOOM is required,
%    which should normaly provided in the "private"-Directory of SHOWIMG.
%
% By Click and Drag with the Left-MouseKey you can Zoom in.
% A Click with the Middle-MouseKey Zoom's out,
% A DoubleClick with the  Left-MouseKey Zoom's to the origin.
%
%--------------------------------------------------------------------------
%
% see also: IMFINFO, IMREAD, READXPM, READPPM, 
%                    IMWRITE, WRTXPM,  WRTPPM, SAVEIMG
%

  
Msg = '';
fig = [];
txt = '';

inf = struct( 'Filename'    , { '' } , ...
              'FileModDate' , { '' } , ...
              'FileSize'    , { [] } , ...
              'Format'      , { '' } , ...
              'Width'       , { [] } , ... 
              'Height'      , { [] } , ...
              'BitDepth'    , { [] } , ...
              'ColorType'   , { '' } , ...
              'Colormap'    , { [] } , ...
              'Index'       , { [] }       );

c    = zeros(0,0);
cmap = zeros(0,3);

off = [ 10  30  10 ];  % Offset [ Left  Top  FloorRound ]
pap = 1;               % PaperOffset [inches]
grey = uint8(214);     % NaN-Color

%--------------------------------------------------

  nl = char(10);

Nin  = nargin;
Nout = nargout;

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
    fcn = fcn(max(find(fsp))+1:end);
end
app = upper(fcn);


%***********************************************************
% Basic InputCheck

if Nin < 1
   file = '';
else
   file = varargin{1};
   if isempty(file)
      file = '';
   end
end

if ~chkstr(file) & chkhndl(file) & ~isequal(file,0)
    f = file;
    ok = chkhndl(f,'figure');
    if ~ok
       [m,f] = recpar(f,'figure');
       ok = ( isempty(m) & ~isempty(f) );
       if ok
          f = f(1);
       end
    end
    if ok
       ok = isappdata(f,app);
    end
    if ok
       ptr = get(f,'pointer');
             set(f,'pointer','watch');
       try
          Msg = clbfcn(f,file,app,fcn,varargin{2:end});
       end
       if ishandle(f), set(f,'pointer',ptr); end
       return
    else
       Msg = sprintf('HandleInput must origin from a %s-Figure',app); 
    end
end

%***********************************************************
% Check for Rotation and Flip/Flop

rot = [];
flp = [];

if Nin >= 2

   ok = ones(1,Nin);

   for ii = 2 : Nin
       if isnumeric(varargin{ii})
          nv = prod(size(varargin{ii}));
          if ( nv == 1 )
             if real(varargin{ii}) == 0
                ok(ii) = 0;
                rot = imag(varargin{ii});
             end
          elseif ( nv == 2 )
             ok(ii) = 0;
             flp = varargin{ii}(:)';
             flp = ~( flp == 0 );
          end
       end
   end

   Nin = sum(ok);
   varargin = varargin(find(ok));

   if ~isempty(rot)
       varargin(ii) = [];
       Nin = Nin - 1;
       if ~( mod(rot,90) == 0 )
           Msg = [ Msg  nl(1:(end*(~isempty(Msg)))) ...
                   'Rotation must be a multiple of 90.' ];
       end
   end

end

%--------------------------------------------------
% Check FigureName

if Nin < 2
   name = '';
else
   name = varargin{2};
   if isempty(name)
      name = '';
   end
end

if ~chkstr(name) 
    if chkhndl(name,'figure',app);
       fig  = name;
       name = '';
    else
        Msg = [ Msg  nl(1:(end*(~isempty(Msg)))) ...
                sprintf('2. Input must be a FigureName or %s-FigureHandle.',app) ];
    end
end

%--------------------------------------------------
% Check FigureTag or Handle

if Nin < 3
   tag = '';
else
   tag = varargin{3};
   if isempty(tag)
      tag = '';
   end
end

if ~chkstr(tag) 
    if chkhndl(tag,'figure',app) & isempty(fig)
       fig = tag;
       tag = '';
    else
       Msg = [ Msg  nl(1:(end*(~isempty(Msg)))) ...
               sprintf('3. Input must be a FigureTag (or %s-FigureHandle).',app) ];
    end
end

%-----------------------------------------------------------

if ~isempty(Msg)
    if Nout == 0
       error(Msg)
    end
    return
end

%***********************************************************

is_file = ~isequal(file,0);    % !!!!!

new_fig = isempty(fig);

%***********************************************************
if is_file
%***********************************************************

    %-------------------------------------------------------
    % Get Data

    if chkstr(file)

       [Msg,inf,c,cmap,txt] = loadimg(file,inf,grey);

    else

       [Msg,inf,c] = getimg(file,inf,grey);

    end

    if ~isempty(Msg) | isempty(c)
        if isequal(Msg,0)
           Msg = '';
        elseif ~isempty(Msg) 
           if isempty(file) & ( Nout == 0 ) & ( Nin == 0 )
              msgbox(Msg,'Error','warn')
              clear Msg
           elseif Nout == 0
              error(Msg)
           end
        end
        if isempty(Msg) & Nout == 0
           clear Msg
        end
        return
    end

    %-------------------------------------------------------

    form = '%.0fx%.0f %s-Image, %.0f Bit, %s';

    txt = sprintf( form , inf.Width , inf.Height , upper(inf.Format) , ... 
                          inf.BitDepth , inf.ColorType ) ;

    auto = isempty(name);

    if auto
       if chkstr(file)
          [p,name,ext] = fileparts(inf.Filename);
          name = [ name  ext ];
       else
          name = 'Image';
       end
    end

    %---------------------------------------------------------------
    % Check Size of Image

    ssi = get(0,'screensize');

    if  inf.Width*inf.Height > 2*ssi(3)*ssi(4)

        wrn = sprintf('%.0f x %.0f',inf.Width,inf.Height);
        wrn = sprintf('Large Image:  %s\n\n   Show it anyway?\n',wrn);

        ok = questdlg(wrn,'Warning','Show','Cancel','Cancel');

        if isequal(ok,'Cancel')
           if Nout == 0
              clear Msg
           end
           return
        end

    end

%***********************************************************
end
%***********************************************************

%--------------------------------------------
% Create New Figure
if new_fig
   fig = newfig(tag,app,fcn,inf,off,pap);
end

%--------------------------------------------
% Set ImageData
if is_file

   if ~isempty(tag)
       set(fig,'tag',tag);
   end

   ud = getappdata(fig,app);

   set_rot(fig,app,ud.rot,ud.flp);

   try
      setimg(fig,app,rot,flp,c,cmap,name,auto,txt,inf);
   catch
      Msg = sprintf('Error set Image.\n%s',lasterr);
   end

   if isempty(Msg)

      ud = getappdata(fig,app);

      set(ud.img,'userdata',c);   % !!! Use for ColorManipulation
    
   else

       if new_fig
          try
             delete(fig);
             fig = [];
          end
       end
       if isempty(file) & ( Nout == 0 ) & ( Nin == 0 )
          msgbox(Msg,'Error','warn')
          clear Msg
       elseif Nout == 0
          error(Msg)
       end

   end

end


if Nout == 0
   clear Msg
end

if is_file | new_fig | ( Nout < 3 )
   return
end

%--------------------------------------------
% Get ImageData

ud = getappdata(fig,app);

txt = ud.txt;
inf = ud.inf;

if Nout < 5
   return
end

c = get(ud.img,'cdata');

cmap = ud.cmap;


%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function fig = newfig(tag,app,fcn,inf,off,pap);


scr_uni = get(0,'units');      set(0,'units','pixels')
scr_si  = get(0,'ScreenSize'); set(0,'units',scr_uni);

figpos        = zeros(1,4);
figpos([3 4]) = 100;

figpos([1 2]) = scr_si([3 4])-off([1 2])-figpos([3 4]);

%-----------------------------------------------------------

 fig = figure( 'paperunits'         , 'inches'   , ...
               'paperorientation'   , 'portrait' , ...
               'units'              , 'pixels'   , ...
               'position'           ,  figpos    , ...
               'color'              , [1 1 1]    , ...
               'numbertitle'        , 'off'      , ...
               'menubar'            , 'none'     , ...
               'toolbar'            , 'none'     , ...
               'name'               , ''         , ...
               'tag'                , tag        , ...
               'colormap'           , ones(1,3)  , ...
               'visible'            , 'off'      , ...
               'handlevisibility'   , 'callback' , ...
               'createfcn'          , ''         , ...
               'resize'             , 'off'             );

%-----------------------------------------------------------

axe = axes( 'parent'   , fig            , ...
            'units'    , 'normalized'   , ...
            'position' , [ 0  0  1  1 ] , ...
            'xlim'     , [ 0  1 ]       , ...
            'ylim'     , [ 0  1 ]       , ...
            'xtick'    , []             , ...
            'ytick'    , []             , ...
            'color'    , 'none'         , ...
            'xgrid'    , 'off'          , ...
            'ygrid'    , 'off'          , ...
            'xdir'     , 'normal'       , ...
            'ydir'     , 'reverse'      , ...
            'view'     , [ 0  90 ]      , ...
     'dataaspectratio' , [ 1  1  1 ]    , ...
            'visible'  , 'off'          , ...
               'tag'   , 'AXES'         , ...
    'handlevisibility' , 'callback'         );


img = image( 'parent' , axe , ...
             'xdata'  , 0.5 , ...
             'ydata'  , 0.5 , ...
             'cdata'  , NaN , ...
       'cdatamapping' , 'direct'  , ...
           'clipping' , 'off' , ...
           'visible'  , 'on'  , ...
              'tag'   , 'IMAGE'             );

lin =  line( 'parent' , axe , ...
             'xdata'  , NaN , ...
             'ydata'  , NaN , ...
             'color'  , 'k' , ...
          'linestyle' , '--'   , ...
          'linewidth' , 0.5    , ...
             'marker' , 'none' , ...
          'erasemode' , 'xor'  , ...
           'clipping' , 'off'  , ...
           'visible'  , 'off'  , ...
              'tag'   , 'LINE'             );

%------------------------------------------------------------
% ContexMenu

form = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 );
fg   = sprintf(form,fig);
tg   = 'Context';

CB = sprintf('%s(gcbo,''%s'');',fcn,'INFO');

hc = uicontextmenu( 'parent' , fig   , ...
                    'tag'    , tg    , ...
                  'callback' , CB    , ...
             'interruptible' , 'off' , ...
                'BusyAction' , 'cancel'      );

hm = mkmenu(hc,'Menu',guimenu,tg,fcn);

set([axe,img],'uicontextmenu',hc);

smt = getappdata(hm.Smooth,'Children');
set( smt.SMT01,'checked','on');

%------------------------------------------------------------

ud = struct( 'auto' , {   1   } , ...
              'txt' , {  ''   } , ...
              'inf' , { inf   } , ...
             'name' , {  ''   } , ...
             'cmap' , { zeros(0,3) } , ...
              'off' , { off   } , ...
              'pap' , { pap   } , ...
              'axe' , { axe   } , ... 
              'img' , { img   } , ...
             'line' , { lin   } , ...
          'context' , { hc    } , ...
             'menu' , { hm    } , ...
              'rot' , {  0    } , ...
              'flp' , { [0 0] } , ...
              'inv' , {  0    } , ...
              'hue' , {  0    } , ...
              'smt' , {  1    } , ...
              'lim' , { [0 1 0 1] } , ...
             'area' , { [0 1 0 1] }        );

setappdata(fig,app,ud);

%------------------------------------------------------------
% AXE_ZOOM

if isempty(which('axe_zoom'))
   set(hm.Help,'visible','off');
   return
end

fcn = { fcn fig 'ZOOM' };
try
  Msg1 = axe_zoom('new',axe,fcn);
  if isempty(Msg1)
     Msg1 = axe_zoom('on',axe);
  end
catch
  Msg1 = lasterr;
end

if isempty(Msg1)
   btd = sprintf('%s(gcbo,''DOWN'');',fcn{1});
   set( [ axe img ] , 'buttondownfcn' , btd );
else
   set(hm.Help,'visible','off');
end


%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [rot,flp] = set_rot(fig,app,rot,flp);

ud = getappdata(fig,app);

hm = ud.menu;

hr = get(hm.Rotation,'children');

hl = char(get(hr,'label'));

jf = ( hl(:,1) == 'F' );     % True for Flip/Flop-Menus

hf = hr(find( jf));
hr = hr(find(~jf));  % Remove Flip/Flop !!!

%----------------------------------------------------------
% Set Rotation

if ~isempty(rot)

    r0t = get(hr,'userdata');

    r0t = cat(1,r0t{:});

    set(hr,'checked','off');

    rot = rot - 360 * floor(rot/360);
    r0t = r0t - 360 * floor(r0t/360);

    ii = ( r0t == rot );

    if any(ii)
       ii = find(ii);
       set(hr(ii),'checked','on');
       ud.rot = rot;
       setappdata(fig,app,ud);
    else
       rot = [];
    end

end

%----------------------------------------------------------
% Get/Set Flip/Flop

nf = prod(size(hf));

jf = get(hf,'userdata');

jf = cat(2,jf{:});

if  isempty(flp)

    flp = zeros(1,nf);

    chk = strcmp( get(hf,'checked') , 'on' );

    flp(jf) = chk;

else

    set(hf,'checked','off');

    for ii = 1 : nf

        if flp(jf(ii))
           set(hf(ii),'checked','on');
        end

    end

    ud.flp = flp;

    setappdata(fig,app,ud);

end

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function Config = guimenu

act = 'PRINT';

%------------------------------------------------------------------

ClipBoard = struct( 'Meta' , { { 'MetaFile'  0  act  {'meta'}   } } , ...
                    'Bmp'  , { { 'BitMap'    0  act  {'bitmap'} } }       );

%------------------------------------------------------------------
% PostScript

 ps = struct( 'RGB'  , { { 'RGB'  0  act  { 'ps' '-loose'        } } } , ...
              'CMYK' , { { 'CMYK' 0  act  { 'ps' '-loose' '-cmyk'} } }        );

eps = struct( 'RGB'  , { { 'RGB'  0  act  {'eps' '-loose'        } } } , ...
              'CMYK' , { { 'CMYK' 0  act  {'eps' '-loose' '-cmyk'} } }        );

PostScript = struct(  'PS' , { {  'PS'  0   ps   [] } } , ...
                     'EPS' , { { 'EPS'  0   eps  [] } }       );


%------------------------------------------------------------------
% Image

 jpg = struct( 'Q100' , { { 'Quality 100%'  0  act  { 'jpg' 100 } } } , ...
               'Q095' , { { 'Quality  95%'  0  act  { 'jpg'  95 } } } , ...
               'Q090' , { { 'Quality  90%'  0  act  { 'jpg'  90 } } } , ...
               'Q085' , { { 'Quality  85%'  0  act  { 'jpg'  85 } } } , ...
               'Q080' , { { 'Quality  80%'  0  act  { 'jpg'  80 } } } , ...
               'Q075' , { { 'Quality  75%'  0  act  { 'jpg'  75 } } }         );

 img = struct(  'JPG' , { { '&JPEG'  0   jpg  []      } } , ...
                'BMP' , { { '&BMP'   0   act  {'bmp'} } } , ...
                'TIF' , { { '&TIFF'  0   act  {'tif'} } } , ...
                'PNG' , { { 'P&NG'   0   act  {'png'} } } , ...
                'PCX' , { { 'P&CX'   0   act  {'pcx'} } } , ...
                'XPM' , { { '&XPM'   0   act  {'xpm'} } } , ...
                'PPM' , { { '&PPM'   0   act  {'ppm'} } }      );
                                 
%------------------------------------------------------------------
% Print-Menu

if strcmp( upper(computer) , 'PCWIN' )

  prn = struct( 'PrintDlg'   , { { '&Drucken'        0   act        {'dlg'} } } , ...
                'ClipBoard'  , { { '&Zwischenablage' 1   ClipBoard          } } , ...
                'PostScript' , { { '&PostScript'     1   PostScript         } } , ...
                'Image'      , { { '&Image'          0   img                } }        );

else

  prn = struct( 'PrintDlg'   , { { '&Drucken'        0   act    {'dlg'} } } , ...
                'PostScript' , { { '&PostScript'     1   PostScript     } } , ...
                'Image'      , { { '&Image'          0   img            } }        );

end

%------------------------------------------------------------------
% Info-Menu

inf = struct( 'HEX' , { { 'HE&X'  0  '' 'HE&X: %s' } } , ...
              'RGB' , { { '&RGB'  0  '' '&RGB: %s' } } , ...
              'HSV' , { { '&HSV'  0  '' '&HSV: %s' } }       );

%------------------------------------------------------------------
% Help-Menu

lft = 'Click and Drag with &Left-MouseKey to Zoom in.';
dbl = '&DoubleClick with Left-MouseKey to Zoom to origin.';
mdl = 'Click with &Middle-MouseKey to Zoom out.';

hlp = struct( 'lft' , { { lft  0  ''  [] } } , ...
              'dbl' , { { dbl  0  ''  [] } } , ...
              'mdl' , { { mdl  0  ''  [] } }        );


%------------------------------------------------------------------
% RotationMenu

act = 'ROTATE';
flp = 'FLP';

g = char(176);

rot = struct( 'R270' , { { ['-090' g]    0  act  270 } } , ...
              'R000' , { { ['+000' g]    0  act   00 } } , ...
              'R090' , { { ['+090' g]    0  act   90 } } , ...
              'R180' , { { ['+180' g]    0  act  180 } } , ...
              'Flip' , { { 'Fl&ip   |'   1  flp    1 } } , ...
              'Flop' , { { 'Fl&op --'    0  flp    2 } }       );
         
%------------------------------------------------------------------
% ColorMenu

act = 'SET';

rs = ( 0 : 30 : 360-30 );

n  = size(rs,2);
r0 = cell(2,n);
for ii = 1 : n
    lab = sprintf('+%3.3d%s',rs(ii),char(176));
    r0{1,ii} = sprintf('S%3.3d',rs(ii));
    r0{2,ii} = { { lab  0  act   rs(ii) } };
end

act = { 'Right' 'Left' };
 
rr = [ 1  2  5  10  20  30  45  60  90 ];

n  = size(rr,2);
r1 = cell(2,n);
r2 = cell(2,n);
for ii = 1 : n
    lab = sprintf('+%2.2d%s',rr(ii),char(176));
    r1{1,ii} = sprintf('R%3.3d',rr(ii));
    r1{2,ii} = { { lab  0  act{1}   rr(ii) } };
    lab = sprintf('-%2.2d%s',rr(ii),char(176));
    r2{1,ii} = sprintf('R%3.3d',rr(ii));
    r2{2,ii} = { { lab  0  act{2}  -rr(ii) } };
end

col = struct( 'Set'  , { { '&Set'            0  struct(r0{:}) [] } } , ...
              'Rotl' , { { 'Rotate &right'   0  struct(r1{:}) [] } } , ...
              'Rotr' , { { 'Rotate &left'    0  struct(r2{:}) [] } }       );


act = 'Light';

ls = [ 50  20  10  0  -10 -20 -50 ];
n  = size(ls,2);
l0 = cell(2,n);
for ii = 1 : n
    lab = sprintf('%s%2.2d%%',char(44-sign(ls(ii))),abs(ls(ii)));
    lab = strrep(lab,char(44),'0');
    l0{1,ii} = sprintf('L%s%3.3d',char('C'-sign(ls(ii))),abs(ls(ii)));
    l0{2,ii} = { { lab  0  act   ls(ii) } };
end

act = { 'Brite'  'Dim' };
 
ll = [ 1  2  5  10  20  50 ];

n  = size(ll,2);
l1 = cell(2,n);
l2 = cell(2,n);
for ii = 1 : n
    lab = sprintf('+%2.2d%%',ll(ii));
    l1{1,ii} = sprintf('B%3.3d',ll(ii));
    l1{2,ii} = { { lab  0  act{1}   ll(ii) } };
    lab = sprintf('-%2.2d%%',ll(ii));
    l2{1,ii} = sprintf('D%3.3d',ll(ii));
    l2{2,ii} = { { lab  0  act{2}  -ll(ii) } };
end

lgt = struct( 'Set'   , { { '&Set'     0  struct(l0{:}) [] } } , ...
              'Brite' , { { '&Brite'   0  struct(l1{:}) [] } } , ...
              'Dim'   , { { '&Dim'     0  struct(l2{:}) [] } }       );

rud = ['&Rotation %.1f' g ' '];
lud = '&Light   %s%2.2d%% ';

col = struct( 'Invert' , { { '&Invert'        0   'INV'   []  } } , ...
              'Rotate' , { { '&Rotate ... '   1    col    rud } } , ...
              'Light'  , { { '&Light ...'     1    lgt    lud } }    );

%------------------------------------------------------------------
% SmoothMenu

act = 'SMOOTH';

smt = struct( 'SMT01' , { { ' &none'   0  act   01 } } , ...
              'SMT03' , { { '0&3 px'  0  act   03 } } , ...
              'SMT05' , { { '0&5 px'  0  act   05 } } , ...
              'SMT07' , { { '0&7 px'  0  act   07 } } , ...
              'SMT09' , { { '0&9 px'  0  act   09 } } , ...
              'SMT11' , { { '&11 px'  0  act   11 } }       );

%*****************************************************************

Config = struct( 'Info' , { { '&Info'          0  inf   '&Pixel: %3.3d , %3.3d ' } } , ...
             'Rotation' , { { '&Rotation ... ' 1  rot      []  } } , ...
                'Color' , { { '&Color ... '    1  col      []  } } , ...
               'Smooth' , { { '&Smooth ... '   1  smt      []  } } , ...
                 'Save' , { { '&Save as ... '  1  prn      []  } } , ...
               'Export' , { { '&Export ... '   0  img      []  } } , ...
                 'New'  , { { '&New ... '      1  'New'    []  } } , ...
                 'Help' , { { '&Help ... '     1  hlp      []  } } , ...
                'Close' , { { '&Close'         1  'Close'  []  } } );

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function setimg(fig,app,rot,flp,c,cmap,name,auto,txt,inf)

Nin = nargin;

ud = getappdata(fig,app);

zpp = 'AXE_ZOOM';

is_zoom = isappdata(ud.axe,zpp);

new_img = ( Nin >= 5 );

if new_img

   ud.txt  = txt;
   ud.inf  = inf; 
   ud.cmap = cmap;
   ud.name = name;

   set_name =  ( auto | ud.auto ); 

   ud.auto = auto;

else

   c    = get(ud.img,'cdata');
   cmap = ud.cmap;
   inf  = ud.inf;
   auto = ud.auto;
   name = ud.name;

   set_name = auto;

end

%********************************************************
% Flip first if NewImage  or FlipInput

is_flp = ( new_img | ~isempty(flp) );

fl0 = ud.flp;

if is_flp

   if isempty(flp);
      flp = fl0;
      swp = flp;
   else
      ud.flp = flp;
      swp = xor(flp,fl0);
   end

   if swp(1), c = c(end:-1:1,:,:); end
   if swp(2), c = c(:,end:-1:1,:); end

   set_rot(fig,app,[],flp);

end

%********************************************************
% Check for Rotation

r0t = ud.rot;

if ~isempty(rot)
    rot = set_rot(fig,app,rot,[]);
end

is_rot = ~isempty(rot);
if is_rot
   is_rot = ~( rot == r0t );
end

if is_rot
   [nr,c] = rotimg( c , r0t-rot );
   is_rot = ~( nr == 0 );
else
   nr = 0;
end

%********************************************************
% Check for Flip

no_flp = isempty(flp);

if no_flp
   flp = fl0;
end

if mod(nr,2) == 1
   flp = flp([2 1]);
   set_rot(fig,app,[],flp);
   ud.flp = flp;
end

if ~( new_img | is_rot | is_flp )
    return
end

if is_rot
   ud.rot = rot;
end

%********************************************************

scr_uni = get(0,'units');      set(0,'units','pixels')
scr_si  = get(0,'ScreenSize'); set(0,'units',scr_uni);

%--------------------------------------------
% Maximum ImageSize

off    = ud.off;  % [ Left Top  FloorRound ];

scr_si = scr_si([3 4]) - 2 * off([1 2]);    
scr_si = off(3) * floor( scr_si / off(3) );
scr_si = max(scr_si,off(3));

si = [ size(c,2) size(c,1) ];

fsc = min( 1 , min( scr_si ./ si ) );   % Scale Image

figpos        = zeros(1,4);
figpos([3 4]) = fsc * si;

set(fig,'units','pixels');
fp = get(fig,'position');

%--------------------------------------------
% Fit to upper Left

figpos(1) = fp(1);
figpos(2) = fp(2) + fp(4) - figpos(4);

%--------------------------------------------
% [ Bottom Left ] positive

figpos([1 2]) = figpos([1 2]) - min(0,figpos([1 2]));

%--------------------------------------------
% Fit to ScreenSize

offs = figpos([1 2]) + figpos([3 4]) - ( scr_si([1 2]) + off([1 2]) );

figpos([1 2]) = figpos([1 2]) - max(0,offs);

figpos([1 2]) = max( figpos([1 2]) , off([1 2]) );

figpos(1) = min( figpos(1) , scr_si(1)-200 );

%--------------------------------------------
% PaperPosition

mode = { 'portrait'  'landscape' };


set( fig , 'paperunits'       , 'inches' , ...
           'paperorientation' , 'portrait' );

pap_si = get(fig,'papersize');          % [ Width  Height ]  
ppi    = get(0,'screenpixelsperinch');

img_si = pap_si - 2*ud.pap;             % Max. [ Width Height ]
fig_si = figpos([3 4]) / ppi;           % Act. [ Width Height ]

psc = min( img_si./fig_si , 1 );

flip = ( psc(1) < psc(2) );

psc = min(psc);

if flip

   pscf = min( min( img_si([2 1])./fig_si , 1 ) );

   flip = ( pscf > psc );

   if flip
      img_si = img_si([2 1]);
      pap_si = pap_si([2 1]);
      psc    = pscf;
   end

end

mode = mode{ 1 + flip };

pappos([3 4]) = fig_si .* psc;
pappos([1 2]) = ( pap_si - pappos([3 4]) ) / 2;


%--------------------------------------------
if set_name & auto

   wh = sprintf( '%.0fx%.0f' , inf.Width , inf.Height);

   if ~( fsc == 1 )
       wh = sprintf( '%s (%.0fx%.0f / %.0f%%)', wh , ...
                     figpos([3 4]) , round(fsc*100) );
   end

   name = sprintf('%s   %s   %s' , name , wh , upper(inf.Format) );

end

%--------------------------------------------

set( fig , 'units'            , 'pixels'   , ...
           'position'         , figpos     , ...
           'paperorientation' , mode       , ...
           'paperunits'       , 'inches'   , ...
           'paperposition'    , pappos     , ...
           'colormap'         , cmap       , ...
           'visible'          , 'on'             );

if set_name
   set(fig,'name',name);
end

%-----------------------------------------------------------

offs = 1/2 * [ -1  1 ];   % Will added to AXE_ZOOM-History below !!!

xl = [ 1  si(1) ];
yl = [ 1  si(2) ];

set( ud.axe , 'xlim'     , xl+offs , ...
              'ylim'     , yl+offs , ...
   'dataaspectratiomode' , 'auto'   );

set( ud.img , 'xdata'  , xl  , ...
              'ydata'  , yl' , ...
              'cdata'  , c         );

%-----------------------------------------------------------

ud.lim = [ xl  yl ] + offs([1 2 1 2]);

if new_img
   ud.area = ud.lim;
end

setappdata(fig,app,ud);

%-----------------------------------------------------------
% Check for indexed ==> XPM

is_ind = ~( isempty(inf.Colormap) | isempty(inf.Index) );

sets = { 'off'  'on' };
sets = sets{ 1 + is_ind };

hm = getappdata(ud.menu.Save,'Children');
hm = getappdata(hm.Image,'Children');

set( [ hm.PCX  hm.XPM ] , 'visible' , sets );

%-----------------------------------------------------------
% Check FigurePositon

drawnow

set(fig,'position', 2*figpos - get(fig,'position') );
 
%------------------------------------------------------------
% AXE_ZOOM ==>  check History

if is_zoom 
   aud = getappdata(ud.axe,zpp);
   %--------------------------------------------
   if new_img          % New History
   %--------------------------------------------
      aud.History = [ xl  yl ] + offs([1 2 1 2]);
   %--------------------------------------------
   elseif is_rot       % Rotate History
   %--------------------------------------------
      if mod(nr,2) == 1
         aud.History = aud.History(:,[3 4 1 2]);
      end
      if nr >= 2
         aud.History(:,[1 2]) = xl(2) - aud.History(:,[2 1]) + xl(1);
      end
      if nr <= 2
         aud.History(:,[3 4]) = yl(2) - aud.History(:,[4 3]) + yl(1);
      end
      set( ud.axe , 'xlim' , aud.History(end,[1 2]) , ...
                    'ylim' , aud.History(end,[3 4])       );
   %--------------------------------------------
   end
   %--------------------------------------------
   setappdata(ud.axe,zpp,aud);
end

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [nr,c] = rotimg(c,rot);

% ROTIMG  Rotate image

nr = round(rot/90);

nr = nr - 4 * floor(nr/4);

if ( nr == 0 );
   return
end

for ii = 1 : nr
    c = permute(c,[2 1 3]);
    c = c(end:-1:1,:,:);
end

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function hueimg(fig,app,hue,c)

% HUEIMG  Switch HUE-Value of Image

ud = getappdata(fig,app);

if isequal(hue,ud.hue)
   return
end

if isempty(hue)
   hue = ud.hue;
end

%-------------------------------------------------------

h = real(hue);
l = imag(hue);

l = sign(l) * min(abs(l),100);

hok = ~( mod(h,360) == 0 );
lok = ~any( abs(l) == [ 0  100 ] );

%-------------------------------------------------------
% Get / Smooth ImageData

if isempty(c)

   c = get(ud.img,'userdata');

   ssi = get(0,'screensize');

   csi = [ size(c,2) size(c,1) ];

   if  ( prod(csi) > 0.5*ssi(3)*ssi(4) ) & ( hok | lok )

        wrn = sprintf('%.0f x %.0f',csi);
        wrn = sprintf('Large Image:  %s\n\n   Transform anyway?\n',wrn);

        ok = questdlg(wrn,'Warning','Transform','Cancel','Cancel');

        if isequal(ok,'Cancel')
           return
        end

    end
   
    [ok,c] = smtimg(fig,app,[],c);

end

%---------------------------------------------------------------
% Rotate and Flip/Flop

[nr,c] = rotimg( c , -ud.rot );

if ud.flp(1), c = c(end:-1:1,:,:); end
if ud.flp(2), c = c(:,end:-1:1,:); end

%---------------------------------------------------------------

if hok | lok

   if ud.inv
      c = 1 - double(c)/255;
   end

   c = rgb2hsv(c);

   if hok

      c(:,:,1) = c(:,:,1) + h/360;

      c(:,:,1) = c(:,:,1) - floor(c(:,:,1));

   end

   if lok

      ci = 2 + ( l < 0 );  % [ Sat | Value ]
      cv = 1 - abs(l/100);

      c(:,:,ci) = c(:,:,ci) * cv;
      if ci == 2
         c(:,:,3)  = c(:,:,3) + ( 1 - c(:,:,3) ) * ( 1 - cv );
      end

   end

   c = uint8( round( 255 * hsv2rgb(c) ) );

elseif abs(l) == 100

   c = uint8( 255*(1+sign(l))/2 + zeros(size(c)) );

elseif ud.inv

   c = uint8( 255 - double(c) );

end

set( ud.img , 'cdata' , c );

ud.hue = hue;

setappdata(fig,app,ud);

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,c] = smtimg(fig,app,smt,c)

% SMTIMG  Smooth Image

ok = 1;

ud = getappdata(fig,app);

if isequal(smt,ud.smt)
   return
end

if isempty(smt)
   smt = ud.smt;
end

smt = 2 * floor(smt/2) + 1;

%---------------------------------------------------------------
% Get ImageData, Check in Size

if isempty(c)

   c = get(ud.img,'userdata');

   ssi = get(0,'screensize');

   csi = [ size(c,2) size(c,1) ];

   if  ( prod(csi) > exp(-1)*ssi(3)*ssi(4) ) & ( smt > 1 )

        wrn = sprintf('%.0f x %.0f',csi);
        wrn = sprintf('Large Image:  %s\n\n   Smooth anyway?\n',wrn);

        ok = questdlg(wrn,'Warning','Smooth','Cancel','Cancel');

        ok = ~isequal(ok,'Cancel');
        if ~ok
            return
        end

    end

end

%---------------------------------------------------------------
% Smooth

if smt > 1
   for ii = 1 : size(c,3)
       c(:,:,ii) = uint8(round(meanind2(double(c(:,:,ii)),smt,'cos')));
   end
end

%---------------------------------------------------------------
% Set CData using HUEIMG

if nargout < 2

   ud.smt = smt;

   setappdata(fig,app,ud);

   hueimg(fig,app,[],c);

end

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,inf,c] = getimg(c,inf,grey)

Msg = '';

cmap = zeros(0,3);
clim = zeros(0,2);

if isstruct(c)

   if ~( prod(size(c)) == 1 )
       Msg = 'Length of StructArray must be 1.';
       c = zeros(0,0);
       return
   end

   f = fieldnames(c);
   fl = lower(f);

   ok = strcmp( fl , 'cdata' );
   if ~any(ok)
       Msg = 'StructArray must contain Field "cdata".';
       c = zeros(0,0);
       return
   end

   %-----------------------------------------------
   % Check for ColorMap and ColorLimit first

   cm = ( strcmp( fl , 'colormap' ) | ...
          strcmp( fl , 'cmap' ) );
   if any(cm)

      cm   = find(cm);
      cmap = getfield(c,f{cm(1)});

      cm = ( strcmp( fl , 'colorlimit' ) | ...
             strcmp( fl , 'clim' ) );

      if any(cm)
         cm   = find(cm);
         clim = getfield(c,f{cm(1)});
      end

   end

   %-----------------------------------------------
   % Get CData

   ok = find(ok);
    
    c = getfield(c,f{ok});

end

%**************************************************

inf.Width  = size(c,2);
inf.Height = size(c,1);

if isempty(c)
   return
end

%**************************************************
% Check CData

if ~isnumeric(c)
    Msg = 'CData must be numeric.';
    return
end

s3 = size(c,3);

if ~( any( ndims(c) == [ 2  3 ] ) & ...
      any(       s3 == [ 1  3 ] )       );
    Msg = 'CData must be a 2-D indexed or 3-D TrueColor-Matrix (RGB).';
    return
end

%**************************************************
% Check ColorMap

if  isempty(cmap) | ( s3 == 3 )
    cmap = zeros(0,3);
elseif ~isempty(cmap)
    if ~( ( size(cmap,2) == 3 ) & ( ndims(cmap) == 2 ) & ...
          all( abs(cmap(:)-0.5) <= 0.5 ) )
        Msg = 'ColorMap must have 3 Columns with Values between 0 and 1.';
        return
    end
end

%**************************************************
% Check ColorLimit

if  isempty(clim) | ( s3 == 3 )
    clim = zeros(0,3);
elseif ~isempty(clim)
    ok = ( ( prod(size(clim)) == 2 ) & all(isfinite(clim)) );
    if ok
       ok = ( clim(1) < clim(2) );
    end
    if ~ok
        Msg = 'ColorLimit must have 2 finite increasing Values.';
        return
    end
end

%**************************************************

[c,inf.Index,cmap,inf.ColorType] = chkimg(c,cmap,grey,clim);

inf.BitDepth = get(0,'screendepth');


%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [c,ci,cmap,ct] = chkimg(c,cmap,grey,clim);

if nargin < 4
   clim = zeros(0,2);
end

ci = [];

%**************************************************
% Check Class of CData

cl = class(c);

if ~any(strcmp(cl,{'double' 'uint8'}))
    cl = 'double';
    c  = feval(cl,c);
end

s3 = size(c,3);

%*************************************************
% Check 2D-CData: Direct | GrayScale
if ( s3 == 1 )

   nc = size(cmap,1);

   if strcmp(cl,'uint8')
      is_nan = [];
   else
      is_nan = find(isnan(c));
      if ~isempty(is_nan)
          c(is_nan) = 1;
      end
   end

   if strcmp(cl,'double') & ~isempty(clim)
      dc = ( clim(2) - clim(1) ) / nc;
      c  = floor( ( c - clim(1) ) / dc ) + 1;
      c  = min(max(c,1),nc);
   end
   
   if isempty(cmap)

      if strcmp(cl,'double')
         c = uint8(round(c*255));
      end

      if ~isempty(is_nan)
          c(is_nan) = grey;
      end

      c  = c(:,:,ones(1,3));

      ct = 'grayscale';

   else

      if strcmp(cl,'uint8');
         c = double(c) + 1;
      end

      c = min(max(c,1),nc);

      ci = c;

      if ~isempty(is_nan)
          ci(is_nan) = NaN;
          cmap = cat(1,cmap,double(grey(ones(1,3)))/255);
          c(is_nan) = nc+1;
          nc        = nc+1;
      end

      cm = uint8(round(cmap*255));

      ct = 'indexed';

      c = cm(cat(3,c,c+nc,c+2*nc));

   end

%*************************************************
% Check 3D-CData
elseif ( s3 == 3 )

   if  strcmp(cl,'double')
       c(find(isnan(c))) = 0;
       c = uint8(round(min(max(c,0),1)*255));
   end

   ct = 'truecolor';

end

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,inf,c,cmap,txt] = loadimg(file,inf,grey)

Msg = '';

c    = zeros(0,0);
cmap = zeros(0,3);
txt  = '';

%--------------------------------------------------

ok = isempty(file);
if ~ok
    ok = any( file == '*' );
end

if ok

   if isempty(file)
      file = '*.*';
   end

   [f,p] = uigetfile(file,'Select an ImageFile');
   if isequal(f,0)
      return
   end

   file = fullfile(p,f);

elseif ~( exist(file,'file') == 2 );

    Msg = sprintf('File not found: %s',file);
    return

end

%--------------------------------------------------

 d  = dir(file);
 ok = ( prod(size(d)) == 1 );
 if ~ok
     ff = which(file);
      d = dir(ff);
     ok = ( prod(size(d)) == 1 );
     if ok
        file = ff;
     end
 end

 if ~ok
     Msg = sprintf('Can''t list File, using DIR(''%s'').',file);
     return
 end

%***********************************************************

inf.Filename    = d.name;
inf.FileModDate = d.date;
inf.FileSize    = d.bytes;


%***********************************************************

[p,n,ext] = fileparts(file);

try

  %------------------------------------------
  % Try IMREAD

             s = imfinfo(file);

  [ c , cmap ] = imread(file);

  [Msg1,inf] = structcmp( inf , s );

  ok = 1;

  [c,inf.Index,cmap] = chkimg(c,cmap,grey);

catch

  switch lower(ext)

   %-------------------------------------------
   case '.xpm'

    %------------------------------------------
    % Try READXPM for XPM
 
    [ Msg , c , cmap , ccm ] = readxpm(file);

    ok = isempty(Msg);

    if ok

       inf.Format    = 'xpm';
       inf.Width     = size(c,2);
       inf.Height    = size(c,1);
       inf.BitDepth  = 8 * size(ccm,2);
       inf.ColorType = 'indexed'; 
       inf.Colormap  = cmap;

       if size(c,3) > 1
          c = c(:,:,1);
       end

       [c,inf.Index,cmap] = chkimg(c,cmap,grey);

    end

   %-------------------------------------------
   case '.ppm'

    %------------------------------------------
    % Try READPPM for PPM
 
    [ Msg , c ] = readppm(file);

    ok = isempty(Msg);

    if ok

       inf.Format    = 'ppm';
       inf.Width     = size(c,2);
       inf.Height    = size(c,1);
       inf.BitDepth  = 24;
       inf.ColorType = 'truecolor'; 
       inf.Colormap  = cmap;

       if size(c,4) > 1
          c = c(:,:,:,1);
       end

    end

   %-------------------------------------------
   case { '.ps' '.psc' '.eps' '.epsc' }

    %------------------------------------------
    % Try GV/GHOSTVIEW for PS

     ini = '%!PS-Adobe';

     i0  = 3;   % Start in <ini> for Format

     ini = double(ini(:));

     ni  = size(ini,1);
 
     fid = fopen(file,'r');

     ok = ~( fid == -1 );

     if ok

        % Read the first <ni> Characters
        bb = fread( fid , ni , 'char' );

        bb = double(bb(:));

        ok = isequal( bb , ini );

        if ok
   
           bb = fread( fid , 20 , 'char' );
 
           bb(find(bb==13)) = [];
           
           bb = cat( 1 , bb(:) , 10 );

           jj = find( bb == 10 ) - 1;

           ini = cat( 1 , ini(i0:ni) , bb( 1 : jj(1) ) );

           inf.Format = char(ini');

        end

        fclose(fid);

     end

     if ok

       %-------------------------------------------------------------
       if isunix

          for cc = { 'gv'  'ghostview' }
              [s,w] = unix( cat( 2 , 'which ' , cc{1} ) );
              if s == 0
                 [s,w] = unix( cat( 2 , cc{1} , ' ' , file , ' &' ) );
                 if ( s == 0 )
                    break
                 end
              end
          end
        
       %-------------------------------------------------------------
       elseif strcmp( upper(computer) , 'PCWIN' )

           [s,w] = dos( cat( 2 , 'gsview32 ' , file ) );

       end

       txt = cat( 2 , 'PostScript: ' , inf.Format );

       ok = ( s == 0 );

       if ok
          Msg = 0;
       end

     end
     % PS ok

   %-------------------------------------------
   otherwise

     ok = 0;

  end
  % switch

end
% try

%------------------------------------------------

if ~ok

  Msg = 'Unable to determine Image-File-Format.';

end

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function Msg = clbfcn(fig,h,app,fcn,varargin)

Msg = '';
nl = char(10);

Msg0 = sprintf('%s: ',app);

Nin = nargin - 3;
if Nin < 1
   Msg = sprintf('%sNot enough InputArguments.',Msg0);
end

action = '';
if chkstr(varargin{1},1);
  action = varargin{1};
end

if ~isempty(action)
    Msg0 = sprintf('%s(%s): ',app,upper(action));
end

ud = getappdata(fig,app);

%***********************************************       
switch upper(action)

%***********************************************       
case 'DOWN'

     sel = get(fig,'selectiontype');

     if ~strcmp(sel,'alt')
         axe_zoom('Down',ud.axe,1);
     end

case 'ZOOM'

      if Nin < 3
         Msg = sprintf('%sNot enough InputArguments.',Msg0);
         return
      end

      axe = varargin{2};
      lim = varargin{3};

      if ~isequal(ud.axe,axe)
          Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
                  'Invalid AxesHandle.' ];
      end

      if ~isequal(size(lim),[1 4])
          Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
                  'Limit must be a [ 1 x 4 ] - Vector.' ];
      end

      if ~isempty(Msg)
          Msg = sprintf('%sInvalidInputs.',Msg0,Msg);
          return
      end

      op  = [ -1  1  -1  1 ];
      lim = op .* ceil( op .* (lim+0.5) ) - 0.5;
      lim = op .* min( lim.*op , ud.lim.*op );

      ud.area = lim;

      d0 = ud.lim([2 4]) - ud.lim([1 3]);
      d1 =    lim([2 4]) -    lim([1 3]);

      d1 = max( d1 , d0./d0([2 1]).*d1([2 1]) );

      lim = ( lim([1 1 3 3]) + lim([2 2 4 4]) + op.*d1([1 1 2 2]) ) / 2;

      dl = op.*ud.lim - op.*lim;
      dl = op .* min( dl , 0 );

      dl = dl( [1 1 3 3] + ( dl([1 1 3 3]) == 0 ) );

      lim = lim + dl;

      lim = op .* min( lim.*op , ud.lim.*op );

       set( axe , 'xlim' , lim([1 2]) , ...
                  'ylim' , lim([3 4]) , ...
     'dataaspectratiomode' , 'auto'   );

      udz = getappdata(axe,'AXE_ZOOM');
      udz.History(end,:) = lim;
      setappdata(axe,'AXE_ZOOM',udz);

      if ~( isequal(ud.lim,ud.area) | ( size(udz.History,1) == 1 ) )
          set( ud.line , 'xdata' , ud.area([1 2 2 1 1]) , ...
                         'ydata' , ud.area([3 3 4 4 3]) , ...
                       'visible' , 'on' );
      else
          set( ud.line , 'xdata' , NaN , ...
                         'ydata' , NaN , ...
                       'visible' , 'off' );
      end

      setappdata(fig,app,ud);

%***********************************************       
case 'INFO'

      cp = get(ud.axe,'currentpoint');
      cp = floor( cp(1,[1 2]) + 0.5 );

      cc = get(ud.img,'cdata');

      sc = [ size(cc,2) size(cc,1) ];

      cp = min( cp , sc );

      nr = round( ud.rot / 90 );
      nr = nr - 4 * floor(nr/4);

      pp = cp;
      for ii = 1 : nr
          pp(1) = sc(1) - pp(1) + 1;
          sc = sc([2 1]);
          pp = pp([2 1]);
      end

      set( ud.menu.Info , 'label' , sprintf(get(ud.menu.Info,'userdata'),pp) );

      cc = get(ud.img,'cdata');
      cc = double(permute(cc(cp(2),cp(1),:),[1 3 2]));

      ch = getappdata(ud.menu.Info,'Children');

      hex = cat(2,char(32*ones(3,2)),dec2hex(cc));
      hex = permute(hex,[2 1]);
      hex = hex(:)';

      rgb = sprintf('%3.3d ',cc);
      hsv = sprintf('%3.3d ',round(255*rgb2hsv(cc/255)));

      set( ch.HEX , 'label' , sprintf(get(ch.HEX,'userdata'),hex) );
      set( ch.RGB , 'label' , sprintf(get(ch.RGB,'userdata'),rgb) );
      set( ch.HSV , 'label' , sprintf(get(ch.HSV,'userdata'),hsv) );


      ch = getappdata(ud.menu.Color,'Children');

      lr = sprintf(get(ch.Rotate,'userdata'),real(ud.hue));

      set( ch.Rotate , 'label' , lr );

      vv = imag(ud.hue);

      ll = sprintf(get(ch.Light,'userdata'),char(44-sign(vv)),abs(vv));
      ll = strrep(ll,char(44),'0');

      set( ch.Light , 'label' , ll );

      ch = getappdata(ch.Light,'Children');

      sets = { 'on'  'off' };

      for ff = { 'Brite'  'Dim' }

          hh = getappdata(getfield(ch,ff{1}),'Children');

          for f = fieldnames(hh)'

              h = getfield(hh,f{1});
              v = get(h,'userdata');

              set( h , 'enable' , sets{ 1 + ( abs(vv+v) > 100 ) } );
 
          end
  
      end

%***********************************************       
case 'PRINT'

     [m,par] = recpar(h,'uimenu');

     if ~isempty(m)
         return
     end

     prn = ~any(strmatch('&Export',get(par,'label')));

     if prn
        op = [ -1  1 ];
        op = op([1 2 1 2]);
        area = op .* floor( ud.area .* op );
        area = { (area(3):area(4))  (area(1):area(2))  ':' };
     end

     opt = get(h,'userdata');
     ok = ( ~isempty(opt) & iscell(opt) );
     if ok
        ok = chkstr(opt{1},1);
     end
     if ~ok
         return
     end

     typ = opt{1};
     opt = opt(2:end);

     is_win = strcmp( upper(computer) , 'PCWIN' );

     msg = '';

     str = 'Error save Image';

     %----------------------------------------------
     % PRINTDLG

     if strcmp(typ,'dlg')
        vis = get(ud.line,'visible');
              set(ud.line,'visible','off')
        try
           f = printdlg(fig);
           if ~is_win
               waitfor(f);
           end
        catch
           msg = sprintf('Error call PRINTDLG.\n%s',lasterr);
        end
        if ~isempty(msg)
            msgbox(msg,str,'warn');
        end
        set(ud.line,'visible',vis)
        return
     end

     %----------------------------------------------

     form = sprintf( '-f%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 );
     fg   = sprintf(form,fig);

     dev = sprintf('-d%s',typ);
     if any(strcmp(typ,{'ps' 'eps'}))
        dev = sprintf('%sc',dev);
     end

     %----------------------------------------------
     % Windows-ClipBoard

     if any(strcmp(typ,{'meta'  'bitmap'}))
        if ~is_win
            msg = sprintf('%s only supported under WindowsSystem.',upper(dev));
        else
           vis = get(ud.line,'visible');
                 set(ud.line,'visible','off')
           try
              print( fg , dev );
           catch
              msg = sprintf('Error call PRINT(''%s'').\n%s',dev,lasterr);
           end
           set(ud.line,'visible',vis)
        end
        if ~isempty(msg)
            msgbox(msg,str,'warn');
        end
        return
     end

     %----------------------------------------------
     % Get Filename

     flt = sprintf('*.%s',typ);

     [file,pfad]= uiputfile( flt , 'Save Image ... ' );

     if isequal(file,0) | isequal(pfad,0) | isempty(file)
        return
     end

     file = fullfile(pfad,file);

     %----------------------------------------------
     switch typ
 
       %-----------------------------------------------------------------------
       case  { 'ps'  'eps' }       

          vis = get(ud.line,'visible');
                set(ud.line,'visible','off')

          try
             print(fg,dev,opt{:},file);
          catch
              msg = sprintf('Error call PRINT(''%s'',''%s'').\n%s',dev,file,lasterr);
          end

          set(ud.line,'visible',vis)

       %-----------------------------------------------------------------------
       case  { 'pcx' 'xpm' }

          c = ud.inf.Index;
          m = ud.inf.Colormap;

          if isempty(m) | isempty(c)

             msg = sprintf('%s-Images only supported for indexed Images.',upper(typ));

          else
 
             if prn
                c = c(area{:});
             end

             [nr,c] = rotimg( c , -ud.rot );
             if ud.inv
                m = 1 - m;
             end
             m = rgb2hsv(m);
             m(:,1) = m(:,1) + ud.hue/360;
             m(:,1) = m(:,1) - floor(m(:,1));
             m = hsv2rgb(m);

             if strcmp(typ,'xpm')
                try
                   msg = wrtxpm(c,m,file);
                catch
                   msg = sprintf('Error call WRTXPM(''%s'').\n%s',file.lasterr);
                end
             else
                try
                   imwrite(c,m,file,typ);
                catch
                   msg = sprintf('Error call IMWRITE(''%s'').\n%s',file.lasterr);
                end
             end

          end

       %-----------------------------------------------------------------------
       case  'ppm'

          c = get(ud.img,'cdata');

          if prn
             c = c(area{:});
          end

          try
             msg = wrtppm(c,file);
          catch
             msg = sprintf('Error call IMWRITE(''%s'',''%s'').\n%s',file,typ,lasterr);
          end

       %-----------------------------------------------------------------------
       otherwise

          if strcmp(typ,'jpg') & ~isempty(opt)
             opt = cat(2,{'Quality'},opt);
          end

          c = get(ud.img,'cdata');

          if prn
             c = c(area{:});
          end

          try
             imwrite(c,file,typ,opt{:});
          catch
             msg = sprintf('Error call IMWRITE(''%s'',''%s'').\n%s',file,typ,lasterr);
          end

     end

     if ~isempty(msg)
         msgbox(msg,str,'warn');
     end

%***********************************************       
case 'ROTATE'

     rot = get(h,'userdata');

     if isequal(rot,ud.rot)
        set(h,'checked','on');
        return
     end

     setimg(fig,app,rot,[]);

%***********************************************       
case 'FLP'

     ii = get(h,'userdata');

     ud.flp(ii) = 1 - strcmp(get(h,'checked'),'on');

     setimg(fig,app,[],ud.flp);

%***********************************************       
case 'INV'

    sets = { 'off'  'on' };

      jj = 3 - find( strcmp( get(h,'checked') , sets ) );

    set( h , 'checked' , sets{jj} );

    set( ud.img , 'cdata' , uint8( 255 - double(get(ud.img,'cdata')) ) );

    ud.inv = jj - 1;

    ud.hue = real(ud.hue) - i * imag(ud.hue);

    setappdata(fig,app,ud);

%***********************************************       
case 'SET'

     v = get(h,'userdata');

     if isequal(v,real(ud.hue))
        return
     end

     hue = v + i * imag(ud.hue);

     hueimg(fig,app,hue,[]);

%***********************************************       
case { 'RIGHT'  'LEFT' }

     rot = get(h,'userdata');

     hue = real(ud.hue)+rot + i * imag(ud.hue);

     hueimg(fig,app,hue,[]);

%***********************************************       
case 'LIGHT'

     v = get(h,'userdata');

     if isequal(v,imag(ud.hue))
        return
     end

     hue = real(ud.hue) + i * v;

     hueimg(fig,app,hue,[]);

%***********************************************       
case { 'BRITE'  'DIM' }

     rot = get(h,'userdata');

     hue = real(ud.hue) + i * ( imag(ud.hue) + rot );

     hueimg(fig,app,hue,[]);

%***********************************************       
case 'SMOOTH'

     smt = get(h,'userdata');

     ok = smtimg(fig,app,smt,[]);

     if ok

        p = get(h,'parent');
        set(get(p,'children'),'checked','off');

        set(h,'checked','on');

     end

%***********************************************       
case 'NEW'

     p1  = fileparts(ud.inf.Filename);
     pok = ~isempty(p1);
     if pok
        p0 = cd;
        try
           cd(p1);
        catch
           pok = 0;
        end
     end

     [file,pfad]= uigetfile( '*.*' , 'Select Image ... ' );

     if pok
        cd(p0);
     end

     if isequal(file,0) | isequal(pfad,0) | isempty(file)
        return
     end

     file = fullfile(pfad,file);

     try
        msg = feval(fcn,file,fig);
     catch
        msg = sprintf('Error call %s.\n%s',upper(fcn),lasterr);
     end

     if ~isempty(msg)
         msgbox(msg,'Error','warn');
     end

%***********************************************       
case 'CLOSE'


     close(fig);

end

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function HM = mkmenu(par,group,cnf,tag,fcn);

% MKMENU  Creates UIMenus from Structure
%
% MKMENU( ParentHandle ,  Group , Config , Tag );
%
%
%  FieldValues of Config for UIMenu-Definition:
%
%    { Label  SeparatorOn  CallbackFcn  UserData }
%    {          0 | 1      SubStructure          }
%
%  FieldValues for UIContextMenu-Definition:
%
%    { NaN  MenuStructure  CallbackFcn  UserData }
%


%-----------------------------------------------------------

HM  = []; 

if isempty(cnf)
  return
end


%************************************************************


field = fieldnames(cnf);
field = field(:)';

    n = size(field,2);

%------------------------------------------------------------
% HandleStructure

HM      = cell( 2 , n );
HM(1,:) = field;
HM(2,:) = { {[]} };

HM      = struct( HM{:} );

HC      = zeros( n , 1 );

%------------------------------------------------------------

sets = { 'off'  'on' };

if ~isempty(group)
   tag = cat( 2 , tag , '.' , group );
end

for ii = 1 : n

  cc = getfield( cnf , field{ii} );


  CB    = '';
  usd   = [];

  if chkstr(cc{3},1)
     CB = sprintf('%s(gcbo,''%s'',''%s'',''%s'',1);',fcn,cc{3},group,field{ii});
  end

  if prod(size(cc)) >= 4
     usd = cc{4};
  end

  tag1 = cat( 2 , tag , '.' , field{ii} );

  if ischar(cc{1})

    %----------------------------------
    % UIMenu

    enb = sets{ 1 + (~isempty(rmblank(cc{1},2))) };
    sep = sets{ 1 + isequal(cc{2},1) };
 
    HC(ii) = uimenu( ...
       'parent'          , par      , ...
       'callback'        , CB       , ...
       'userdata'        , usd      , ...
       'label'           , cc{1}    , ...
       'tag'             , tag1     , ...
       'enable'          , enb      , ...
       'separator'       , sep      , ...
       'checked'         , 'off'    , ...
       'Interruptible'   , 'off'    , ...
       'BusyAction'      , 'cancel'  );

     ch = cc{3};

   else

     %----------------------------------
     % UIContextMenu

     HC(ii) = uicontextmenu( ...
       'parent'          , par      , ...
       'callback'        , CB       , ...
       'userdata'        , usd      , ...
       'tag'             , tag1     , ...
       'Interruptible'   , 'off'    , ...
       'BusyAction'      , 'cancel'  );
 
     ch = cc{2};

   end 

  if isstruct(ch)
  % Children

     ch = mkmenu(HC(ii),field{ii},ch,tag,fcn);

  else

     ch = [];

  end

  setappdata( HC(ii) , 'Children' , ch );

  HM = setfield( HM , field{ii} , HC(ii) );

end

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,v0] = structcmp(v0,v1);

% STRUCTCMP  Compares 2 Structures
%
% [Msg,V] = STRUCTCMP( V0 , V1 )
%
%  Check if Fieldnames of V1 contained in V0 and
%    SubStructures too.
%

Msg = '';

if nargin < 2
  Msg = 'Not enough Input Arguments.';
end


nl = char(10);

%----------------------------------------------------
% Check Class

c0 = class(v0);

if ~isequal(class(v1),c0)
   Msg = [ bl  'Variable must be of Class "'  c0  '".' ];
   return
end

%----------------------------------------------------
% Check DimensionLength and Size

   n0 = ndims(v0);
   n1 = ndims(v1);

   ok = ( n0 == n1 );

   if ~ok
      nd  = sprintf('%.0f',n0);
      Msg = [ 'Variable must have DimensionLength ' nd '.' ];
      return
   end

   s0 = size(v0);
   s1 = size(v1);

   ok = ( all( s0 <= s1 )  |  strcmp(c0,'char') );

   if ~ok

       si = sprintf('%.0f x ',s0);
       si = si( 1 : end-3 );

       Msg = [ 'Variable must have Size: [ ' si ' ].' ];

       return

   end

%----------------------------------------------------
   
if ~strcmp(c0,'struct')
   v0 = v1;
   return
end

if ~( ( prod(size(v0)) == 1 ) & ( prod(size(v1)) == 1 ) )
    Msg = 'Structures must be single.';
    return
end

%----------------------------------------------------

f0 = fieldnames(v0);
f1 = fieldnames(v1);

n0 = prod(size(f0));
n1 = prod(size(f1));

ok0 = zeros(n0,1);
ok1 = zeros(n1,1);

for ii = 1 : n1

    ok1(ii) = any( strcmp( f0 , f1{ii} ) );

    if ok1(ii)

       jj = find( strcmp( f0 , f1{ii} ) );

       jj = jj(1);

       ok0(jj) = 1;

       [Msg1,p] = structcmp( getfield(v0,f0{jj}) , getfield(v1,f1{ii}) );

       if isempty(Msg1)
 
           v0 = setfield( v0 , f0{jj} , p );

       else

          Msg = [ Msg  nl(1:(end*(~isempty(Msg))))  ...
                  'Invalid SubStructure for "' f1{ii} '".' nl  Msg1 ];

       end

    end

end

%***************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  ok = chkhndl(h,typ,app);

% CHKHNDL(H,typ,ApplicationData)  Checks, if H is a Handle of specified Type
%

Nin = nargin;

ok = ( isnumeric(h) & ( prod(size(h)) == 1 ) );

if ~ok
   return
end

ok = ishandle(h);
if ~ok | ( Nin < 2 )
   return
end

%-------------------------------------------------------------------------
% Check with Type

ok = strcmp( get(h,'type') , typ  );

if ~ok | ( Nin < 3 )
   return
end

%-------------------------------------------------------------------------
% Check with ApplicationData

ok = isappdata(h,app);


%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,c,t] = recpar(h,typ);

% RECPAR  returns ParentHistory of Handle
%
% [Msg,HandleHist,TypeHist] = RECPAR( Handle , StopType )
%
%  recurse trough the Parents unto ( ParentType == StopType )
%
%  StopType starts with UpperCase, returns History excluding
%     Handle with ( ParentType == StopType )
%
%  default: StopType = 'Root'  (recurse unto 'figure')
%
%  HandleHist(end) == Handle
%    TypeHist(end) == HandleType
%

Msg = '';
 c  = zeros(0,1);
 t  =  cell(0,1);

if nargin < 1
   Msg = 'Input Handle is missing.';
   return
end

if nargin < 2
   typ = 'Root';
end

%-----------------------------------------------

if isempty(h)
   return
end

ok = ( isnumeric(h) &  ( prod(size(h)) == 1 ) );
if ok
   ok = ishandle(h);
end

if ~ok
   Msg = 'First Input must be a Single Handle.';
   return
end

if ~( ischar(typ) & ~isempty(typ) & ...
      ( prod(size(typ)) == size(typ,2) ) )
   Msg = 'Type must be a String';
   return
end

%-----------------------------------------------

c = h;
t = { get(h,'type') };

z = 1;

t0 = lower(typ);

while ~( c(1) == 0 )  &  ( ~strcmp(t{1},t0) | ( z == 1 ) )

   z = z + 1;

   c = cat( 1 ,         get(c(1),'parent')  , c );

   t = cat( 1 , { lower(get(c(1),'type')) } , t );

end


if strcmp( t{1} , t0 )

  n = 1 + strcmp( typ(1) , upper(typ(1)) );

  c = c(n:z);
  t = t(n:z);

else

   Msg = [ 'Handle has no Parents with Type '''  typ '''.' ]; 

end

%-----------------------------------------------
   
%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = rmblank(str,dim,cc)

% RMBLANK  Remove Blanks, NewLines at Begin and End of CharacterArrays
%
% String = RMBLANK( CharArray )
%
% CharArray  2-dimensional CharacterArray
%
% further Options:
%
% String = RMBLANK( CharArray , DIM , CHAR )
%
%  
%  DIM  specifies Dimension to work, 
%       default: 2
%
%    A positive complex Value for DIM, to removes Blanks only from Start,
%    A negative complex Value for DIM, to removes Blanks only from End.
%       
%  CHAR specifies BlankCharacters to remove
%       default:  [ 160  32  13  10  9  0 ];  % [ NBSP Space CR LF TAB ZERO ]
%

  
msg = '';
 nl = char(10);

Nin = nargin;

if Nin < 1
  error('Not enough Input Arguments.')
else
  if ischar(str)
    str = double(str);
  end
  ok = isnumeric(str);
  if ok
    ok = all( ( mod(str(:),1) == 0 )  & ...
              ( str(:) >= 0 ) & isfinite(str(:))  );
  end
  if ~ok
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CharArray must be a String or ASCII-Codes.'];
  end
  if size(str,1)*size(str,2) ~= prod(size(str))
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CharArray must be 2-dimensional.'];
  end     
end

if Nin < 2
  dim = 2;
else
  if ~isnumeric(dim)
    msg = [ msg nl(1:(end*(~isempty(msg)))) ...
            'Input DIM must be numeric.' ];
  elseif ~isempty(dim)
    dim = dim(:);
    if ~all( ( abs(dim) == 1 ) |  ( abs(dim) == 2 ) )
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
             'Values for Input DIM must define 1. or 2. Dimension.' ];
    end
  end 
end

if Nin < 3
  cc = [ 160  32  13  10  9  0 ];  % [ NBSP  Space CR LF TAB ZERO ]
else
  if ischar(cc)
    cc = double(cc);
  end
  ok = isnumeric(cc);
  if ok & ~isempty(cc)
    cc = cc(:)';
    ok = all( ( mod(cc,1) == 0 )  & ...
              ( cc >= 0 ) & isfinite(cc)  );
  end
  if ~ok
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CHAR must be a String or ASCII-Codes.'];
  end
end

if ~isempty(msg)
  error(msg)
end


if isempty(str)
 str = '';
 return
end

if isempty(dim) | isempty(cc)
  str = double(str);
  return
end


  blank  = 0*str;

  for ii = cc
    blank = ( blank | ( str == ii ) );
  end

  si = size(str);

  for ii = 1 : size(dim,1)

    d = dim(ii);

    s = sign(imag(d));  % Remove from wich Side:  1  0  -1 
 
    d = abs(d);

    jj = find( sum(blank,3-d) == si(3-d) );  % Columns with full Blanks

    if ~isempty(jj) 

         p  = [ 3-d  d ];
        str = permute(str,p);

         jj = jj(:)';
         nb = size(jj,2);

        %--------------------------------------------
        % Blank at Begin

        ind = ( 1 : nb );
        jj1 = find( ( ( jj == ind ) & ( s >= 0 ) ) );

        %--------------------------------------------
        % Blank at End

        ind = ind + si(d) - nb;
        jj2 = find( ( ( jj == ind ) & ( s <= 0 ) ) );

        %--------------------------------------------

        str(:,jj([jj1 jj2])) = [];

        str = permute(str,p);

    end
    
  end

  str = char(str);

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [data,nn,dev] = meanind2(dat,n,m,p)

% MEANIND2  2-Dimensional Running Mean of Matrice
%
% [ ZM , NN , STDEV ] = MEANIND2( Z , N , Mode , Potenz );
%
% N       = Length of window, oddinary Integer(s): [ Ny Nx ]
%
% Mode    = Mode of window:  'linear' | 'triangle' | 'cosine' | ...
%                            'gauss'  | 'binomial'
%
% Mode    = positive Numeric for Power of Window: ( -1 .. 1 )
%             Mode == 0    ==>  Impuls
%             Mode == 1    ==>  Triangle
%             Mode == Inf  ==>  Linear (constant)
%
%  defaults:  N = [ 1  1 ]
%             Mode = 'linear'
%             Potenz = 1
%
% NaN's in Z make no Problems.
%
 
Nin  = nargin;
Nout = nargout;

%---------------------------------------------------
% Check Inputs

if Nin < 1
   error('Not enough InputArguments.');
end

si = size(dat);
ps = prod(si);

if ( Nin < 2 ) | ( ps <= 1 )
   data = dat;
   if Nout > 1
    dev =  zeros(si);
     nn = ~isnan(dat);
   end
   return
end

if Nin < 3
   m = 'linear';
end

if Nin < 4
   p = 1;
end

%--------------------------------------
% Check for old version with SET_NAN

if isnumeric(m)
   set_nan = m(1);
   m = 'linear';
else
   set_nan = NaN;
end

%--------------------------------------
% Check N

ok = ( isnumeric(n)  &  ( prod(size(n)) <= 2 ) );
if ok
   ok = all( ( mod(n,2) == 1 ) & ( n > 0 ) );
end
if ~ok
   error('N must be odd Integer(s).');
end

%--------------------------------------

n = n(:)';
n = cat( 2 , n , n );
n = n([1 2]);

try
  c1 = window(n(1),m,p);
  c2 = window(n(2),m,p);
catch
  error(lasterr);
end

if isequal(n,[1 1])
   data = dat;
   if Nout > 1
    dev =  zeros(si);
     nn = ~isnan(dat);
   end
   return
end

%---------------------------------------------------
% Expand dat with -n2  ..  +n2 - Elements

n2  = ( n - 1 ) / 2;

s1 = size(dat,1);
s2 = size(dat,2);

i1 = n2(1);
i2 = n2(2);

o1 = ones(1,i1);
o2 = ones(1,i2);

dat  = dat( cat(2,o1,(1:s1),s1*o1) , cat(2,o2,(1:s2),s2*o2) );

dat( cat(2,(1:i1),(1:i1)+(s1+i1)) , : ) = NaN;
dat( : , cat(2,(1:i2),(1:i2)+(s2+i2)) ) = NaN;

%---------------------------------------------------
% Prepare Variables

quo = ~isnan(dat);

dat( find(~quo) ) = 0;

data = zeros(s1,s2);
quot = zeros(s1,s2);

if Nout >= 2
   nn = zeros(s1,s2);
end

i1 = ( 1 : s1 );
i2 = ( 1 : s2 );

for ii1 = 1 : n(1)
    for ii2 = 1 : n(2)
        data = data + c1(ii1)*c2(ii2)*dat(i1+(ii1-1),i2+(ii2-1));
        quot = quot + c1(ii1)*c2(ii2)*quo(i1+(ii1-1),i2+(ii2-1));
        if Nout > 1
           nn = nn + quo(i1+(ii1-1),i2+(ii2-1));
        end
    end
end

quo = find(~quo(i1+n2(1),i2+n2(2)));

data = data ./ ( quot + ( quot == 0 ) );

%------------------------------------------
if Nout > 2

   dev = zeros(s1,s2);

   for ii1 = 1 : n(1)
       for ii2 = 1 : n(2)
           dev = dev + (dat(i1+(ii1-1),i2+(ii2-1))-data).^2;
       end
   end

   qadd  = 2*(nn==0)+(nn==1);

   dev = sqrt( dev ./ (nn+qadd-1) );

   dev = dev .* (~(qadd==1));

   dev(quo) = NaN;

end
%------------------------------------------

data(quo) = set_nan;

if Nout > 1
   nn(quo) = 0;
end  

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [c,x] = window(n,m,p)

% WINDOW  returns an symetric, normalized window
%
% [C,X] = WINDOW( N , Mode , Potenz )
%
% N       = Length of window, oddinary Integer
% N < 0   = Right part of Window is ZERO
% 
% Mode    = Mode of window: 'linear' | 'triangle' | 'cosine' | ...
%                                      'gauss'    | 'binomial'
%
% Mode    = positive Numeric for Power of Window: ( -1 .. 1 )
%             Mode == 0    ==>  Impuls
%             Mode == 1    ==>  Triangle
%             Mode == Inf  ==>  Linear (constant)
% 
% Potenz  = Potenz for pre-normalized ( max == 1 ) window
%
%
% C       = normalized window, [ 1 by N ], sum(C) == 1 
%
% X       = X-Coordinates, [ 1 by N ], [ -1 .. 1 ]
%
%
%  defaults:  N = 1
%             Mode = 'linear'
%             Potenz = 1
%
% For Demonstration see WINDEMO, type: >> windemo 
%

%---------------------------------------------------
% Check Inputs

msg = '';

nl = char(10);

Nin = nargin;

%---------------------------------------------------
% Check Length

if Nin < 1

   n = 1;

else

  ok = ( isnumeric(n)  &  ( prod(size(n)) == 1 ) );
  if ok
     ok = ( mod(n,2) == 1 );
  end

  if ~ok
     msg = 'N must be an odd Integer.';
  end

end

%---------------------------------------------------
% Check Mode

e = 1;

if Nin < 2

  m = 'l';

else

  ok = ( ischar(m) & ( prod(size(m)) == size(m,2) ) & ~isempty(m) );

  if ~ok
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
      msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                 'Mode must be a String or positive Numeric.' );
  end

end

%---------------------------------------------------
% Check Potenz

if Nin < 3

  p = 1;

else

  ok = ( isnumeric(p)  &  ( prod(size(p)) == 1 ) );
  if ok
     ok = isfinite(p);
  end

  if ~ok
     msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                'Potenz must be a finite Numeric.'    );
  end

end


%---------------------------------------------------

if ~isempty(msg)
   error(msg)
end

%---------------------------------------------------
% Build Window

is2 = ( n < 0 );
n   = abs(n);

if n == 1
   x = 0;
   c = 1;
   return
end

x = ( 1 : n );

switch lower(m(1))

  %-------------------------------------------------
  % Binomial

  case 'b'

    nc = n - 1;

    c = ones(2,nc+1);

    c(1,2:nc+1) = nc - ( 0 : nc-1 );
    c(2,2:nc+1) = ( 1 : nc );

    c = round( cumprod( c(1,:) ./ c(2,:) ) );

  %-------------------------------------------------
  % Cosine

  case 'c'

    c = 1 - cos( 2*pi * x/(n+1) );

  %-------------------------------------------------
  % Gauss

  case 'g'
 
    c = pi;  % exp(1);

    c = c * ( 2 * x / (n+1) - 1 );  % linspace(-c,c,n+2)

    c = exp( (-1) * c.^2 / 2 );

  %-------------------------------------------------
  % Triangle | Potenz

  case { 't'  'p' }

    if isinf(e)

       c = ones(1,n);

    else

       c = 2 * x / (n+1) - 1;

       if e == 0
          c = double( c == 0 );
       elseif e == 1
          c = 1 - abs(c);
       else    
          c = 1 - abs(c) .^ e;
       end

    end

  %-------------------------------------------------
  % linear 
 
  otherwise


    c = ones(1,n);


end

x = 2 * (x-1) / (n-1) - 1;  % linspace(-1,1,n)

%---------------------------------------------------
% Half Window

if is2
   c( (n+1)/2+1 : n ) = 0;
end

%---------------------------------------------------
% Normalize

if ( p == 1 )
  c = c / sum(c);
  return
end

%---------------------------------------------------
% Potenz and Normalize

c = ( c / max(c) ) .^ p;

c = c / sum(c);


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,c,m] = chkcdata(c,m,grey)

% CHKCDATA  Checks Input CData (and ColorMap) for valid ImageData
%
% returns RGB-ColorMatrice or indexed ColorMatrice and ColorMap
%
%-----------------------------------------------------------------------
%     RGB-Input,     RGB-Output: [Msg,CData] = CHKCDATA(CData)
% Indexed-Input,     RGB-Output: [Msg,CData] = CHKCDATA(CData,ColorMap)
%
%     RGB-Input, Indexed-Output: [Msg,CData,ColorMap] = CHKCDATA(CData)
% Indexed-Input, Indexed-Output: [Msg,CData,ColorMap] = CHKCDATA(CData,ColorMap)
%
%-----------------------------------------------------------------------
% RGB-Input
%
%   The RGB-Input CData can be:
%     - a 3-dimensional RGB-ColorMatrice or 
%     - a 2-dimensional GrayScaled Matrice
%    of class:  UINT8, UINT16 or DOUBLE ( 0 <= CData <= 1 ) 
%
%   Multiple Frames (Nf) can be added to the 4. (3.) Dimension of CData.
%
%   The Dimension of the Input RGB-CData can be:
%
%     [ Ny  by  Nx  by  3   by  Nf ] RGB-ColorMatrice
%     [ Ny  by  Nx  by  1   by  Nf ] GrayScaledMatrice
%     [ Ny  by  Nx  by  Nf  by  1  ] GrayScaledMatrice, ( Nf ~= 3 )
%
%-----------------------------------------------------------------------
% Indexed-Input
%
%   The Indexed-Input CData should be a 2-dimensional indexed Matrice,
%    refering to ColorMap, which is a 3-Column RGB-Matrice
%    of class:  UINT8, UINT16 or DOUBLE ( 0 <= ColorMap <= 1 ).
%
%   The Class of the indexed CData can be  UINT8, UINT16 or DOUBLE.
%
%   Multiple Frames (Nf) can be added to the 3. or 4. Dimension of CData.
%
%   The Dimension of the Input indexed CData can be:
%
%     [ Ny  by  Nx  by  Nf  by  1  ]  indexed
%     [ Ny  by  Nx  by  1   by  Nf ]  indexed
%   
%-----------------------------------------------------------------------
% RGB-Output
%
%   The RGB-Output returns a 3-dimensional UINT8-ColorMatrice,
%    contains the RGB-ColorValues ( 0 <= C <= 255 ).
%
%   Multiple Frames (Nf) are added to the 4. Dimension:
%
%     [ Ny  by  Nx  by  3   by  Nf ]  UINT8
%
%-----------------------------------------------------------------------
% Indexed-Output
%
%   The Indexed Output returns a 2-dimensional indexed Matrice,
%    refering to ColorMap, which is a 3-Column RGB-Matrice
%    of class DOUBLE ( 0 <= ColorMap <= 1 ).
%
%   Multiple Frames (Nf) are added to the 3. Dimension:
%
%     [ Ny  by  Nx  by  Nf ]  DOUBLE
%
%-----------------------------------------------------------------------
% NaN-Colors
%
% The optional GREY defines the Color to use for NaN-Colors.
%
%      CHKCDATA( ... , GREY )
%
%  GREY must be a single UINT8-Value.
%
%  Colors, which have in all Values of RGB a NaN will set to GREY.
%  Colors, which have any but not all Values NaN will set to ZERO in NaN-Values.
%
%  default for RGB-Output: GREY = uint8(214)
%
%-----------------------------------------------------------------------
%
% See also: IND2RGB, MAT2IND
%

msg = '';

Nin  = nargin;
Nout = nargout;

if Nin == 0
   c = [];
   m = zeros(0,3);
   return
end

if Nin < 2
   m = [];
end

if Nin < 3
   grey = [];
end

msg = cell(0,1);

%**************************************************
% Check GreyColor

if strcmp(class(m),'uint8') & ( prod(size(m)) == 1 )
   grey = m;
   m    = [];
end

is_grey = ~isempty(grey);

if ~is_grey
    grey = uint8(214);
elseif ~( strcmp(class(grey),'uint8') & ( prod(size(grey)) == 1 ) )
    msg = cat(1,msg,{'Value for GREY must be a single UINT8.'});
end
   
%**************************************************
% Check ColorMap

is_map = ~isempty(m);

cm = class(m);

if any(strcmp(cm,{'uint8' 'uint16'}))
   m = double(m);
   if is_map
      p = 8 * ( 1 + strcmp(cm,'uint16') );
      m = m / ( 2^p - 1 );
   end
   cm = 'double';
end

if ~strcmp(cm,'double')
    msg = cat(1,msg,{'ColorMap must be of class DOUBLE, UINT16 or UINT8.'});
elseif ~is_map
    m = zeros(0,3);
elseif ~( ( ndims(m) == 2 ) & ( size(m,2) == 3 ) & all( abs(m(:)-0.5) <= 0.5 ) )
    msg = cat(1,msg,{'ColorMap must have 3 Columns with Values between 0 and 1.'});
end

%**************************************************
% Check CData

cl = class(c);

if ~any(strcmp(cl,{'uint8' 'uint16' 'double'}))
    msg = cat(1,msg,{'CData must be of class DOUBLE, UINT16 or UINT8.'});
elseif ~isempty(c)
    s3 = size(c,3);
    s4 = size(c,4);
    nd = ndims(c);
    if ( nd > 4 ) | ( ( s4 > 1 ) & ~any( s3 == [ 1  3 ] ) )
        str = sprintf('%s\n%s','CData must have max. 4 Dimensions:', ...
                '2(3|4) Indexed | GrayScale  or  3(4) TrueColor (RGB).');
        msg = cat(1,msg,{str});
    elseif strcmp(cl,'double')
        if ~all( isfinite(c(:)) | isnan(c(:)) )
            msg = cat(1,msg,{'Values of CData must be finite or NaN.'});
        end
    end
end

%**************************************************

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    return
end

msg = '';

if isempty(c) | ( Nout == 1 )
   return
end

%*************************************************
% Input: Indexed CData

s31 = ( s3 == 1 );
s41 = ( s4 == 1 );

if ( is_map & ( s31 | s41 ) )

   is_nan = [];

   if ~strcmp(cl,'double');
      c = double(c) + 1;
   else
      is_nan = find(isnan(c));
   end

   nc = size(m,1);

   c = min(max(round(c),1),nc);

   %----------------------------------------------
   % Output: Indexed CData

   if Nout == 3

      if ~isempty(is_nan)
          if is_grey
             m = cat( 1 , m , double(grey([1 1 1]))/255 );
             c(is_nan) = nc+1;
          else
             c(is_nan) = NaN;
          end
      end

      if ~s41
          c = permute(c,[1 2 4 3]);
      end

      return

   end

   %----------------------------------------------
   % Output: TrueColor CData

   m = uint8(round(m*255));

   if ~isempty(is_nan)
       m  = cat(1,m,grey(ones(1,3)));
       nc = nc+1;
       c(is_nan) = nc;
   end

   if ~s31
       c = permute(c,[1 2 4 3]);
   end

   c = m(cat(3,c,c+nc,c+2*nc));

   return

end

%*************************************************
% Check Class of CData

if strcmp(cl,'uint16')
   c  = uint8( round( 255 * double(c) / (  2^16 - 1 ) ) );
   cl = 'uint8';
elseif strcmp(cl,'double')
   if ~all( ( ( 0 <= c(:) ) & ( c(:) <= 1 ) ) | isnan(c(:)) )
       msg = 'TrueColor CData of class DOUBLE must have Values between 0 and 1.';
       return
   end
end

%*************************************************
% Input: GrayScale CData

if ( s4 == 1 ) & ~any( s3 == [ 1  3 ] )
   c = permute(c,[1 2 4 3]);
   s4 = s3;
   s3 =  1;
end

is_gray = ( s3 == 1 );

if is_gray
   c = c(:,:,ones(1,3),:);
end

%*************************************************
% Output: TrueColor CData

if Nout == 2

   if strcmp(cl,'double')
      is_nan = isnan(c);
      c = uint8(round(c*255));
      if any(is_nan(:))
         if is_gray
              is_nan  = find(is_nan);
            c(is_nan) = grey;
         else
                 nan3 = ( sum(is_nan,3) == 3 );
              is_nan  = find(is_nan);
            c(is_nan) = uint8(0);
            if any(nan3(:))
                 is_nan  = find(nan3(:,:,[1 1 1],:));
               c(is_nan) = grey;
            end    
         end
      end
   end

   return

end

%*************************************************
% Output: Indexed CData

if strcmp(cl,'uint8')
   c = double(c) / 255;
end

s1 = size(c,1);
s2 = size(c,2);

c = permute(c,[1 2 4 3]);
c = reshape(c,s1,s2*s4,3);

[c,m] = mat2ind(c);

is_nan = isnan(m);

m = min(max(m,0),1);

c = reshape(c,s1,s2,s4);

if any(is_nan(:))

   m(find(is_nan)) = 0;

   nan3 = ( sum(is_nan,2) == 3 );

   if any(nan3)

      jj = find(nan3);

      if is_grey
         m(jj,:) = double(grey)/255;
      else
         m(jj,:) = [];
         c(find(c==jj)) = NaN;
         c = c - ( c > jj );
      end
 
   end

end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [c,m] = mat2ind(c)

% MAT2IND  Converts Matrice to Indexed
%
% [ Y , Map ] = MAT2IND( X )
%
%   X   numeric Matrice with max. 3 Dimensions
%
%   Y   Indexed Matrice, refering to Values of Map
%
%  Map  Unique Values of X: [ N  by  size(X,3) ] 
%
% Convert back:
%
%   s1 = size(Y,1); s2 = size(Y,2); s3 = size(Map,2);
%   X  = feval(class(Map),zeros(s1,s2,s3));
%   for ii = 1 : s3
%       X(:,:,ii) = reshape(Map(Y,ii),s1,s2);
%   end
%

if ~( isnumeric(c) & ( ndims(c) <= 3 ) )
    error('Input must be numeric with max. 3 Dimensions.');
elseif isempty(c)
    c = zeros(size(c));
    m = zeros(0,size(c,3));
    return
end


cl = class(c);
if ~strcmp(cl,'double')
    c = double(c);
end

if ~all( isfinite(c(:)) | isnan(c(:)) )
    error('Input must have finite Values or NaN.');
end

s1 = size(c,1);
s2 = size(c,2);
s3 = size(c,3);

n = s1 * s2;

c = reshape(c,n,s3);

si = ( 1 : n )';

for ii = s3 : -1 : 1
    [m,jj] = sort(c(si,ii));
       si  = si(jj);
end

c = c(si,:);

is_nan = find(isnan(c));
if ~isempty(is_nan)
    nn = ceil(max(c(:))) + 1;
    c(is_nan) = nn+1;
    is_nan = 1;
end

c(2:n,:) = c(2:n,:) - c(1:(n-1),:);

if is_nan    
    n = ~( sum( abs(c) < 1e3*eps , 2 ) == s3 );
else
    n = ~( sum(     c  < 1e3*eps , 2 ) == s3 );
end

n(1) = 1;

m = cumsum(c(find(n),:),1);

if is_nan
   m( find( m > nn ) ) = NaN;
end

n = cumsum(n);

c = reshape(n,s1,s2);

c(si) = c;

if ~strcmp(cl,'double')
    m = feval(cl,m);
end

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function Msg = wrtppm(c,m,file);

% WRTPPM  Writes PPM-ImageFile (truecolor) from ImageData
%
% Indexed   ImageData:  Msg = WRTPPM( CData , ColorMap , file );
%  
% TrueColor ImageData:  Msg = WRTPPM( CData , file );
%
% Multiple Frames can be added 
%    to the 3. Dimension of CData in case of indexed   ImageData
%    to the 4. Dimension of CData in case of truecolor ImageData 
%
%    to the 3. or 4. Dimension of CData in case of grayscaled ImageData 
%
% Description:  Writes image in PPM format, a widely-used, simple, 
%      portable, but non-compressed format.  PPM images can be converted
%      to gif, jpg, tif, bmp, pict, and nearly every other image format
%      know to man (or nerd).   Look for the 'netpbm' distribution on 
%      the internet.
%
%
%  see also:  READPPM, WRTXPM, READXPM, IMREAD, IMWRITE
%
%

Nin = nargin;

Msg = '';

nl  = char(10);

grey = uint8(214);     % NaN-Color

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
Msg0 = sprintf('%s: ',upper(fcn));

if Nin < 2
   Msg = sprintf('%sNot enough InputArguments.',Msg0);
   return
end

Msg = cell(0,1);

%**************************************************
% Check File

if ischar(m)
   file = m;
   m    = [];
   Nin  = 3;
elseif Nin < 3
   Msg = cat(1,Msg,{'Input File is missing.'});
end

if Nin == 3 
   if ~( ischar(file) & ~isempty(file) & ...
         ( prod(size(file)) == size(file,2) ) );
       Msg = cat(1,Msg,{'File must be a String.'});
   end
end

%**************************************************
% Check CData and ColorMap

[mc,c] = chkcdata(c,m);

if ~isempty(mc)
    Msg = cat(1,Msg,{mc});
end

if ~isempty(Msg)
    Msg = sprintf('%s\n',Msg{:});
    Msg = sprintf('%sInvalid Inputs.\n%s',Msg0,Msg);
    return
end

Msg = '';

if isempty(c)
   return
end

%**************************************************

w = size(c,2);
h = size(c,1);

c = double(permute(c,[3 2 1 4]));

%**************************************************
% Write file

dt = clock;
dt = datenum(dt(1),dt(2),dt(3),dt(4),dt(5),dt(6));
dt = datestr(dt,0);

fid=fopen(file,'w');

if fid == -1
   Msg = sprintf('%sError open File "%s" for writing.',Msg0,file);
   return
end

fprintf(fid,'P6\n');
fprintf(fid,'# %s\n',dt);
fprintf(fid,'# %s: Matlab %s %s\n',fcn,version,computer);

nf = size(c,4);
if nf > 1
   fprintf(fid,'# %.0f Frames\n',nf);
end

fprintf(fid,'%4d %4d\n',w,h);
fprintf(fid,'%4d\n',255);

fwrite(fid,c,'uchar');

fclose(fid);

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,c,b] = readppm(file);

% READPPM  Reads ImageData from PPM-Format
%
% [ Msg , C ] = READPPM(FileName)
%
%    returns the 3-dimensional UINT8-ColorMatrice C, contains
%     the RGB-ColorValues ( 0 <= C <= 255 ).
%
%    Multiple frames will added to the 4. Dimension of C.
%
%
% [ Msg , C , ColorMap ] = READPPM(FileName)
%
%    returns the 2-dimensional Indexed-ColorMatrice C, 
%     and the corresponding RGB-ColorMap with ColorValues ( 0 <= C <= 1 ).
%
%    Multiple frames will added to the 3. Dimension of C.
%
%
% Description:  PPM format, a widely-used, simple, 
%      portable, but non-compressed format.  PPM images can be converted
%      to gif, jpg, tif, bmp, pict, and nearly every other image format
%      know to man (or nerd).   Look for the 'netpbm' distribution on 
%      the internet.
%
%

Nin  = nargin;
Nout = nargout;

Msg = '';
c   = [];
b   = [];

nl  = char(10);

lmax = 100;    % Max Number of HeaderLines
cm   = '#';    % CommentMark
hd   = { 'P6' 'P3' };   % First Line

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
Msg0 = sprintf('%s: ',upper(fcn));

if Nin < 1
   Msg = sprintf('%sNot enough InputArguments.',Msg0);
   return
end

if ~( ischar(file) & ~isempty(file) & ...
      ( prod(size(file)) == size(file,2) ) );
    Msg = sprintf('%sFile must be a String.',Msg0);
    return
end

%**************************************************
% Open File

fid = fopen(file,'r');

if fid == -1
   Msg = sprintf('%sError open File "%s" for reading.',Msg0,file);
   return
end

%*****************************************************************
% Get Header: hd  Width Height MaxColorValue

str0 = sprintf('%sEnd of File reached during read Header.',Msg0);
str1 = sprintf('%sTo much CommentLines.',Msg0);

cc = 0;

inf = zeros(0,1);

%************************************
% FirstLine

while ( cc < lmax )  &  ( size(inf,1) < 3 )

   m  = 3 - size(inf,1);

   cc = cc + 1;
   bb = fgetl(fid);
   if isequal(bb,-1)
      Msg = str0;
      fclose(fid);
      return
   end
   
   bb = rmblank(bb,2*i);

   %------------------------------------
   % Check FirstLine

   if cc == 1

      bb = rmblank(bb,2*i);
      ok = ( size(bb,2) >= 2 );
      if ok
         ok = any(strcmp(bb(1:2),hd));
      end
      if ~ok
          Msg = sprintf('%sFirst HeaderLine "%s" or "%s" expected in PPM-File "%s".',Msg0,hd{:},file);
          fclose(fid);
          return
      end

      hd = bb(1:2);

      bb = rmblank(bb(3:end),2);

   end

   if ~isempty(bb)
       if ~strcmp(bb(1),cm);
           try
              v = eval(sprintf('[%s]',bb));
           catch
              Msg = sprintf('%sError evaluate HeaderValues: "%s".',Msg0,bb);
              fclose(fid);
              return
           end
           if ~isempty(v)
               ok = isnumeric(v);
               if ok
                  v = v(:);
                  v = v( 1 : min(size(v,1),m) );
                  ok = all( isfinite(v) & ( v > 0 ) & ...
                            ( round(v) == v ) );
               end
               if ~ok
                   Msg = sprintf('%sInvalid HeaderValues: "%s".',Msg0,bb);
                   fclose(fid);
                   return
               end
               inf = cat( 1 , inf , v );
           end
       end
   end

end

if ( size(inf,1) < 3 )
   Msg = str1;
   fclose(fid);
   return
end

%*****************************************************************

sz = inf([1 2]);   % [ Width Height ]
mc = inf(3);       % MaxColorValue

b = 8 * ( 1 + ( mc >=256 ) );

fmt = sprintf('ubit%.0f',b);

%*****************************************************************
% Read Image

is_plain = strcmp(hd,'P3');

if is_plain
   c = fscanf(fid,'%f');
else
   c = fread(fid,fmt);
end

fclose(fid);

if isequal(c,-1)
   Msg = sprintf('%sEnd of File reached during read ImageData.',Msg0);
   return
end

p = 3 * prod(sz);
s = size(c,1);

m = mod(s,p);

if ~( m == 0 )
    wrn = sprintf('Size of %.0f-bit ImageData doesn''t correspond with HeaderSize: %.0f %.0f.', ...
                    b , sz );
    warning(wrn);
    if s < p
       c = cat( 1 , c , zeros(p-s,1) );
    else
       c = cat( 1 , c , zeros(floor(s/p)-s,1) );       
    end
end

n = size(c,1) / p;  % Number of Frames

c = reshape(c,3,sz(1),sz(2),n);

c = permute(c,[3 2 1 4]);
if ~( mc == 255 )
    c = 255 * c / mc;
end

if Nout < 3
   c = uint8(round(c));
   return
end

%*****************************************************************
% Transform to Indexed

c = permute(c,[1 2 4 3]);
c = reshape(c,sz(2),sz(1)*n,3);

[c,b] = mat2ind(c);

c = reshape(c,sz(2),sz(1),n);

b = b / 255;


%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,ccm,ctxt] = wrtxpm(c,m,file);

% WRTXPM  Writes XPM-ImageFile (indexed) from ImageData
%
% Indexed   ImageData:  Msg = WRTXPM( CData , ColorMap , file );
%  
% TrueColor ImageData:  Msg = WRTXPM( CData , file );
%
% CData can be a Indexed ColorMatrice, refering to ColorMap
%   with NaN's for Color None; or a TrueColor-Matrice:
%     - a 3-dimensional RGB-ColorMatrice or 
%     - a 2-dimensional GrayScaled Matrice
%    of class:  UINT8, UINT16 or DOUBLE ( 0 <= CData <= 1 ) 
%
% Note: The XPM-Format is recommended for small indexed Images like Icons etc.
%       The conversion of truecolor to indexed CData can spend a lot of Time !!!
%
%   Multiple Frames can be added to the 3. (4.) Dimension of CData.
%
%
% ColorMap  [ R G B ], 0 <= ColorMap <= 1
%
% [Msg,CharacterColorMap,ColorMapText] = WRTXPM( ... )
%
%    returns the Colormap, translated into Characters (without Color "None").
% 
%    CharacterColorMap = [ size(ColorMap,1) by CharacterPerPixel ]
%    ColorMapText      = [ size(ColorMap,1) by NString           ]
%
% WRTXPM uses follwing special Characters:
%
%     '.'   Color NaN    None        comes first in List !!!
%     'o'   Color White  #FFFFFF
%     '#'   Color Black  #000000
%
%   The Characters ' "%\' will not used.
%
% 
%  see also:  READXPM, WRTPPM, READPPM, IMREAD, IMWRITE
%


Nin = nargin;

Msg  = '';

ccm  = [];
ctxt = '';

nl   = char(10);

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
Msg0 = sprintf('%s: ',upper(fcn));

%*****************************************************************

ColorDepth = 2;   % [ 1  |  2  |  4 ]

ColorScale = 16^ColorDepth - 1;

% Special Characters

ColorNaN   = '.';
ColorWhite = 'o';
ColorBlack = '#';

% Forbidden Characters

ColorForbid = ' "%\';

% Comment Characters

cm0 = '/*';
cm1 = '*/';

%*****************************************************************
% Check Inputs

if Nin < 2
  Msg = sprintf('%sNot enough InputArguments.',Msg0);
  return
end

Msg = cell(0,1);

%**************************************************
% Check File

if ischar(m)
   file = m;
   m    = [];
   Nin  = 3;
elseif Nin < 3
   Msg = cat(1,Msg,{'Input File is missing.'});
end

if Nin == 3 
   if ~( ischar(file) & ~isempty(file) & ...
         ( prod(size(file)) == size(file,2) ) );
       Msg = cat(1,Msg,{'File must be a String.'});
   end
end

%**************************************************
% Check CData and ColorMap

[mc,c,m] = chkcdata(c,m);

if ~isempty(mc)
    Msg = cat(1,Msg,{mc});
end

if ~isempty(Msg)
    Msg = sprintf('%s\n',Msg{:});
    Msg = sprintf('%sInvalid Inputs.\n%s',Msg0,Msg);
    return
end

Msg = '';

if isempty(c)
   return
end

%-------------------------------------------------------
% Check for NaN-Colors in ColorMatrice

c(find(isnan(c))) = 0;

is_nan = any( c(:) == 0 );

if is_nan
   c = c + 1;
   m = cat( 1 , NaN*ones(1,3) , m );
end

nc = size(m,1);


%-------------------------------------------------------
% ColorCharacters

%                 a ..  z       A .. Z       0 .. 9

ccm = cat( 2 , ( 97 : 122 ) , ( 65 : 90 ) , ( 48 : 57 ) , ...
               ( 33 :  47 ) , ( 58 : 64 ) , ( 91 : 94 ) , ...
               (123 : 126 ) );
  
% abcdefghijklmnopqrstuvwxyz ABCDEFGHIJKLMNOPQRSTUVWXYZ 0123456789
% !"#$%&'()*+,-./   :;<=>?@   [\]^   {|}~

ccm = char(ccm);

ccm = ccm(:);


% Remove Special and Forbidden Colors

for cf = [ ColorNaN ColorWhite ColorBlack ColorForbid ]

  ccm( find( double(ccm) == double(cf) ) ) = [];

end

nini = size(ccm,1);

ccm = cat( 1 , ColorNaN  , ColorWhite , ColorBlack , ccm );

%--------------------------------
% Check for Special Color

in  = find( isnan(sum(m,2)) );
isn = ~isempty(in);

m(in,:) = NaN;
 
iw  = find( sum(m,2) == 3 );
isw = ~isempty(iw);

ik  = find( sum(m,2) == 0 );
isk = ~isempty(ik);

isp = [ prod(size(in)) prod(size(iw)) prod(size(ik)) ];

sp  = sum(isp);

%--------------------------------
% Get Char_Per_Pixel

cpp = 1;
while  nc-sp > nini^cpp 
  cpp = cpp + 1; 
end

%--------------------------------
% ColorCharacters

ok     = ones(nc,1);
ok(in) = 0;
ok(iw) = 0;
ok(ik) = 0;
ok     = find(ok);

ii     = zeros(nc,cpp);

ii(ok,:) = comb(nini,cpp,nc-sp) + 3;

ii(in,:) = 1;
ii(iw,:) = 2;
ii(ik,:) = 3;

ccm     = ccm(ii);


%--------------------------------------------------------  
% Transform ColorMap to Character

if isn
   m(in,:) = 1;   %  NaN  --->  1   for DEC2HEX
end

m = round([ ones(1,3) ; m ] * ColorScale );
m = m(:);

m = sprintf(['%0' int2str(ColorDepth) 'X'],m);  

% cmap = [ 3*(nc+1) by ColorDepth ]
% cmap = [ R1 ; R2 ; ... ; G1 ; G2 ; ... ; B1 ; ... BN ] 

m = permute(m,[2 1]);  

% [ R1  .. RN , G1 .. GN , B1 .. BN ]
m = reshape(m,ColorDepth,nc+1,3);
m = permute(m,[2 1 3]);

m = reshape(m,nc+1,3*ColorDepth);

m = m(2:nc+1,:);   % Remove First Color, added before dec2hex

 
i1   = ones(nc,1);

c1 = double('"');
c2 = [  9       double('c #') ];
c3 = [ double('",')  10        ];

m = cat( 2 , char(i1*c1)  , ...
             ccm          , ...
             char(i1*c2)  , ...
             m            , ...
             char(i1*c3)         );      

% NaN-Color  -->  None

if isn

   cn = 'None';

   cn = cat( 2 , cn , char(32*ones(1,3*ColorDepth-size(cn,2)+1)) );

   i1 = ones( prod(size(in)) , 1 );

   cn = cat( 2 , char(i1*c1)   , ...
                 ccm(in,:)     , ...
                 char(i1*c2(1:end-1))  , ...
                 cn(i1,:)              , ...
                 char(i1*c3)                 );      

   m(in,:) = cn;

end
 
ctxt = m;

m = permute(m,[2 1]);
m = permute(m(:),[2 1]);               


%****************************************************
% Transform Image to Character

% Matrice for CharacterIndize

nf = size(c,3);   % Number of Frames

sc = size(c);

%-----------------------------------------------
% Multiple Frames

if nf > 1

   c = reshape( permute( c , [ 2 1 3 ] ) , sc(2) , sc(1)*nf );

   c = permute( c , [ 2  1 ] );

end

%-----------------------------------------------

cc = zeros(sc(1)*nf,sc(2)*cpp);

ccm = permute(ccm,[2 1]);   % [ cpp by  nc(+1) ] 

for ii = 1 : nc

 [i1,i2] = find( c == ii );

 if ~isempty(i1)

    jj = ones(size(i1,1),1);
    ip = ones(1,cpp);


   % Indize in Row of cc
    ic = ((i2-1)*cpp) * ip  +  jj * ( 1 : cpp );

   % Indize in cc
    ic = (ic-1)*sc(1)*nf + i1 * ip;

   % Index of Color ii in ccm
    ind = (ii-1)*cpp + (1:cpp);

    cc(ic) = jj * ind;

 end
end


% Matrice of ColorCharacters

cc = ccm(cc);


if ~isequal( size(cc) , [ sc(1)*nf  sc(2)*cpp ] )
  cc = permute(cc,[2 1]);
  if ~isequal( size(cc) , [ sc(1)*nf  sc(2)*cpp ] )
    Msg = [ Msg0 'Can''t determine CharacterMatrice.' nl ...
                 'File: '  file  '  is empty.'               ];
    fclose(fid);
    return
  end
end

i1   = ones(sc(1)*nf,1);

cc = cat( 2 , char(i1*c1)  , ...
              cc           , ...
              char(i1*c3)         );      

%-------------------------------------
% CharacterColorMap

ccm = permute(ccm,[2 1]);

if is_nan

   % Remove added NaN-Color

    ccm(1,:) = [];  
   ctxt(1,:) = [];

end

%-------------------------------------
% Add FrameNumber as Comment

if nf > 1

   ind = ones( (sc(1)+1)*nf , 1 );

    ic = ( 1 : sc(1)+1 : (sc(1)+1)*nf-sc(1) );  % Index for added Commentary

   ind( ic + 1 ) = 0;

   ind = cumsum(ind);

   cc = cc(ind,:);

   cm = cell(nf,2);    % [ Dummy  True ]
   s2 = size(cc,2)-1;  % Without NewLine

   for ii = 1 : nf

       cm{ii,1} = sprintf('%%%3.3d',ii);
       cm{ii,1} = cat( 2 , cm{ii,1} , char( 32*ones(1,s2-size(cm{ii,1},2))) );

       cm{ii,2} = sprintf( '%s Frame %3.0f of %.0f %s' , cm0 , ii , nf ,cm1 );      

       cc(ic(ii),1:s2) = cm{ii,1};

   end

end

%-------------------------------------
% Build String

cc = permute(cc,[2 1]);
cc = permute(cc(:),[2 1]);               

%-------------------------------------
% Replace DummyComment for FrameNumber

if nf > 1

   for ii = 1 : nf

       cc = strrep(cc,cm{ii,1},cm{ii,2});

   end

end

%****************************************************
% Write File

% File starts with:

% /* XPM */
% /* DATE */
% /* MFile VERSION Computer */
% static char *go_left[] = {
% /* width height Ncolor chars_per_pixel Nframe*/
% "    70    70   7       1              1",
% /* colors */
% ". c None",
% "a c #003910",
% /* pixels */

%-------------------------------------

dt = clock;
dt = datenum(dt(1),dt(2),dt(3),dt(4),dt(5),dt(6));
dt = datestr(dt,0);

fid = fopen(file,'wt');

if fid == -1
  Msg = [ Msg0 'Can''t open file: '  file '   for writing.' ];
  return
end

[pfad,name]=fileparts(file);


fprintf(fid,'%s %s %s\n',cm0,'XPM',cm1);
fprintf(fid,'%s %s %s\n',cm0,dt,cm1);
fprintf(fid,'%s %s: Matlab %s %s %s\n',cm0,fcn,version,computer,cm1);
fprintf(fid,'static char *%s[] = {\n',name);
fprintf(fid,'%s width height Ncolor chars_per_pixel Nframe %s\n',cm0,cm1);
fprintf(fid,'"%s ",\n',sprintf('   %3.0f',[ sc([2 1])  nc  cpp  nf]));
fprintf(fid,'%s %s %s\n',cm0,'colors',cm1);
fprintf(fid,'%s',m);
fprintf(fid,'%s %s %s\n',cm0,'pixels',cm1);

fprintf( fid ,  cc(1:end-2) );

fprintf(fid,'\n};\n');

fclose(fid);

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function c = comb(n,p,l)

% C = COMB(N,P,L);
%
% Combinations of N Elements to P
%
% Returns C with Max. Length L
%

 c = zeros(n,n^(p-1),p);

 c(:,1,1) = (1:n)';
 
 for ii = 2 : p
    z = n^(ii-2);  % RowNumber of existing Combinations
                   %  c(:,1:z,:)
   jj = ( 1 : z*n );

   jj = jj - z*(ceil(jj/z)-1); % [ 1 .. z   1 .. z   1 .. z  ... ]

   c( : , 1 : (z*n) , 1 : (ii-1) ) = c(:,jj,1:(ii-1));

   jj                  =  ones(n,z*n);
   jj(1,:)             = 0;
   jj(1,z*(0:(n-1))+1) = 1;
   jj(1,:)             = cumsum(jj(1,:),2);
   jj                  = cumsum(jj,1);
   jj                  = jj - n * ( jj > n );

   c(:,1:z*n,ii)       = jj;

 end

 c = permute(c,[3 2 1]);  % [ p   by  n^(p-1)  by  n ]
 c = reshape(c,p,n^p);    % [ p   by  n^p ]
 c = permute(c,[2 1]);    % [ n^p by  p   ]

 c = c( 1:min(n^p,l) , [ 1  (p:-1:2) ] );

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,c,cmap,ccm] = readxpm(file);


% READXPM  Reads Image Data from XPM-Format
%
% [ Msg , C , Colormap , CharacterColorMap ] = READXPM(FileName)
%
%  returns the Indexed Colormatrice C and the corresponding ColorMap.
%
%
% [ Msg , C ] = READXPM(FileName)
%
%   returns the 3-dimensional True-ColorMatrice C, contains
%    the ColorValues  ( 0 <= C <= 1 );
%
% READXPM supports multiple Frames, concatenated in the XPM-File.
%
%  The Frames will added to the 3. Dimension of C, in case of indexed Colors,
%                        or the 4. Dimension of C, in case of RGB-True-Colors.
%
%
% Note:  READXPM  requires Hexadecimal ColorDefinitions,
%                  like:    '#ffffff'
%
%                 or X11-ColorNames,
%                  like:    'white'
%
%  The X11-ColorNames are defined by the function RGB at the end of READXPM
%
%  see also: WRTXPM, IMREAD, IMWRITE
%
%----------------------------------------------------------------------------
%  example for an XPM-File:
%
%/* XPM */
%/* MagnifyGlass */
%static char*glass[]={
%/* width hight num_colors chars_per_pixel */
%"16 16 4 1",
%/* colors */
%". c None",
%"# c #000000",
%"o c #ffffff",
%/* pixels */
%"....###.........",
%"..##o.o##.......",
%".#.o...o.#......",
%".#o..#..o#......",
%"#o...#...o#.....",
%"#..#####..#.....",
%"#o...#...o#.....",
%".#o..#..o#......",
%".#.o...o.#o.....",
%"..##o.o####o....",
%"....###.o###o...",
%".........o###o..",
%"..........o###o.",
%"...........o###.",
%"............o#o.",
%"................"};
%


Msg    = '';
c      = [];
cmap   = zeros(0,3);
ccm    = [];

nl  = char(10);

Msg0 = 'READXPM: ';

Nout = nargout;


if nargin < 1
  Msg = [ Msg0  'Input FileName is undefined.' ];
  return
end

if ~ischar(file) | isempty(file) | ( prod(size(file)) ~= size(file,2) )
  Msg = [ Msg0  'FileName must be a string.' ];
  return
end


% File starts with:

% /* XPM */
% static char *go_left[] = {
% /* width height num_colors chars_per_pixel */
% "    70    70        7            1",
% /* colors */
% ". c #c0c0c0",
% 

% ColorIdentifer
cid = 'c ';

%********************************************************************
% Open & Read File

fid = fopen(file,'r');

if fid == -1
 file1 = cat( 2 , file(1:(end-strcmp(file(end),'.'))) , '.xpm' );
 fid = fopen(file1,'r');
 if fid == -1
   Msg = [ Msg0  'Can''t open file: '  file ];
   return
 end
 file = file1;
end

%-------------------------------------------------
% Get First Line and Check Header

bb = fgetl(fid);

if ~ischar(bb)
 if  bb == -1
   Msg = [ Msg0  'Empty file: '  file ];
   return
 end
end
 
header = '/* XPM */';
nh     = size(header,2);
nb     = size(bb,2);
n      = min(nh,nb);

if ~strcmp(bb(1:n),header(1:n))
 Msg = [ Msg0  'Invalid Header, '''  header  ''' expected.' ];
 return
end
  
%-------------------------------------------------
% Read Following Data

bb = fread(fid,'char');

fclose(fid);

%********************************************************************

bb = bb(:)';

%-----------------------------------------------------------------
% Remove CR

bb( find( bb == 13 ) ) = [];

nb = size(bb,2);

%-----------------------------------------------------------------
% Remove Blanks, leading NewLine

jj = 1;

while ~isempty(jj)

   jj = find( ( bb(1:end-1) == 32 ) & ( bb(2:end) == 10 ) );

   bb(jj) = [];

end

nb = size(bb,2);

%-----------------------------------------------------------------
% Get Start-EndIndex of DataLines: " .. "

i01 = find( bb == double('"') );

ni = prod(size(i01));

ok = zeros(1,ni);

%-------------------------------------------------------------
% Check Previous Character of Start-EndIndex: NewLine

is1 = ( i01(1) == 1 );

if is1
   ok(1) = ok(1) + 1;
end

ind = ( (1+is1) : ni );

ok(ind) = ( bb(i01(ind)-1) == 10 );  % Previous is NewLine

%-------------------------------------------------------------
% Check Following Characters of Start-EndIndex: NewLine or ',\n'

% Following is NewLine

ind = ( 1 : (ni-1) );

ok(ind) = ok(ind) + ( bb(i01(ind)+1) == 10 ); 


% Following is ',\n'

is2 = ( i01(ni-1) == nb-1 );

ind = ( 1 : (ni-1-is2) );

ok(ind) = ok(ind) + ( ( bb(i01(ind)+1) == double(',') ) & ( bb(i01(ind)+2) == 10 ) );  

% The should must be allways ok

ok(ni) = 1;  

%-------------------------------------------------------------
% Only ONE ok is allowed

ok = ( ok == 1 );

i01 = i01( find( ok == 1 ) );

%-------------------------------------------------------------

if ( mod(i01,2) ~= 0 )  |  isempty(i01)
 Msg = [ Msg0  'Syntax Error.' ];
 return
end


i01 = reshape(i01,2,length(i01)/2)';

i01(:,1) = i01(:,1)+1;  % Start  of '"'
i01(:,2) = i01(:,2)-1;  % End    of '"'

%********************************************************************

%-----------------------------------------------------------------
% Read Initialisation from 1. Line
% /* width height num_colors chars_per_pixel */

msg = '';
try
  eval(['ini = ['  char(bb(i01(1,1):i01(1,2)))  '];'] );
catch
  msg = lasterr;
end

ok = isempty(msg);
if ok 
   ini = ini(:)';
   ok = ( isnumeric(ini) & ( size(ini,2) >= 4 ) & ~isempty(ini) );
   if ok
      ok = all( isfinite(ini) & ( ini > 0 ) & ( mod(ini,1) == 0 ) );
   end
end


if ~ok
  Msg = [ Msg0 'Invalid Initalisation for [ Width Hight NumColors CharPerPixel ]. ' ];
  return
end 

si  = ini([2 1]);  % [ Height  Width ]
nc  = ini(3);      % ColorNumber
cpp = ini(4);      % CharPerPixel

%********************************************************************
% Check Color- and PixelNumber

npr = size(i01,1) - 1 - nc;   % Number pf Pixel-Rows 

if npr < si(1)
  Msg = [ Msg0 'Invalid Number of Colors and Pixel. ' ];
  return
end 

nf = floor( npr / si(1) );  % Number of Frames

if nf*si(1) < npr

   ww = warnstat;

   warning('on');

   fprintf(nl)
   warning(['Incomplete Frames in XPM-File: '  file  ]);
   fprintf(nl)

   warning(ww);

   i01 = i01( 1 : 1+nc+nf*si(1) , : );
end

npr = nf * si(1);

ind = 1 + nc + ( 1 : npr );

if ~all( i01(ind,2)-i01(ind,1)+1 == cpp*si(2) )
  Msg = [ Msg0 'Invalid Number of Pixel per Row. ' ];
  return
end 


p = 10 .^ ( 3 * (0:cpp-1) );  % Potenz for add multiple PixelCharacter


%-----------------------------------------------------------------
% Get Colors

if ~all( i01(2:nc+1,2)-i01(2:nc+1,1)+1 >= cpp+1 )
  Msg = [ Msg0 'Error reading Colormap.' ];
  return
end

cmap = NaN*zeros(nc,3);
cini = NaN*zeros(nc,2);   % [ IniValue  CmapIndex ]

crgb = rgb;  % { ColorName    ColorHex     ColorDec        }
             %   'white'      '#FFFFFF'   [ 255 255 255 ] 

crgb(:,1) = lower(crgb(:,1));

% CharacterMap

ccm = bb( ones(nc,1)*(1:cpp) + i01((1:nc)+1,ones(1,cpp)) - 1 );

if isequal(size(ccm),[cpp nc])
   ccm = permute(ccm,[2 1]);
end

cini(:,1) = sum( ccm.*p(ones(nc,1),:) , 2);
 
for ii = 1 : nc

  ind = ( i01(ii+1,1)+cpp : i01(ii+1,2) );
 

  % Search for HexColor
  jj = find( bb(ind) == double('#') );

  if isempty(jj)
 
     cc = cat(2,32,bb(ind),32);    % ColorString, Blanks added
    
     % Remove Single Characters

      ch = ~( cc == 32 );  % NonBlanks
     dch = diff(ch);
     is1 = ( ( dch(1:end-1) ==  1 )  &  ... 
             ( dch(2:end-0) == -1 )         );
     cc(find(is1)+1) = [];
        
     % Remove Blanks from left and right
     jj = find( ( cc == 32 )  |   ( cc == 9 ) );
     if ~isempty(jj) 
        jj1 = find( jj == ( 1 : length(jj) ) );
        jj2 = find( jj == ( size(cc,2)-length(jj)+1 : size(cc,2) ) );
        cc(jj([jj1 jj2])) = []; 
     end

     % Check ColorName  or 'none'
     kk  = find( strcmp( lower(char(cc)) , crgb(:,1) ) );

     if isempty(kk)  &  isempty( findstr( lower(char(cc)) , 'none' ) )

       Msg = [ Msg nl(1:(end*(~isempty(Msg))))  ...
               'Color "'  char(bb(ind)) '" not supported.' ];

     elseif ~isempty(kk)

        kk = kk(1);

        cc = crgb{kk,2}(2:end);
        s  = size(cc,2)/3;

        cmax = hex2dec(char(double('f')*ones(1,s)));
        for jj = 1 : 3
          cmap(ii,jj) = hex2dec( char( cc( (jj-1)*s+1 : jj*s ) )  ) / cmax;       
        end

        cini(ii,2) = ii;

     end


  else

     cc = bb(ind(jj)+1:ind(end));
     
     % Remove Blanks from left and right
     jj = find( ( cc == 32 )  |   ( cc == 9 ) );
     if ~isempty(jj) 
        jj1 = find( jj == ( 1 : length(jj) ) );
        jj2 = find( jj == ( size(cc,2)-length(jj)+1 : size(cc,2) ) );
        cc(jj([jj1 jj2])) = []; 
     end

     % Step for single Color
     s = size(cc,2) / 3;

     if ~any( s == [ 1  2  4 ] )
       Msg = [ Msg0 'Invalid Color. ' ];
       return
     end

     cmax = hex2dec(char(double('f')*ones(1,s)));
     for jj = 1 : 3
       cmap(ii,jj) = hex2dec( char( cc( (jj-1)*s+1 : jj*s ) )  ) / cmax;       
     end

     cini(ii,2) = ii;

  end

end

if ~isempty(Msg)
  Msg = [ Msg0 'Error reading Colormap.'  nl  Msg  ];
  return
end

jj = find(~isnan(cini(:,2)));
cini = cini(jj,:);
cmap = cmap(jj,:);
ccm  =  ccm(jj,:);

cini(:,2) = cumsum(ones(size(cini,1),1),1);


%-----------------------------------------------------------------
% Get Pixel

ind = 1 + nc + ( 1 : npr );

ii      = ones(npr,cpp*si(2));
ii(:,1) = i01(ind,1);

ii = cumsum(ii,2);

c = bb(ii);

c = permute(c,[2 1]);            % [ cpp*si(2)    by npr ]
c = reshape(c,cpp,si(2),npr);    % [ cpp by si(2) by npr ]
c = permute(c,[3 2 1]);          % [ npr by si(2) by cpp ]

p = permute(p(:),[2 3 1]);          % [ 1 by 1 by cpp ]

c = sum( c .* p(ones(1,npr),ones(1,si(2)),:) , 3 );

c0 = c;
c  = NaN*c;

for ii = 1 : size(cini,1)

  c(find(c0==cini(ii,1))) = cini(ii,2);

end

%-----------------------------------------------
% RGB-Image

if Nout < 3

 cmap = cat( 1, cmap , NaN*ones(1,3) );

 n = size(cmap,1);

 c(find(isnan(c))) = n;

 c = cmap(cat(3,c,c+n,c+2*n));

end

%-----------------------------------------------
% Multiple Frames

if nf > 1

   sc = size(c);

   nd = size(sc,2);

   perm = cat( 2 , ( 2 : nd ) , 1 );
   
   sc   = cat( 2 , sc(2:nd) , si(1) , nf );

   c = reshape( permute( c , perm ) , sc );

   perm = cat( 2 , nd , ( 1 : nd-1 ) , nd+1 );

   c = permute( c , perm );

end

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  c = rgb

% RGB  returns X-Windows RGB-Colors 
%
% RGB is a 3-Column CellArray: 
%             { ColorName   HexString   DecColors }
%  example:     'white'      '#FFFFFF'   [ 255 255 255 ] 

 c = { ... 

  'snow'                    '#FFFAFA'    [ 255 250 250 ]
  'ghost white'             '#F8F8FF'    [ 248 248 255 ]
  'GhostWhite'              '#F8F8FF'    [ 248 248 255 ]
  'white smoke'             '#F5F5F5'    [ 245 245 245 ]
  'WhiteSmoke'              '#F5F5F5'    [ 245 245 245 ]
  'gainsboro'               '#DCDCDC'    [ 220 220 220 ]
  'floral white'            '#FFFAF0'    [ 255 250 240 ]
  'FloralWhite'             '#FFFAF0'    [ 255 250 240 ]
  'old lace'                '#FDF5E6'    [ 253 245 230 ]
  'OldLace'                 '#FDF5E6'    [ 253 245 230 ]
  'linen'                   '#FAF0E6'    [ 250 240 230 ]
  'antique white'           '#FAEBD7'    [ 250 235 215 ]
  'AntiqueWhite'            '#FAEBD7'    [ 250 235 215 ]
  'papaya whip'             '#FFEFD5'    [ 255 239 213 ]
  'PapayaWhip'              '#FFEFD5'    [ 255 239 213 ]
  'blanched almond'         '#FFEBCD'    [ 255 235 205 ]
  'BlanchedAlmond'          '#FFEBCD'    [ 255 235 205 ]
  'bisque'                  '#FFE4C4'    [ 255 228 196 ]
  'peach puff'              '#FFDAB9'    [ 255 218 185 ]
  'PeachPuff'               '#FFDAB9'    [ 255 218 185 ]
  'navajo white'            '#FFDEAD'    [ 255 222 173 ]
  'NavajoWhite'             '#FFDEAD'    [ 255 222 173 ]
  'moccasin'                '#FFE4B5'    [ 255 228 181 ]
  'cornsilk'                '#FFF8DC'    [ 255 248 220 ]
  'ivory'                   '#FFFFF0'    [ 255 255 240 ]
  'lemon chiffon'           '#FFFACD'    [ 255 250 205 ]
  'LemonChiffon'            '#FFFACD'    [ 255 250 205 ]
  'seashell'                '#FFF5EE'    [ 255 245 238 ]
  'honeydew'                '#F0FFF0'    [ 240 255 240 ]
  'mint cream'              '#F5FFFA'    [ 245 255 250 ]
  'MintCream'               '#F5FFFA'    [ 245 255 250 ]
  'azure'                   '#F0FFFF'    [ 240 255 255 ]
  'alice blue'              '#F0F8FF'    [ 240 248 255 ]
  'AliceBlue'               '#F0F8FF'    [ 240 248 255 ]
  'lavender'                '#E6E6FA'    [ 230 230 250 ]
  'lavender blush'          '#FFF0F5'    [ 255 240 245 ]
  'LavenderBlush'           '#FFF0F5'    [ 255 240 245 ]
  'misty rose'              '#FFE4E1'    [ 255 228 225 ]
  'MistyRose'               '#FFE4E1'    [ 255 228 225 ]
  'white'                   '#FFFFFF'    [ 255 255 255 ]
  'black'                   '#000000'    [   0   0   0 ]
  'dark slate gray'         '#2F4F4F'    [  47  79  79 ]
  'DarkSlateGray'           '#2F4F4F'    [  47  79  79 ]
  'dark slate grey'         '#2F4F4F'    [  47  79  79 ]
  'DarkSlateGrey'           '#2F4F4F'    [  47  79  79 ]
  'dim gray'                '#696969'    [ 105 105 105 ]
  'DimGray'                 '#696969'    [ 105 105 105 ]
  'dim grey'                '#696969'    [ 105 105 105 ]
  'DimGrey'                 '#696969'    [ 105 105 105 ]
  'slate gray'              '#708090'    [ 112 128 144 ]
  'SlateGray'               '#708090'    [ 112 128 144 ]
  'slate grey'              '#708090'    [ 112 128 144 ]
  'SlateGrey'               '#708090'    [ 112 128 144 ]
  'light slate gray'        '#778899'    [ 119 136 153 ]
  'LightSlateGray'          '#778899'    [ 119 136 153 ]
  'light slate grey'        '#778899'    [ 119 136 153 ]
  'LightSlateGrey'          '#778899'    [ 119 136 153 ]
  'gray'                    '#BEBEBE'    [ 190 190 190 ]
  'grey'                    '#BEBEBE'    [ 190 190 190 ]
  'light grey'              '#D3D3D3'    [ 211 211 211 ]
  'LightGrey'               '#D3D3D3'    [ 211 211 211 ]
  'light gray'              '#D3D3D3'    [ 211 211 211 ]
  'LightGray'               '#D3D3D3'    [ 211 211 211 ]
  'midnight blue'           '#191970'    [  25  25 112 ]
  'MidnightBlue'            '#191970'    [  25  25 112 ]
  'navy'                    '#000080'    [   0   0 128 ]
  'navy blue'               '#000080'    [   0   0 128 ]
  'NavyBlue'                '#000080'    [   0   0 128 ]
  'cornflower blue'         '#6495ED'    [ 100 149 237 ]
  'CornflowerBlue'          '#6495ED'    [ 100 149 237 ]
  'dark slate blue'         '#483D8B'    [  72  61 139 ]
  'DarkSlateBlue'           '#483D8B'    [  72  61 139 ]
  'slate blue'              '#6A5ACD'    [ 106  90 205 ]
  'SlateBlue'               '#6A5ACD'    [ 106  90 205 ]
  'medium slate blue'       '#7B68EE'    [ 123 104 238 ]
  'MediumSlateBlue'         '#7B68EE'    [ 123 104 238 ]
  'light slate blue'        '#8470FF'    [ 132 112 255 ]
  'LightSlateBlue'          '#8470FF'    [ 132 112 255 ]
  'medium blue'             '#0000CD'    [   0   0 205 ]
  'MediumBlue'              '#0000CD'    [   0   0 205 ]
  'royal blue'              '#4169E1'    [  65 105 225 ]
  'RoyalBlue'               '#4169E1'    [  65 105 225 ]
  'blue'                    '#0000FF'    [   0   0 255 ]
  'dodger blue'             '#1E90FF'    [  30 144 255 ]
  'DodgerBlue'              '#1E90FF'    [  30 144 255 ]
  'deep sky blue'           '#00BFFF'    [   0 191 255 ]
  'DeepSkyBlue'             '#00BFFF'    [   0 191 255 ]
  'sky blue'                '#87CEEB'    [ 135 206 235 ]
  'SkyBlue'                 '#87CEEB'    [ 135 206 235 ]
  'light sky blue'          '#87CEFA'    [ 135 206 250 ]
  'LightSkyBlue'            '#87CEFA'    [ 135 206 250 ]
  'steel blue'              '#4682B4'    [  70 130 180 ]
  'SteelBlue'               '#4682B4'    [  70 130 180 ]
  'light steel blue'        '#B0C4DE'    [ 176 196 222 ]
  'LightSteelBlue'          '#B0C4DE'    [ 176 196 222 ]
  'light blue'              '#ADD8E6'    [ 173 216 230 ]
  'LightBlue'               '#ADD8E6'    [ 173 216 230 ]
  'powder blue'             '#B0E0E6'    [ 176 224 230 ]
  'PowderBlue'              '#B0E0E6'    [ 176 224 230 ]
  'pale turquoise'          '#AFEEEE'    [ 175 238 238 ]
  'PaleTurquoise'           '#AFEEEE'    [ 175 238 238 ]
  'dark turquoise'          '#00CED1'    [   0 206 209 ]
  'DarkTurquoise'           '#00CED1'    [   0 206 209 ]
  'medium turquoise'        '#48D1CC'    [  72 209 204 ]
  'MediumTurquoise'         '#48D1CC'    [  72 209 204 ]
  'turquoise'               '#40E0D0'    [  64 224 208 ]
  'cyan'                    '#00FFFF'    [   0 255 255 ]
  'light cyan'              '#E0FFFF'    [ 224 255 255 ]
  'LightCyan'               '#E0FFFF'    [ 224 255 255 ]
  'cadet blue'              '#5F9EA0'    [  95 158 160 ]
  'CadetBlue'               '#5F9EA0'    [  95 158 160 ]
  'medium aquamarine'       '#66CDAA'    [ 102 205 170 ]
  'MediumAquamarine'        '#66CDAA'    [ 102 205 170 ]
  'aquamarine'              '#7FFFD4'    [ 127 255 212 ]
  'dark green'              '#006400'    [   0 100   0 ]
  'DarkGreen'               '#006400'    [   0 100   0 ]
  'dark olive green'        '#556B2F'    [  85 107  47 ]
  'DarkOliveGreen'          '#556B2F'    [  85 107  47 ]
  'dark sea green'          '#8FBC8F'    [ 143 188 143 ]
  'DarkSeaGreen'            '#8FBC8F'    [ 143 188 143 ]
  'sea green'               '#2E8B57'    [  46 139  87 ]
  'SeaGreen'                '#2E8B57'    [  46 139  87 ]
  'medium sea green'        '#3CB371'    [  60 179 113 ]
  'MediumSeaGreen'          '#3CB371'    [  60 179 113 ]
  'light sea green'         '#20B2AA'    [  32 178 170 ]
  'LightSeaGreen'           '#20B2AA'    [  32 178 170 ]
  'pale green'              '#98FB98'    [ 152 251 152 ]
  'PaleGreen'               '#98FB98'    [ 152 251 152 ]
  'spring green'            '#00FF7F'    [   0 255 127 ]
  'SpringGreen'             '#00FF7F'    [   0 255 127 ]
  'lawn green'              '#7CFC00'    [ 124 252   0 ]
  'LawnGreen'               '#7CFC00'    [ 124 252   0 ]
  'green'                   '#00FF00'    [   0 255   0 ]
  'chartreuse'              '#7FFF00'    [ 127 255   0 ]
  'medium spring green'     '#00FA9A'    [   0 250 154 ]
  'MediumSpringGreen'       '#00FA9A'    [   0 250 154 ]
  'green yellow'            '#ADFF2F'    [ 173 255  47 ]
  'GreenYellow'             '#ADFF2F'    [ 173 255  47 ]
  'lime green'              '#32CD32'    [  50 205  50 ]
  'LimeGreen'               '#32CD32'    [  50 205  50 ]
  'yellow green'            '#9ACD32'    [ 154 205  50 ]
  'YellowGreen'             '#9ACD32'    [ 154 205  50 ]
  'forest green'            '#228B22'    [  34 139  34 ]
  'ForestGreen'             '#228B22'    [  34 139  34 ]
  'olive drab'              '#6B8E23'    [ 107 142  35 ]
  'OliveDrab'               '#6B8E23'    [ 107 142  35 ]
  'dark khaki'              '#BDB76B'    [ 189 183 107 ]
  'DarkKhaki'               '#BDB76B'    [ 189 183 107 ]
  'khaki'                   '#F0E68C'    [ 240 230 140 ]
  'pale goldenrod'          '#EEE8AA'    [ 238 232 170 ]
  'PaleGoldenrod'           '#EEE8AA'    [ 238 232 170 ]
  'light goldenrod yellow'  '#FAFAD2'    [ 250 250 210 ]
  'LightGoldenrodYellow'    '#FAFAD2'    [ 250 250 210 ]
  'light yellow'            '#FFFFE0'    [ 255 255 224 ]
  'LightYellow'             '#FFFFE0'    [ 255 255 224 ]
  'yellow'                  '#FFFF00'    [ 255 255   0 ]
  'gold'                    '#FFD700'    [ 255 215   0 ]
  'light goldenrod'         '#EEDD82'    [ 238 221 130 ]
  'LightGoldenrod'          '#EEDD82'    [ 238 221 130 ]
  'goldenrod'               '#DAA520'    [ 218 165  32 ]
  'dark goldenrod'          '#B8860B'    [ 184 134  11 ]
  'DarkGoldenrod'           '#B8860B'    [ 184 134  11 ]
  'rosy brown'              '#BC8F8F'    [ 188 143 143 ]
  'RosyBrown'               '#BC8F8F'    [ 188 143 143 ]
  'indian red'              '#CD5C5C'    [ 205  92  92 ]
  'IndianRed'               '#CD5C5C'    [ 205  92  92 ]
  'saddle brown'            '#8B4513'    [ 139  69  19 ]
  'SaddleBrown'             '#8B4513'    [ 139  69  19 ]
  'sienna'                  '#A0522D'    [ 160  82  45 ]
  'peru'                    '#CD853F'    [ 205 133  63 ]
  'burlywood'               '#DEB887'    [ 222 184 135 ]
  'beige'                   '#F5F5DC'    [ 245 245 220 ]
  'wheat'                   '#F5DEB3'    [ 245 222 179 ]
  'sandy brown'             '#F4A460'    [ 244 164  96 ]
  'SandyBrown'              '#F4A460'    [ 244 164  96 ]
  'tan'                     '#D2B48C'    [ 210 180 140 ]
  'chocolate'               '#D2691E'    [ 210 105  30 ]
  'firebrick'               '#B22222'    [ 178  34  34 ]
  'brown'                   '#A52A2A'    [ 165  42  42 ]
  'dark salmon'             '#E9967A'    [ 233 150 122 ]
  'DarkSalmon'              '#E9967A'    [ 233 150 122 ]
  'salmon'                  '#FA8072'    [ 250 128 114 ]
  'light salmon'            '#FFA07A'    [ 255 160 122 ]
  'LightSalmon'             '#FFA07A'    [ 255 160 122 ]
  'orange'                  '#FFA500'    [ 255 165   0 ]
  'dark orange'             '#FF8C00'    [ 255 140   0 ]
  'DarkOrange'              '#FF8C00'    [ 255 140   0 ]
  'coral'                   '#FF7F50'    [ 255 127  80 ]
  'light coral'             '#F08080'    [ 240 128 128 ]
  'LightCoral'              '#F08080'    [ 240 128 128 ]
  'tomato'                  '#FF6347'    [ 255  99  71 ]
  'orange red'              '#FF4500'    [ 255  69   0 ]
  'OrangeRed'               '#FF4500'    [ 255  69   0 ]
  'red'                     '#FF0000'    [ 255   0   0 ]
  'hot pink'                '#FF69B4'    [ 255 105 180 ]
  'HotPink'                 '#FF69B4'    [ 255 105 180 ]
  'deep pink'               '#FF1493'    [ 255  20 147 ]
  'DeepPink'                '#FF1493'    [ 255  20 147 ]
  'pink'                    '#FFC0CB'    [ 255 192 203 ]
  'light pink'              '#FFB6C1'    [ 255 182 193 ]
  'LightPink'               '#FFB6C1'    [ 255 182 193 ]
  'pale violet red'         '#DB7093'    [ 219 112 147 ]
  'PaleVioletRed'           '#DB7093'    [ 219 112 147 ]
  'maroon'                  '#B03060'    [ 176  48  96 ]
  'medium violet red'       '#C71585'    [ 199  21 133 ]
  'MediumVioletRed'         '#C71585'    [ 199  21 133 ]
  'violet red'              '#D02090'    [ 208  32 144 ]
  'VioletRed'               '#D02090'    [ 208  32 144 ]
  'magenta'                 '#FF00FF'    [ 255   0 255 ]
  'violet'                  '#EE82EE'    [ 238 130 238 ]
  'plum'                    '#DDA0DD'    [ 221 160 221 ]
  'orchid'                  '#DA70D6'    [ 218 112 214 ]
  'medium orchid'           '#BA55D3'    [ 186  85 211 ]
  'MediumOrchid'            '#BA55D3'    [ 186  85 211 ]
  'dark orchid'             '#9932CC'    [ 153  50 204 ]
  'DarkOrchid'              '#9932CC'    [ 153  50 204 ]
  'dark violet'             '#9400D3'    [ 148   0 211 ]
  'DarkViolet'              '#9400D3'    [ 148   0 211 ]
  'blue violet'             '#8A2BE2'    [ 138  43 226 ]
  'BlueViolet'              '#8A2BE2'    [ 138  43 226 ]
  'purple'                  '#A020F0'    [ 160  32 240 ]
  'medium purple'           '#9370DB'    [ 147 112 219 ]
  'MediumPurple'            '#9370DB'    [ 147 112 219 ]
  'thistle'                 '#D8BFD8'    [ 216 191 216 ]
  'snow1'                   '#FFFAFA'    [ 255 250 250 ]
  'snow2'                   '#EEE9E9'    [ 238 233 233 ]
  'snow3'                   '#CDC9C9'    [ 205 201 201 ]
  'snow4'                   '#8B8989'    [ 139 137 137 ]
  'seashell1'               '#FFF5EE'    [ 255 245 238 ]
  'seashell2'               '#EEE5DE'    [ 238 229 222 ]
  'seashell3'               '#CDC5BF'    [ 205 197 191 ]
  'seashell4'               '#8B8682'    [ 139 134 130 ]
  'AntiqueWhite1'           '#FFEFDB'    [ 255 239 219 ]
  'AntiqueWhite2'           '#EEDFCC'    [ 238 223 204 ]
  'AntiqueWhite3'           '#CDC0B0'    [ 205 192 176 ]
  'AntiqueWhite4'           '#8B8378'    [ 139 131 120 ]
  'bisque1'                 '#FFE4C4'    [ 255 228 196 ]
  'bisque2'                 '#EED5B7'    [ 238 213 183 ]
  'bisque3'                 '#CDB79E'    [ 205 183 158 ]
  'bisque4'                 '#8B7D6B'    [ 139 125 107 ]
  'PeachPuff1'              '#FFDAB9'    [ 255 218 185 ]
  'PeachPuff2'              '#EECBAD'    [ 238 203 173 ]
  'PeachPuff3'              '#CDAF95'    [ 205 175 149 ]
  'PeachPuff4'              '#8B7765'    [ 139 119 101 ]
  'NavajoWhite1'            '#FFDEAD'    [ 255 222 173 ]
  'NavajoWhite2'            '#EECFA1'    [ 238 207 161 ]
  'NavajoWhite3'            '#CDB38B'    [ 205 179 139 ]
  'NavajoWhite4'            '#8B795E'    [ 139 121  94 ]
  'LemonChiffon1'           '#FFFACD'    [ 255 250 205 ]
  'LemonChiffon2'           '#EEE9BF'    [ 238 233 191 ]
  'LemonChiffon3'           '#CDC9A5'    [ 205 201 165 ]
  'LemonChiffon4'           '#8B8970'    [ 139 137 112 ]
  'cornsilk1'               '#FFF8DC'    [ 255 248 220 ]
  'cornsilk2'               '#EEE8CD'    [ 238 232 205 ]
  'cornsilk3'               '#CDC8B1'    [ 205 200 177 ]
  'cornsilk4'               '#8B8878'    [ 139 136 120 ]
  'ivory1'                  '#FFFFF0'    [ 255 255 240 ]
  'ivory2'                  '#EEEEE0'    [ 238 238 224 ]
  'ivory3'                  '#CDCDC1'    [ 205 205 193 ]
  'ivory4'                  '#8B8B83'    [ 139 139 131 ]
  'honeydew1'               '#F0FFF0'    [ 240 255 240 ]
  'honeydew2'               '#E0EEE0'    [ 224 238 224 ]
  'honeydew3'               '#C1CDC1'    [ 193 205 193 ]
  'honeydew4'               '#838B83'    [ 131 139 131 ]
  'LavenderBlush1'          '#FFF0F5'    [ 255 240 245 ]
  'LavenderBlush2'          '#EEE0E5'    [ 238 224 229 ]
  'LavenderBlush3'          '#CDC1C5'    [ 205 193 197 ]
  'LavenderBlush4'          '#8B8386'    [ 139 131 134 ]
  'MistyRose1'              '#FFE4E1'    [ 255 228 225 ]
  'MistyRose2'              '#EED5D2'    [ 238 213 210 ]
  'MistyRose3'              '#CDB7B5'    [ 205 183 181 ]
  'MistyRose4'              '#8B7D7B'    [ 139 125 123 ]
  'azure1'                  '#F0FFFF'    [ 240 255 255 ]
  'azure2'                  '#E0EEEE'    [ 224 238 238 ]
  'azure3'                  '#C1CDCD'    [ 193 205 205 ]
  'azure4'                  '#838B8B'    [ 131 139 139 ]
  'SlateBlue1'              '#836FFF'    [ 131 111 255 ]
  'SlateBlue2'              '#7A67EE'    [ 122 103 238 ]
  'SlateBlue3'              '#6959CD'    [ 105  89 205 ]
  'SlateBlue4'              '#473C8B'    [  71  60 139 ]
  'RoyalBlue1'              '#4876FF'    [  72 118 255 ]
  'RoyalBlue2'              '#436EEE'    [  67 110 238 ]
  'RoyalBlue3'              '#3A5FCD'    [  58  95 205 ]
  'RoyalBlue4'              '#27408B'    [  39  64 139 ]
  'blue1'                   '#0000FF'    [   0   0 255 ]
  'blue2'                   '#0000EE'    [   0   0 238 ]
  'blue3'                   '#0000CD'    [   0   0 205 ]
  'blue4'                   '#00008B'    [   0   0 139 ]
  'DodgerBlue1'             '#1E90FF'    [  30 144 255 ]
  'DodgerBlue2'             '#1C86EE'    [  28 134 238 ]
  'DodgerBlue3'             '#1874CD'    [  24 116 205 ]
  'DodgerBlue4'             '#104E8B'    [  16  78 139 ]
  'SteelBlue1'              '#63B8FF'    [  99 184 255 ]
  'SteelBlue2'              '#5CACEE'    [  92 172 238 ]
  'SteelBlue3'              '#4F94CD'    [  79 148 205 ]
  'SteelBlue4'              '#36648B'    [  54 100 139 ]
  'DeepSkyBlue1'            '#00BFFF'    [   0 191 255 ]
  'DeepSkyBlue2'            '#00B2EE'    [   0 178 238 ]
  'DeepSkyBlue3'            '#009ACD'    [   0 154 205 ]
  'DeepSkyBlue4'            '#00688B'    [   0 104 139 ]
  'SkyBlue1'                '#87CEFF'    [ 135 206 255 ]
  'SkyBlue2'                '#7EC0EE'    [ 126 192 238 ]
  'SkyBlue3'                '#6CA6CD'    [ 108 166 205 ]
  'SkyBlue4'                '#4A708B'    [  74 112 139 ]
  'LightSkyBlue1'           '#B0E2FF'    [ 176 226 255 ]
  'LightSkyBlue2'           '#A4D3EE'    [ 164 211 238 ]
  'LightSkyBlue3'           '#8DB6CD'    [ 141 182 205 ]
  'LightSkyBlue4'           '#607B8B'    [  96 123 139 ]
  'SlateGray1'              '#C6E2FF'    [ 198 226 255 ]
  'SlateGray2'              '#B9D3EE'    [ 185 211 238 ]
  'SlateGray3'              '#9FB6CD'    [ 159 182 205 ]
  'SlateGray4'              '#6C7B8B'    [ 108 123 139 ]
  'LightSteelBlue1'         '#CAE1FF'    [ 202 225 255 ]
  'LightSteelBlue2'         '#BCD2EE'    [ 188 210 238 ]
  'LightSteelBlue3'         '#A2B5CD'    [ 162 181 205 ]
  'LightSteelBlue4'         '#6E7B8B'    [ 110 123 139 ]
  'LightBlue1'              '#BFEFFF'    [ 191 239 255 ]
  'LightBlue2'              '#B2DFEE'    [ 178 223 238 ]
  'LightBlue3'              '#9AC0CD'    [ 154 192 205 ]
  'LightBlue4'              '#68838B'    [ 104 131 139 ]
  'LightCyan1'              '#E0FFFF'    [ 224 255 255 ]
  'LightCyan2'              '#D1EEEE'    [ 209 238 238 ]
  'LightCyan3'              '#B4CDCD'    [ 180 205 205 ]
  'LightCyan4'              '#7A8B8B'    [ 122 139 139 ]
  'PaleTurquoise1'          '#BBFFFF'    [ 187 255 255 ]
  'PaleTurquoise2'          '#AEEEEE'    [ 174 238 238 ]
  'PaleTurquoise3'          '#96CDCD'    [ 150 205 205 ]
  'PaleTurquoise4'          '#668B8B'    [ 102 139 139 ]
  'CadetBlue1'              '#98F5FF'    [ 152 245 255 ]
  'CadetBlue2'              '#8EE5EE'    [ 142 229 238 ]
  'CadetBlue3'              '#7AC5CD'    [ 122 197 205 ]
  'CadetBlue4'              '#53868B'    [  83 134 139 ]
  'turquoise1'              '#00F5FF'    [   0 245 255 ]
  'turquoise2'              '#00E5EE'    [   0 229 238 ]
  'turquoise3'              '#00C5CD'    [   0 197 205 ]
  'turquoise4'              '#00868B'    [   0 134 139 ]
  'cyan1'                   '#00FFFF'    [   0 255 255 ]
  'cyan2'                   '#00EEEE'    [   0 238 238 ]
  'cyan3'                   '#00CDCD'    [   0 205 205 ]
  'cyan4'                   '#008B8B'    [   0 139 139 ]
  'DarkSlateGray1'          '#97FFFF'    [ 151 255 255 ]
  'DarkSlateGray2'          '#8DEEEE'    [ 141 238 238 ]
  'DarkSlateGray3'          '#79CDCD'    [ 121 205 205 ]
  'DarkSlateGray4'          '#528B8B'    [  82 139 139 ]
  'aquamarine1'             '#7FFFD4'    [ 127 255 212 ]
  'aquamarine2'             '#76EEC6'    [ 118 238 198 ]
  'aquamarine3'             '#66CDAA'    [ 102 205 170 ]
  'aquamarine4'             '#458B74'    [  69 139 116 ]
  'DarkSeaGreen1'           '#C1FFC1'    [ 193 255 193 ]
  'DarkSeaGreen2'           '#B4EEB4'    [ 180 238 180 ]
  'DarkSeaGreen3'           '#9BCD9B'    [ 155 205 155 ]
  'DarkSeaGreen4'           '#698B69'    [ 105 139 105 ]
  'SeaGreen1'               '#54FF9F'    [  84 255 159 ]
  'SeaGreen2'               '#4EEE94'    [  78 238 148 ]
  'SeaGreen3'               '#43CD80'    [  67 205 128 ]
  'SeaGreen4'               '#2E8B57'    [  46 139  87 ]
  'PaleGreen1'              '#9AFF9A'    [ 154 255 154 ]
  'PaleGreen2'              '#90EE90'    [ 144 238 144 ]
  'PaleGreen3'              '#7CCD7C'    [ 124 205 124 ]
  'PaleGreen4'              '#548B54'    [  84 139  84 ]
  'SpringGreen1'            '#00FF7F'    [   0 255 127 ]
  'SpringGreen2'            '#00EE76'    [   0 238 118 ]
  'SpringGreen3'            '#00CD66'    [   0 205 102 ]
  'SpringGreen4'            '#008B45'    [   0 139  69 ]
  'green1'                  '#00FF00'    [   0 255   0 ]
  'green2'                  '#00EE00'    [   0 238   0 ]
  'green3'                  '#00CD00'    [   0 205   0 ]
  'green4'                  '#008B00'    [   0 139   0 ]
  'chartreuse1'             '#7FFF00'    [ 127 255   0 ]
  'chartreuse2'             '#76EE00'    [ 118 238   0 ]
  'chartreuse3'             '#66CD00'    [ 102 205   0 ]
  'chartreuse4'             '#458B00'    [  69 139   0 ]
  'OliveDrab1'              '#C0FF3E'    [ 192 255  62 ]
  'OliveDrab2'              '#B3EE3A'    [ 179 238  58 ]
  'OliveDrab3'              '#9ACD32'    [ 154 205  50 ]
  'OliveDrab4'              '#698B22'    [ 105 139  34 ]
  'DarkOliveGreen1'         '#CAFF70'    [ 202 255 112 ]
  'DarkOliveGreen2'         '#BCEE68'    [ 188 238 104 ]
  'DarkOliveGreen3'         '#A2CD5A'    [ 162 205  90 ]
  'DarkOliveGreen4'         '#6E8B3D'    [ 110 139  61 ]
  'khaki1'                  '#FFF68F'    [ 255 246 143 ]
  'khaki2'                  '#EEE685'    [ 238 230 133 ]
  'khaki3'                  '#CDC673'    [ 205 198 115 ]
  'khaki4'                  '#8B864E'    [ 139 134  78 ]
  'LightGoldenrod1'         '#FFEC8B'    [ 255 236 139 ]
  'LightGoldenrod2'         '#EEDC82'    [ 238 220 130 ]
  'LightGoldenrod3'         '#CDBE70'    [ 205 190 112 ]
  'LightGoldenrod4'         '#8B814C'    [ 139 129  76 ]
  'LightYellow1'            '#FFFFE0'    [ 255 255 224 ]
  'LightYellow2'            '#EEEED1'    [ 238 238 209 ]
  'LightYellow3'            '#CDCDB4'    [ 205 205 180 ]
  'LightYellow4'            '#8B8B7A'    [ 139 139 122 ]
  'yellow1'                 '#FFFF00'    [ 255 255   0 ]
  'yellow2'                 '#EEEE00'    [ 238 238   0 ]
  'yellow3'                 '#CDCD00'    [ 205 205   0 ]
  'yellow4'                 '#8B8B00'    [ 139 139   0 ]
  'gold1'                   '#FFD700'    [ 255 215   0 ]
  'gold2'                   '#EEC900'    [ 238 201   0 ]
  'gold3'                   '#CDAD00'    [ 205 173   0 ]
  'gold4'                   '#8B7500'    [ 139 117   0 ]
  'goldenrod1'              '#FFC125'    [ 255 193  37 ]
  'goldenrod2'              '#EEB422'    [ 238 180  34 ]
  'goldenrod3'              '#CD9B1D'    [ 205 155  29 ]
  'goldenrod4'              '#8B6914'    [ 139 105  20 ]
  'DarkGoldenrod1'          '#FFB90F'    [ 255 185  15 ]
  'DarkGoldenrod2'          '#EEAD0E'    [ 238 173  14 ]
  'DarkGoldenrod3'          '#CD950C'    [ 205 149  12 ]
  'DarkGoldenrod4'          '#8B6508'    [ 139 101   8 ]
  'RosyBrown1'              '#FFC1C1'    [ 255 193 193 ]
  'RosyBrown2'              '#EEB4B4'    [ 238 180 180 ]
  'RosyBrown3'              '#CD9B9B'    [ 205 155 155 ]
  'RosyBrown4'              '#8B6969'    [ 139 105 105 ]
  'IndianRed1'              '#FF6A6A'    [ 255 106 106 ]
  'IndianRed2'              '#EE6363'    [ 238  99  99 ]
  'IndianRed3'              '#CD5555'    [ 205  85  85 ]
  'IndianRed4'              '#8B3A3A'    [ 139  58  58 ]
  'sienna1'                 '#FF8247'    [ 255 130  71 ]
  'sienna2'                 '#EE7942'    [ 238 121  66 ]
  'sienna3'                 '#CD6839'    [ 205 104  57 ]
  'sienna4'                 '#8B4726'    [ 139  71  38 ]
  'burlywood1'              '#FFD39B'    [ 255 211 155 ]
  'burlywood2'              '#EEC591'    [ 238 197 145 ]
  'burlywood3'              '#CDAA7D'    [ 205 170 125 ]
  'burlywood4'              '#8B7355'    [ 139 115  85 ]
  'wheat1'                  '#FFE7BA'    [ 255 231 186 ]
  'wheat2'                  '#EED8AE'    [ 238 216 174 ]
  'wheat3'                  '#CDBA96'    [ 205 186 150 ]
  'wheat4'                  '#8B7E66'    [ 139 126 102 ]
  'tan1'                    '#FFA54F'    [ 255 165  79 ]
  'tan2'                    '#EE9A49'    [ 238 154  73 ]
  'tan3'                    '#CD853F'    [ 205 133  63 ]
  'tan4'                    '#8B5A2B'    [ 139  90  43 ]
  'chocolate1'              '#FF7F24'    [ 255 127  36 ]
  'chocolate2'              '#EE7621'    [ 238 118  33 ]
  'chocolate3'              '#CD661D'    [ 205 102  29 ]
  'chocolate4'              '#8B4513'    [ 139  69  19 ]
  'firebrick1'              '#FF3030'    [ 255  48  48 ]
  'firebrick2'              '#EE2C2C'    [ 238  44  44 ]
  'firebrick3'              '#CD2626'    [ 205  38  38 ]
  'firebrick4'              '#8B1A1A'    [ 139  26  26 ]
  'brown1'                  '#FF4040'    [ 255  64  64 ]
  'brown2'                  '#EE3B3B'    [ 238  59  59 ]
  'brown3'                  '#CD3333'    [ 205  51  51 ]
  'brown4'                  '#8B2323'    [ 139  35  35 ]
  'salmon1'                 '#FF8C69'    [ 255 140 105 ]
  'salmon2'                 '#EE8262'    [ 238 130  98 ]
  'salmon3'                 '#CD7054'    [ 205 112  84 ]
  'salmon4'                 '#8B4C39'    [ 139  76  57 ]
  'LightSalmon1'            '#FFA07A'    [ 255 160 122 ]
  'LightSalmon2'            '#EE9572'    [ 238 149 114 ]
  'LightSalmon3'            '#CD8162'    [ 205 129  98 ]
  'LightSalmon4'            '#8B5742'    [ 139  87  66 ]
  'orange1'                 '#FFA500'    [ 255 165   0 ]
  'orange2'                 '#EE9A00'    [ 238 154   0 ]
  'orange3'                 '#CD8500'    [ 205 133   0 ]
  'orange4'                 '#8B5A00'    [ 139  90   0 ]
  'DarkOrange1'             '#FF7F00'    [ 255 127   0 ]
  'DarkOrange2'             '#EE7600'    [ 238 118   0 ]
  'DarkOrange3'             '#CD6600'    [ 205 102   0 ]
  'DarkOrange4'             '#8B4500'    [ 139  69   0 ]
  'coral1'                  '#FF7256'    [ 255 114  86 ]
  'coral2'                  '#EE6A50'    [ 238 106  80 ]
  'coral3'                  '#CD5B45'    [ 205  91  69 ]
  'coral4'                  '#8B3E2F'    [ 139  62  47 ]
  'tomato1'                 '#FF6347'    [ 255  99  71 ]
  'tomato2'                 '#EE5C42'    [ 238  92  66 ]
  'tomato3'                 '#CD4F39'    [ 205  79  57 ]
  'tomato4'                 '#8B3626'    [ 139  54  38 ]
  'OrangeRed1'              '#FF4500'    [ 255  69   0 ]
  'OrangeRed2'              '#EE4000'    [ 238  64   0 ]
  'OrangeRed3'              '#CD3700'    [ 205  55   0 ]
  'OrangeRed4'              '#8B2500'    [ 139  37   0 ]
  'red1'                    '#FF0000'    [ 255   0   0 ]
  'red2'                    '#EE0000'    [ 238   0   0 ]
  'red3'                    '#CD0000'    [ 205   0   0 ]
  'red4'                    '#8B0000'    [ 139   0   0 ]
  'DeepPink1'               '#FF1493'    [ 255  20 147 ]
  'DeepPink2'               '#EE1289'    [ 238  18 137 ]
  'DeepPink3'               '#CD1076'    [ 205  16 118 ]
  'DeepPink4'               '#8B0A50'    [ 139  10  80 ]
  'HotPink1'                '#FF6EB4'    [ 255 110 180 ]
  'HotPink2'                '#EE6AA7'    [ 238 106 167 ]
  'HotPink3'                '#CD6090'    [ 205  96 144 ]
  'HotPink4'                '#8B3A62'    [ 139  58  98 ]
  'pink1'                   '#FFB5C5'    [ 255 181 197 ]
  'pink2'                   '#EEA9B8'    [ 238 169 184 ]
  'pink3'                   '#CD919E'    [ 205 145 158 ]
  'pink4'                   '#8B636C'    [ 139  99 108 ]
  'LightPink1'              '#FFAEB9'    [ 255 174 185 ]
  'LightPink2'              '#EEA2AD'    [ 238 162 173 ]
  'LightPink3'              '#CD8C95'    [ 205 140 149 ]
  'LightPink4'              '#8B5F65'    [ 139  95 101 ]
  'PaleVioletRed1'          '#FF82AB'    [ 255 130 171 ]
  'PaleVioletRed2'          '#EE799F'    [ 238 121 159 ]
  'PaleVioletRed3'          '#CD6889'    [ 205 104 137 ]
  'PaleVioletRed4'          '#8B475D'    [ 139  71  93 ]
  'maroon1'                 '#FF34B3'    [ 255  52 179 ]
  'maroon2'                 '#EE30A7'    [ 238  48 167 ]
  'maroon3'                 '#CD2990'    [ 205  41 144 ]
  'maroon4'                 '#8B1C62'    [ 139  28  98 ]
  'VioletRed1'              '#FF3E96'    [ 255  62 150 ]
  'VioletRed2'              '#EE3A8C'    [ 238  58 140 ]
  'VioletRed3'              '#CD3278'    [ 205  50 120 ]
  'VioletRed4'              '#8B2252'    [ 139  34  82 ]
  'magenta1'                '#FF00FF'    [ 255   0 255 ]
  'magenta2'                '#EE00EE'    [ 238   0 238 ]
  'magenta3'                '#CD00CD'    [ 205   0 205 ]
  'magenta4'                '#8B008B'    [ 139   0 139 ]
  'orchid1'                 '#FF83FA'    [ 255 131 250 ]
  'orchid2'                 '#EE7AE9'    [ 238 122 233 ]
  'orchid3'                 '#CD69C9'    [ 205 105 201 ]
  'orchid4'                 '#8B4789'    [ 139  71 137 ]
  'plum1'                   '#FFBBFF'    [ 255 187 255 ]
  'plum2'                   '#EEAEEE'    [ 238 174 238 ]
  'plum3'                   '#CD96CD'    [ 205 150 205 ]
  'plum4'                   '#8B668B'    [ 139 102 139 ]
  'MediumOrchid1'           '#E066FF'    [ 224 102 255 ]
  'MediumOrchid2'           '#D15FEE'    [ 209  95 238 ]
  'MediumOrchid3'           '#B452CD'    [ 180  82 205 ]
  'MediumOrchid4'           '#7A378B'    [ 122  55 139 ]
  'DarkOrchid1'             '#BF3EFF'    [ 191  62 255 ]
  'DarkOrchid2'             '#B23AEE'    [ 178  58 238 ]
  'DarkOrchid3'             '#9A32CD'    [ 154  50 205 ]
  'DarkOrchid4'             '#68228B'    [ 104  34 139 ]
  'purple1'                 '#9B30FF'    [ 155  48 255 ]
  'purple2'                 '#912CEE'    [ 145  44 238 ]
  'purple3'                 '#7D26CD'    [ 125  38 205 ]
  'purple4'                 '#551A8B'    [  85  26 139 ]
  'MediumPurple1'           '#AB82FF'    [ 171 130 255 ]
  'MediumPurple2'           '#9F79EE'    [ 159 121 238 ]
  'MediumPurple3'           '#8968CD'    [ 137 104 205 ]
  'MediumPurple4'           '#5D478B'    [  93  71 139 ]
  'thistle1'                '#FFE1FF'    [ 255 225 255 ]
  'thistle2'                '#EED2EE'    [ 238 210 238 ]
  'thistle3'                '#CDB5CD'    [ 205 181 205 ]
  'thistle4'                '#8B7B8B'    [ 139 123 139 ]
  'gray0'                   '#000000'    [   0   0   0 ]
  'grey0'                   '#000000'    [   0   0   0 ]
  'gray1'                   '#030303'    [   3   3   3 ]
  'grey1'                   '#030303'    [   3   3   3 ]
  'gray2'                   '#050505'    [   5   5   5 ]
  'grey2'                   '#050505'    [   5   5   5 ]
  'gray3'                   '#080808'    [   8   8   8 ]
  'grey3'                   '#080808'    [   8   8   8 ]
  'gray4'                   '#0A0A0A'    [  10  10  10 ]
  'grey4'                   '#0A0A0A'    [  10  10  10 ]
  'gray5'                   '#0D0D0D'    [  13  13  13 ]
  'grey5'                   '#0D0D0D'    [  13  13  13 ]
  'gray6'                   '#0F0F0F'    [  15  15  15 ]
  'grey6'                   '#0F0F0F'    [  15  15  15 ]
  'gray7'                   '#121212'    [  18  18  18 ]
  'grey7'                   '#121212'    [  18  18  18 ]
  'gray8'                   '#141414'    [  20  20  20 ]
  'grey8'                   '#141414'    [  20  20  20 ]
  'gray9'                   '#171717'    [  23  23  23 ]
  'grey9'                   '#171717'    [  23  23  23 ]
  'gray10'                  '#1A1A1A'    [  26  26  26 ]
  'grey10'                  '#1A1A1A'    [  26  26  26 ]
  'gray11'                  '#1C1C1C'    [  28  28  28 ]
  'grey11'                  '#1C1C1C'    [  28  28  28 ]
  'gray12'                  '#1F1F1F'    [  31  31  31 ]
  'grey12'                  '#1F1F1F'    [  31  31  31 ]
  'gray13'                  '#212121'    [  33  33  33 ]
  'grey13'                  '#212121'    [  33  33  33 ]
  'gray14'                  '#242424'    [  36  36  36 ]
  'grey14'                  '#242424'    [  36  36  36 ]
  'gray15'                  '#262626'    [  38  38  38 ]
  'grey15'                  '#262626'    [  38  38  38 ]
  'gray16'                  '#292929'    [  41  41  41 ]
  'grey16'                  '#292929'    [  41  41  41 ]
  'gray17'                  '#2B2B2B'    [  43  43  43 ]
  'grey17'                  '#2B2B2B'    [  43  43  43 ]
  'gray18'                  '#2E2E2E'    [  46  46  46 ]
  'grey18'                  '#2E2E2E'    [  46  46  46 ]
  'gray19'                  '#303030'    [  48  48  48 ]
  'grey19'                  '#303030'    [  48  48  48 ]
  'gray20'                  '#333333'    [  51  51  51 ]
  'grey20'                  '#333333'    [  51  51  51 ]
  'gray21'                  '#363636'    [  54  54  54 ]
  'grey21'                  '#363636'    [  54  54  54 ]
  'gray22'                  '#383838'    [  56  56  56 ]
  'grey22'                  '#383838'    [  56  56  56 ]
  'gray23'                  '#3B3B3B'    [  59  59  59 ]
  'grey23'                  '#3B3B3B'    [  59  59  59 ]
  'gray24'                  '#3D3D3D'    [  61  61  61 ]
  'grey24'                  '#3D3D3D'    [  61  61  61 ]
  'gray25'                  '#404040'    [  64  64  64 ]
  'grey25'                  '#404040'    [  64  64  64 ]
  'gray26'                  '#424242'    [  66  66  66 ]
  'grey26'                  '#424242'    [  66  66  66 ]
  'gray27'                  '#454545'    [  69  69  69 ]
  'grey27'                  '#454545'    [  69  69  69 ]
  'gray28'                  '#474747'    [  71  71  71 ]
  'grey28'                  '#474747'    [  71  71  71 ]
  'gray29'                  '#4A4A4A'    [  74  74  74 ]
  'grey29'                  '#4A4A4A'    [  74  74  74 ]
  'gray30'                  '#4D4D4D'    [  77  77  77 ]
  'grey30'                  '#4D4D4D'    [  77  77  77 ]
  'gray31'                  '#4F4F4F'    [  79  79  79 ]
  'grey31'                  '#4F4F4F'    [  79  79  79 ]
  'gray32'                  '#525252'    [  82  82  82 ]
  'grey32'                  '#525252'    [  82  82  82 ]
  'gray33'                  '#545454'    [  84  84  84 ]
  'grey33'                  '#545454'    [  84  84  84 ]
  'gray34'                  '#575757'    [  87  87  87 ]
  'grey34'                  '#575757'    [  87  87  87 ]
  'gray35'                  '#595959'    [  89  89  89 ]
  'grey35'                  '#595959'    [  89  89  89 ]
  'gray36'                  '#5C5C5C'    [  92  92  92 ]
  'grey36'                  '#5C5C5C'    [  92  92  92 ]
  'gray37'                  '#5E5E5E'    [  94  94  94 ]
  'grey37'                  '#5E5E5E'    [  94  94  94 ]
  'gray38'                  '#616161'    [  97  97  97 ]
  'grey38'                  '#616161'    [  97  97  97 ]
  'gray39'                  '#636363'    [  99  99  99 ]
  'grey39'                  '#636363'    [  99  99  99 ]
  'gray40'                  '#666666'    [ 102 102 102 ]
  'grey40'                  '#666666'    [ 102 102 102 ]
  'gray41'                  '#696969'    [ 105 105 105 ]
  'grey41'                  '#696969'    [ 105 105 105 ]
  'gray42'                  '#6B6B6B'    [ 107 107 107 ]
  'grey42'                  '#6B6B6B'    [ 107 107 107 ]
  'gray43'                  '#6E6E6E'    [ 110 110 110 ]
  'grey43'                  '#6E6E6E'    [ 110 110 110 ]
  'gray44'                  '#707070'    [ 112 112 112 ]
  'grey44'                  '#707070'    [ 112 112 112 ]
  'gray45'                  '#737373'    [ 115 115 115 ]
  'grey45'                  '#737373'    [ 115 115 115 ]
  'gray46'                  '#757575'    [ 117 117 117 ]
  'grey46'                  '#757575'    [ 117 117 117 ]
  'gray47'                  '#787878'    [ 120 120 120 ]
  'grey47'                  '#787878'    [ 120 120 120 ]
  'gray48'                  '#7A7A7A'    [ 122 122 122 ]
  'grey48'                  '#7A7A7A'    [ 122 122 122 ]
  'gray49'                  '#7D7D7D'    [ 125 125 125 ]
  'grey49'                  '#7D7D7D'    [ 125 125 125 ]
  'gray50'                  '#7F7F7F'    [ 127 127 127 ]
  'grey50'                  '#7F7F7F'    [ 127 127 127 ]
  'gray51'                  '#828282'    [ 130 130 130 ]
  'grey51'                  '#828282'    [ 130 130 130 ]
  'gray52'                  '#858585'    [ 133 133 133 ]
  'grey52'                  '#858585'    [ 133 133 133 ]
  'gray53'                  '#878787'    [ 135 135 135 ]
  'grey53'                  '#878787'    [ 135 135 135 ]
  'gray54'                  '#8A8A8A'    [ 138 138 138 ]
  'grey54'                  '#8A8A8A'    [ 138 138 138 ]
  'gray55'                  '#8C8C8C'    [ 140 140 140 ]
  'grey55'                  '#8C8C8C'    [ 140 140 140 ]
  'gray56'                  '#8F8F8F'    [ 143 143 143 ]
  'grey56'                  '#8F8F8F'    [ 143 143 143 ]
  'gray57'                  '#919191'    [ 145 145 145 ]
  'grey57'                  '#919191'    [ 145 145 145 ]
  'gray58'                  '#949494'    [ 148 148 148 ]
  'grey58'                  '#949494'    [ 148 148 148 ]
  'gray59'                  '#969696'    [ 150 150 150 ]
  'grey59'                  '#969696'    [ 150 150 150 ]
  'gray60'                  '#999999'    [ 153 153 153 ]
  'grey60'                  '#999999'    [ 153 153 153 ]
  'gray61'                  '#9C9C9C'    [ 156 156 156 ]
  'grey61'                  '#9C9C9C'    [ 156 156 156 ]
  'gray62'                  '#9E9E9E'    [ 158 158 158 ]
  'grey62'                  '#9E9E9E'    [ 158 158 158 ]
  'gray63'                  '#A1A1A1'    [ 161 161 161 ]
  'grey63'                  '#A1A1A1'    [ 161 161 161 ]
  'gray64'                  '#A3A3A3'    [ 163 163 163 ]
  'grey64'                  '#A3A3A3'    [ 163 163 163 ]
  'gray65'                  '#A6A6A6'    [ 166 166 166 ]
  'grey65'                  '#A6A6A6'    [ 166 166 166 ]
  'gray66'                  '#A8A8A8'    [ 168 168 168 ]
  'grey66'                  '#A8A8A8'    [ 168 168 168 ]
  'gray67'                  '#ABABAB'    [ 171 171 171 ]
  'grey67'                  '#ABABAB'    [ 171 171 171 ]
  'gray68'                  '#ADADAD'    [ 173 173 173 ]
  'grey68'                  '#ADADAD'    [ 173 173 173 ]
  'gray69'                  '#B0B0B0'    [ 176 176 176 ]
  'grey69'                  '#B0B0B0'    [ 176 176 176 ]
  'gray70'                  '#B3B3B3'    [ 179 179 179 ]
  'grey70'                  '#B3B3B3'    [ 179 179 179 ]
  'gray71'                  '#B5B5B5'    [ 181 181 181 ]
  'grey71'                  '#B5B5B5'    [ 181 181 181 ]
  'gray72'                  '#B8B8B8'    [ 184 184 184 ]
  'grey72'                  '#B8B8B8'    [ 184 184 184 ]
  'gray73'                  '#BABABA'    [ 186 186 186 ]
  'grey73'                  '#BABABA'    [ 186 186 186 ]
  'gray74'                  '#BDBDBD'    [ 189 189 189 ]
  'grey74'                  '#BDBDBD'    [ 189 189 189 ]
  'gray75'                  '#BFBFBF'    [ 191 191 191 ]
  'grey75'                  '#BFBFBF'    [ 191 191 191 ]
  'gray76'                  '#C2C2C2'    [ 194 194 194 ]
  'grey76'                  '#C2C2C2'    [ 194 194 194 ]
  'gray77'                  '#C4C4C4'    [ 196 196 196 ]
  'grey77'                  '#C4C4C4'    [ 196 196 196 ]
  'gray78'                  '#C7C7C7'    [ 199 199 199 ]
  'grey78'                  '#C7C7C7'    [ 199 199 199 ]
  'gray79'                  '#C9C9C9'    [ 201 201 201 ]
  'grey79'                  '#C9C9C9'    [ 201 201 201 ]
  'gray80'                  '#CCCCCC'    [ 204 204 204 ]
  'grey80'                  '#CCCCCC'    [ 204 204 204 ]
  'gray81'                  '#CFCFCF'    [ 207 207 207 ]
  'grey81'                  '#CFCFCF'    [ 207 207 207 ]
  'gray82'                  '#D1D1D1'    [ 209 209 209 ]
  'grey82'                  '#D1D1D1'    [ 209 209 209 ]
  'gray83'                  '#D4D4D4'    [ 212 212 212 ]
  'grey83'                  '#D4D4D4'    [ 212 212 212 ]
  'gray84'                  '#D6D6D6'    [ 214 214 214 ]
  'grey84'                  '#D6D6D6'    [ 214 214 214 ]
  'gray85'                  '#D9D9D9'    [ 217 217 217 ]
  'grey85'                  '#D9D9D9'    [ 217 217 217 ]
  'gray86'                  '#DBDBDB'    [ 219 219 219 ]
  'grey86'                  '#DBDBDB'    [ 219 219 219 ]
  'gray87'                  '#DEDEDE'    [ 222 222 222 ]
  'grey87'                  '#DEDEDE'    [ 222 222 222 ]
  'gray88'                  '#E0E0E0'    [ 224 224 224 ]
  'grey88'                  '#E0E0E0'    [ 224 224 224 ]
  'gray89'                  '#E3E3E3'    [ 227 227 227 ]
  'grey89'                  '#E3E3E3'    [ 227 227 227 ]
  'gray90'                  '#E5E5E5'    [ 229 229 229 ]
  'grey90'                  '#E5E5E5'    [ 229 229 229 ]
  'gray91'                  '#E8E8E8'    [ 232 232 232 ]
  'grey91'                  '#E8E8E8'    [ 232 232 232 ]
  'gray92'                  '#EBEBEB'    [ 235 235 235 ]
  'grey92'                  '#EBEBEB'    [ 235 235 235 ]
  'gray93'                  '#EDEDED'    [ 237 237 237 ]
  'grey93'                  '#EDEDED'    [ 237 237 237 ]
  'gray94'                  '#F0F0F0'    [ 240 240 240 ]
  'grey94'                  '#F0F0F0'    [ 240 240 240 ]
  'gray95'                  '#F2F2F2'    [ 242 242 242 ]
  'grey95'                  '#F2F2F2'    [ 242 242 242 ]
  'gray96'                  '#F5F5F5'    [ 245 245 245 ]
  'grey96'                  '#F5F5F5'    [ 245 245 245 ]
  'gray97'                  '#F7F7F7'    [ 247 247 247 ]
  'grey97'                  '#F7F7F7'    [ 247 247 247 ]
  'gray98'                  '#FAFAFA'    [ 250 250 250 ]
  'grey98'                  '#FAFAFA'    [ 250 250 250 ]
  'gray99'                  '#FCFCFC'    [ 252 252 252 ]
  'grey99'                  '#FCFCFC'    [ 252 252 252 ]
  'gray100'                 '#FFFFFF'    [ 255 255 255 ]
  'grey100'                 '#FFFFFF'    [ 255 255 255 ]
  'dark grey'               '#A9A9A9'    [ 169 169 169 ]
  'DarkGrey'                '#A9A9A9'    [ 169 169 169 ]
  'dark gray'               '#A9A9A9'    [ 169 169 169 ]
  'DarkGray'                '#A9A9A9'    [ 169 169 169 ]
  'dark blue'               '#00008B'    [   0   0 139 ]
  'DarkBlue'                '#00008B'    [   0   0 139 ]
  'dark cyan'               '#008B8B'    [   0 139 139 ]
  'DarkCyan'                '#008B8B'    [   0 139 139 ]
  'dark magenta'            '#8B008B'    [ 139   0 139 ]
  'DarkMagenta'             '#8B008B'    [ 139   0 139 ]
  'dark red'                '#8B0000'    [ 139   0   0 ]
  'DarkRed'                 '#8B0000'    [ 139   0   0 ]
  'light green'             '#90EE90'    [ 144 238 144 ]
  'LightGreen'              '#90EE90'    [ 144 238 144 ]

  };

