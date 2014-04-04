function [msg,file,fig,scl] = saveimg(varargin);

% SAVEIMG  Saves Image of Figure
%
% SAVEIMG uses PRINT:
%
%    Msg = SAVEIMG( Fig , File , '.<Format>' , Scale , '-<Option>' , i*Margin )
%
% SAVEIMG uses GETFRAME and IMWRITE:
%
%    Msg = SAVEIMG( Fig , File , '*<Format>' )
%
%---------------------------------------------------------------------
%
%  Fig     Handle of Figure to print, default: GCF
%
%  File    Name of File to save,
%            in case of missing Input or '.' the Name of Figure will used, 
%            in case of empty Input ('') UIPutFile will open.           
%           A plain Name (without a Directory and Extension)
%             will expanded with Format.
%           A Wildcard-Extension ".*" will replaced by Format
%
%          
%  Format  ImageFormat, default: '.png', use any of:
%
%           dlg            Opens PrintDialog, using PRINTDLG
%
%  --- .<Format> PRINT --------------------------------------------
%
%           bmp            Windows Bitmap (BMP)
%           png            Portable Network Graphics (PNG)
%           ppm            Portable Pixmap (plain format)
%           tif            Tagged Image File Format (TIFF)
%           tiff           Tagged Image File Format (TIFF)
%           jpg<Quality>   Joint Photographic Experts Group (JPEG)
%           jpeg<Quality>  Joint Photographic Experts Group (JPEG)
%                           default: Quality = 100
% 
%           ps             PostScript
%           psc            PostScript color
%           ps2            PostScript Level 2
%           psc2           PostScript Level 2 color
%
%           eps            Encapsulated PostScript
%           epsc           Encapsulated PostScript color
%           eps2           Encapsulated PostScript Level 2
%           epsc2          Encapsulated PostScript Level 2 color
%
%          meta            Send Image to Windows clipboard in Metafile format
%          bitmap          Send Image to Windows clipboard in bitmap format
%
%  --- *<Format> GETFRAME + IMWRITE -------------------------------
%
%           bmp            Windows Bitmap (BMP)
%           png            Portable Network Graphics (PNG)
%           tif            Tagged Image File Format (TIFF)
%           tiff           Tagged Image File Format (TIFF)
%           jpg<Quality>   Joint Photographic Experts Group (JPEG)
%           jpeg<Quality>  Joint Photographic Experts Group (JPEG)
%                           default: Quality = 100
%           hdf            Hierarchical Data Format (HDF)
%           pcx            Windows Paintbrush (PCX)
%           xwd            X Window Dump (XWD)
%
%           ppm            Portable Pixmap
%           xpm            X Window Pixmap
%
%
%  Scale   ScaleFactor for Image to print, optional default: 1
%            valid for using PRINT ('.<Format>') only!
%
%  Option  printing Options for using PRINT ('.<Format>')
%
%  Margin  PaperMargin [inches] to check PaperOrientation, 
%            used for PrintDialog and PostScriptFormats.
%          Give this Input as single imaginary Number, larger ZERO.
%          The default Value for Margin is 0.5 [inch].
%
%---------------------------------------------------------------------
%
%  SAVEIMG( ... , AxesHandle , NewFontSize , NoTag , ... )
%
% Scales all absolute FontSizes of Axes- or TextObjects in Figure 
%  to new NewFontSize [points], relative to FontSize of AxesHandle. 
%
% Objects with normalized FontUnits or which Property "Tag" match
%  any of the Strings of the !CellStringArray! NoTag will NOT scaled!
%  The WildCard "*" can be used in the Strings of NoTag.
% 
% The Inputs NewFontSize and NoTag must follow the Input AxesHandle!
%
% defaults: NewFontSize = 11; NoTag = {}
%
%---------------------------------------------------------------------
%
% [Msg,GUI_Struct] = SAVEIMG( '@GUI' , CallBackFcn )
%
%   Creates a GUI_Structure to create UIMenu's using MKMENU,
%    which UserData contains the FormatDefinitions to use SAVEIMG.
%
%   The default CallBackFcn is SAVEIMG.
%
%
% [Msg,GUI_Struct] = SAVEIMG( '@GUI' , ... , Label , FieldName , Scale ... )
% 
%   returns the GUI_Structure as SubMenus of the Menu "Label".
%    Label MUST contain the Character '&' !!!
%
%   The FieldName of the SubMenu has to start with '#' !!!
%   The default FieldName is: UPPER(CallBackFcn)
%
%   Scale gives the Values for the <Scaled>-Menu.
%   In case of empty Input ([]) no <Scaled>-Menu is shown.
%
%
% [Msg,GUI_Struct,HandleStruct] = SAVEIMG( '@GUI' , ... , ParentHandle , ... )
% 
%   Creates the UIMenu's under the ParentHandle, which can be a Handle
%    of a Figure, UIContextMenu or UIMenu.
%
%   The HandleStructure of the UIMenu's is returned.
%
% Example: figure, plot(1:10)
%          [m,c,h0] = saveimg('@GUI',gcf,'&SaveAs');
%          [m,c,h1] = saveimg('@GUI',h0.SAVEIMG,'New&1',[1 2 3 4]); % 100% .. 400%
%          [m,c,h2] = saveimg('@GUI',h0.SAVEIMG,'New&2',[]);        % No Scaled
%
%---------------------------------------------------------------------
% To check the Fit: Points --> Pixels, as for LineWidth's, use:
%
% [Msg,PointFit,PixelFit] = SAVEIMG( '@SCL' , Scale , PointInt , PointMax )
%
%   PointInt  PointIntervall to check
%   PointMax  Maximum Value for Points
%
%   PointWidth to check: ( PointInt : PointInt : PointMax )
%
%   PointFit: [ PointWidth  PixelWidth ]
%   PixelFit: [ PixelWidth  MeanPointWidth MinPointWidth MaxPointWidth ]
%               unique PixelWidth's
%
%   The Inputs Scale, PointInt and PointMax are optional, defaults:
%
%       Scale = 1, PointInt = 0.1; PointMax = 5;
%
% To Show the Result in a Figure use NO Outputs or a NonZero 5. Input:
%
%  SAVEIMG( '@SCL' , Scale , PointInt , PointMax )
%
%  [Msg,PtF,PxF,Fig] = SAVEIMG( '@SCL' , Scale , PointInt , PointMax , 1 );
%
%  The MeanPointWidth per Pixel is displayed at the Labels.
%
%---------------------------------------------------------------------
%
%  Msg     contains ErrorMessage, empty if all was successful
%
% NOTE: All Inputs are optional and can be used in any order.
%       
%       If the Figure contains UIControls, the Image may be NOT
%        printed correctly, use GETFRAME and IMWRITE instead.
%
%       In the Format PPM the Fonts may be not printed equal.
%      
%------------------------------------------------------------------
%
% [ Msg , File , Fig , Scale ] = SAVEIMG( ... )
%
%  returns the used parameter, usefull in case of undefined inputs.
%
%------------------------------------------------------------------
%
% see also: PRINT, IMWRITE, WRTXPM,  WRTPPM, 
%                  IMREAD, READXPM, READPPM
%           SHOWIMG, CUTIMAGE
%

Nout = nargout;

msg = '';

%**************************************************************
% Defaults, Get Inputs

fmt = '.png';
qlt =  100;     % For JPG
scl =    1;
off =  0.5;     % PaperMargin [inches]

no_file = { 'dlg' 'meta' 'bitmap' };

[msg,fig,file,fmt,opt,qlt,typ,scl,off,axe,nfs,nscl] = ...
 checkin(varargin,fmt,qlt,scl,off,no_file,Nout);

if isempty(nfs)
   nfs = 11;    % NewFontSize
end

if ~isempty(msg) & ( Nout == 0 )
    error(msg);
end

no_file = any(strcmp(fmt,no_file));
if ~( isempty(msg) & ~isempty(fig) & ( ~isempty(file) | no_file ) )
    if ~isempty(msg) & ( Nout == 0 )
        error(msg);
    end
    if ( Nout == 0 )
       clear msg
    end
    return
end

ok = ( ishandle(fig) & ( prod(size(fig)) == 1 ) );
if ok
   ok = strcmp( get(fig,'type') , 'figure' );
end

if ~ok
    if ~isempty(msg) & ( Nout == 0 )
        error(msg);
    end
    if ( Nout == 0 )
       clear msg
    end
    return
end

%**************************************************************
% Set AxesModes, Scale FontSize

tm = setmode(fig,'tick');            % { Handle OriginalTicModes }
lm = setmode(fig,'lim');             % { Handle OriginalLimModes }

if ~isempty(axe)
    hf = repfont(fig,axe,nfs,nscl);  % [ Handle OriginalFontSize ]
end

%**************************************************************

ResizeFcn = get( fig , 'ResizeFcn' );
            set( fig , 'ResizeFcn' , '' );


%-------------------------------------------------------
if strcmp(typ,'.') | strcmp(fmt,'dlg')
%-------------------------------------------------------
   msg = printimg(fig,file,fmt,opt,qlt,scl,off);
%-------------------------------------------------------
else
%-------------------------------------------------------
   msg = writeimg(fig,file,fmt,qlt);
%-------------------------------------------------------
end
%-------------------------------------------------------

set( fig , 'ResizeFcn' , ResizeFcn );

%**************************************************************
% Reset AxesModes, Rescale FontSize

setmode(lm,'lim');
setmode(tm,'tick');

if ~isempty(axe)
    repfont(fig,hf);
end

%**************************************************************

if ~isempty(msg) & ( Nout == 0 )
    error(msg);
end

if ( Nout == 0 )
   clear msg
end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function msg = printimg(fig,file,fmt,opt,qlt,scl,papoff);

% PRINTIMG  Saves Figure, using PRINT

msg = '';

%**************************************************************
% Check Format

is_ps  = ( strcmp(fmt(1:min(2,end)),'ps') | ...
           strcmp(fmt(1:min(3,end)),'eps')      );  % PostScript

is_dlg =   strcmp(fmt,'dlg');                       % PrintDialog

is_bmp =   strcmp(fmt,'bitmap');

is_clb = ( strcmp(fmt,'meta') | is_bmp );           % ClipBoard

is_ppm = strcmp(fmt,'ppm');                         % Portable PixMap

is_vis = ~( is_dlg | is_bmp );   % Set FigureVisible OFF

if is_vis
   vis = get(fig,'visible');
   set(fig,'visible','off');
end

%--------------------------------------------------------------
% Get FigureProperties

figuni = get(fig,'units');
papuni = get(fig,'paperunits');
pappos = get(fig,'paperposition');
paport = get(fig,'paperorientation');

set(fig,'units','pixels');

figpos = get(fig,'position');

%--------------------------------------------------------------
% Check FigurePosition for Extension if ClipBoard

if is_clb

   uni = get(0,'units');
         set(0,'units','pixels');
   ssi = get(0,'screensize');
         set(0,'units',uni);

   ssi = ssi([3 4]);

   is_win = strcmp(upper(computer),'PCWIN');

   % [ Left  Bottom  Right Top ]
   figoff = [ is_bmp*[ 10  30  10 ]  60-20*is_win ];  %  [pixels]
   
   figsiz = ssi - figoff([1 2]) - figoff([3 4]); % Max. Extension
    
      scl = scl * min( min(figsiz./(figpos([3 4])*scl)) , 1 );

   sclpos = scl * figpos;

      off = sclpos([1 2]) + sclpos([3 4]) - ( ssi - figoff([3 4]));

   sclpos([1 2]) = sclpos([1 2]) - max(0,off);

   sclpos([1 2]) = max( sclpos([1 2]) , figoff([1 2]) );
   sclpos([1 2]) = max( sclpos([1 2]) , 1 );

end

%--------------------------------------------------------------
% PixelsPerInch

spi = get(0,'screenpixelsperinch');

pp0 = [ 150  spi  72 ];

ppm = 1 + 1 * ( is_dlg | is_clb | is_ps ) + 2 * is_ppm;

pp0 = pp0(ppm);

ppi = pp0 / scl;

is_scl = ( abs( ppi - spi ) > 0.1 );

if is_scl
   cnf = scale_prop(fig,ppi,spi,scl);
end

%--------------------------------------------------------------
% Set PaperProperties

set( fig , 'paperunits' , 'inches' , ...
     'paperorientation' , 'portrait' );

papsiz = get(fig,'papersize');

newpos = figpos / ppi;   % New PaperPosition

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Reduce Height by 1 Pixel !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
newpos(4) = newpos(4) - 1/pp0 * ( ppm == 1 );
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%------------------------------------------------
% Check Orientation if PrintDialog or PostScript

if is_dlg | is_ps

    psc = ( papsiz - 2*papoff ) ./ newpos([3 4]);
    psc = min( psc , 1 );

   flip = ( psc(1) < psc(2) );

    if flip
       set(fig,'paperorientation','landscape');
       papsiz = papsiz([2 1]);
    end

%------------------------------------------------
% Set scaled PaperPositon if ClipBoard

elseif is_clb

    set(fig,'position',sclpos);

end

%------------------------------------------------
% Center Image on Paper

newpos([1 2]) = ( papsiz - newpos([3 4]) ) / 2;

newpos([1 2]) = max( newpos([1 2]) , papoff );

set( fig , 'paperposition' , newpos );

%------------------------------------------------
% Set WindowStyle if ClipBoard-Bitmap

if is_bmp

   figure(fig);

   if is_scl
      set_ws = isempty(cnf.uicontrol.hdl);
   else
      shh = get(0,'showhiddenhandles');
            set(0,'showhiddenhandles','on');
      set_ws = isempty(findobj(fig,'type','uicontrol','visible','on'));
           set(0,'showhiddenhandles','on');
   end

   if set_ws
      ws = get(fig,'windowstyle');
           set(fig,'windowstyle','modal');
   end

   drawnow

end

%**************************************************************
% Print Figure

msg = '';

%-------------------------------------------------------
if is_dlg
%-------------------------------------------------------

   try
      f = printdlg(fig);
      if ~strcmp( upper(computer) , 'PCWIN' )
         waitfor(f);
      end
   catch
      msg = sprintf('Error call PRINTDLG.\n%s',lasterr);
   end

%-------------------------------------------------------
else
%-------------------------------------------------------

   form = sprintf( '-f%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 );
   ff   = sprintf(form,fig);
   dev  = sprintf('-d%s',fmt);

   if strcmp( fmt(1:min(2,end)) , 'jp' ) & ~isempty(qlt)
      dev = sprintf('%s%.0f',dev,qlt);
   end

   if isempty(file)
      file = {};
   else
      file = {file};
   end
   
   try
      print(ff,dev,opt{:},file{:});
	   catch
      str = '';
      if ~( isempty(opt) & isempty(file) )
          str = sprintf(',''%s''',opt{:},file{:});
      end
      msg = sprintf( 'Error call PRINT(''%s'',''%s''''%s'').\n%s' , ...
                      ff , dev , str , lasterr );
   end

%-------------------------------------------------------
end

%**************************************************************
% Reset Settings

if is_bmp
   if set_ws
      set(fig,'windowstyle',ws);
   end
end

if is_clb
   set(fig,'position',figpos);
   drawnow
end

set( fig ,      'units'    , figuni , ...
           'paperunits'    , papuni , ...
        'paperorientation' , paport , ...
           'paperposition' , pappos       );


if is_scl
   scale_prop(cnf);
end

if is_vis
   set(fig,'visible',vis);
end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function cnf = scale_prop(cnf,ppi,spi,scl);

% SCALEPROP  Scales Absolute Properties of Handles

%------------------------------------------
% Reset

if nargin == 1

   fld = fieldnames(cnf);

   for ff = fld(:)'

       c = getfield(cnf,ff{1});

       h = c.hdl;
       p = c.prp;
       s = c.val;
         
       m = prod(size(h));

       for jj = size(p,2)*(m>0) : -1 : 1
           v = getfield(s,p{jj});
           for kk = 1 : m
               set(h(kk),p{jj},v{kk})
           end
       end

   end

   return

end

%------------------------------------------
% Scale

fig = cnf;

shh = get(0,'showhiddenhandles');
      set(0,'showhiddenhandles','on');

%--------------------------------------------------------------
% Initialisation

%%%   { 'uicontrol'  { 'units' 'fontunits' 'fontsize' } 

ini = { 'uicontrol'  { 'units'     'fontunits' } 
        'axes'       { 'units'     'fontunits' 'fontsize'  'linewidth' }
        'text'       { 'position'      'units' 'fontunits' 'fontsize' }
        'line'       { 'linewidth' 'markersize' }
        'patch'      { 'linewidth' 'markersize' }
        'surface'    { 'linewidth' 'markersize' }
        'rectangle'  { 'linewidth' }                };

scl_prp = { 'linewidth' 'markersize' 'fontsize' };  % Properties to Scale

n = size(ini,1);

def = struct( 'hdl' , {[]} , ...
              'prp' , {[]} , ...
              'val' , {[]}       );

cnf      = permute(ini(:,[1 1]),[2 1]);
cnf(2,:) = {def};

wrn = { 'UIControl found, Image may be not printed correctly.'
        '          Use GETFRAME and IMWRITE instead.'
        '          FontSize not scaled'   };

is_uic = strcmp(ini(:,1),'uicontrol');

fu = { 'points'  'normalized' };

%--------------------------------------------------------------
% Scale Properties

for ii = 1 : n

    h = findobj(fig,'type',ini{ii,1},'visible','on');

    s = ini{ii,2}([1 1],:);

    m = prod(size(h));

    no_scl = ( ( scl == 1 ) & is_uic(ii) & any(strcmp(s(1,:),'fontsize')) );

    if 0 % !!! ~isempty(h) & is_uic(ii)
       warning(sprintf('%s\n',wrn{1:(end-(~no_scl))}));
    end

    for jj = 1 : size(s,2)*(m>0)
        p = s{1,jj};
        v = get(h,p);
        if m == 1
           v = {v};
        end
        switch p
           case 'units'
               if strcmp(ini{ii,1},'text')
                  textunit(h,'normalized');
               else
                  set(h,'units','normalized');
               end
           case 'fontunits'
               set(h,'fontunits',fu{1+is_uic(ii)});
           case  scl_prp
               if ~no_scl
                   for kk = 1 : m
                       val = v{kk} * spi/ppi;
                       if strcmp(p,'fontsize')
                          val = 0.5 * round( val / 0.5 ); % !!!
                          val = max(val,1);
                       end
                       set(h(kk),p,val);
                   end
               end
        end
        s{2,jj} = {v};
    end

    cnf{2,ii}.hdl = h;
    cnf{2,ii}.prp = ini{ii,2};
    cnf{2,ii}.val = struct(s{:});

end

cnf = struct(cnf{:});

set(0,'showhiddenhandles',shh);

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function msg = writeimg(fig,file,fmt,qlt);

% WRITEIMG  Saves Figure, using GETFRAME and IMWRITE

msg = '';

%---------------------------------------------------
% Check for UIControls

shh = get(0,'showhiddenhandles');
      set(0,'showhiddenhandles','on');

set_ws = isempty(findobj(fig,'type','uicontrol','visible','on'));

      set(0,'showhiddenhandles','on');

%---------------------------------------------------
% Capture Figure

figure(fig);

if set_ws
   ws = get(fig,'windowstyle');
        set(fig,'windowstyle','modal');
end

drawnow

 c = getframe(fig);

if set_ws
   set(fig,'windowstyle',ws);
end

%---------------------------------------------------
% Function and Options

cmap = {};
if ~isempty( c.colormap )
    cmap = { c.colormap };
end

%---------------------------------------------
if any(strcmp(fmt,{'xpm' 'ppm'}));
%---------------------------------------------

   fcn = [ 'wrt' fmt ];

   try
      msg = feval( fcn , c.cdata , cmap{:} , file );
   catch
      msg = lasterr;
   end

   if ~isempty(msg)
       msg = sprintf( 'Error call %s( CData , ''%s'' ).\n%s' , ...
                       upper(fcn) , file , msg );
   end

%---------------------------------------------
else
%---------------------------------------------

   opt = {};
   if strcmp( fmt(1:min(2,end)) , 'jp' ) & ~isempty(qlt)
      opt = { 'Quality'  qlt };
   end

   try
      imwrite( c.cdata , cmap{:} , file , fmt , opt{:} );
   catch
      msg = sprintf( 'Error call IMWRITE( CData , ''%s'' , ''%s'' ).\n%s' , ...
                      file , fmt , lasterr );
   end

%---------------------------------------------
end
%---------------------------------------------

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ini = setmode(ini,prp)

% SETMODE  Set all AxesModes to manual


cc = 'xyz';
s2 = size(cc,2);

prp = {prp};
prp = prp(ones(1,s2));

for ii = 1 : s2
     prp{ii} = cat(2,cc(ii),prp{ii},'mode');
end

ind = (1:s2) + 1;

%------------------------------------------
% Reset

if iscell(ini)

   prp = prp([1 1],:);

   for ii = 1 : size(ini,1)

       prp(2,:) = ini(ii,ind);

       set( ini{ii,1} , prp{:} );

   end

   return

end

%------------------------------------------
% Get AxesHandles and Modes

shh = get(0,'showhiddenhandles');
      set(0,'showhiddenhandles','on');

axe = findobj(ini,'type','axes');

      set(0,'showhiddenhandles',shh);

n = prod(size(axe));

ini = cell(n,s2+1);

if n == 0
   return
end

for ii = 1 : n

    ini{ii,1} = axe(ii);

    ini(ii,ind) = get( axe(ii) , prp );

end

%------------------------------------------
% Check for Mode "auto"

ok = strcmp(ini(:,ind),'auto');
ok = ~( sum(ok,2) == 0 );

if ~any(ok)
    ini = cell(0,4);
    return
end

if ~all(ok)
    ok  = find(ok);
    ini = ini(ok,:);
    axe = axe(ok);
end

%------------------------------------------
% Set Modes to "manual"

prp      = prp([1 1],:);
prp(2,:) = {'manual'};

set( axe , prp{:} );

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function hf = repfont(fig,axe,fs,no_scl)

% REPFONT  Scale all absolute Axes- and TextFontSizes of Figure


%------------------------------------------
% Reset

if nargin < 3
  
   hf = axe;

   if ~isempty(hf)
       for ii = 1 : size(hf,1)
           set( hf(ii,1) , 'fontsize' , hf(ii,2) );
       end
   end

   return

end

%------------------------------------------
% Get Axes- and TextHandles

shh = get(0,'showhiddenhandles');
      set(0,'showhiddenhandles','on');

hf = cat( 1 , findobj(fig,'type','axes') , ...
              findobj(fig,'type','text')       );

      set(0,'showhiddenhandles',shh);

if isempty(hf)
   
   hf = zeros(0,2);
 
   return

end

%------------------------------------------
% ScaleValue

fu = get(axe,'fontunits');
     set(axe,'fontunits','points');

sc = fs / get(axe,'fontsize');

     set(axe,'fontunits',fu);

%------------------------------------------
% Scale

n = size(hf,1);

hf = hf(:,[1 1]);  % [ Handle  FontSize ]
              
ok = zeros(n,1);   % Ok for NOT 'normalized'

for ii = 1 : n

   fu = get(hf(ii,1),'fontunits');

   ok(ii) = ~strcmp( fu , 'normalized' );

   if ~isempty(no_scl)
       ok(ii) = ( ok(ii) & ~strwcmp( get(hf(ii,1),'tag') , no_scl ) );
   end

   if ok(ii)

     % Original Value
      hf(ii,2) = get(hf(ii,1),'fontsize'); % Original Value

      set( hf(ii,1) , 'fontunits' , 'points' );
      set( hf(ii,1) , 'fontsize'  , get(hf(ii,1),'fontsize')*sc );
      set( hf(ii,1) , 'fontunits' , fu );

   end  

end

hf = hf(find(ok),:);

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,fig,file,fmt,opt,qlt,typ,scl,off,axe,nfs,nscl] = ...
          checkin(vin,dfmt,dqlt,dscl,doff,no_file,Nout);

msg  = '';

fig  = [];
file = -1;
fmt  = -1;
opt  = cell(1,0);
qlt  = [];
typ  = '';
scl  = [];
off  = [];

axe  = [];
nfs  = [];
nscl = {};

nin = prod(size(vin));

ind = ( 1 : nin );

%******************************************************
% Check for GUI or CallBack

if ~( nin == 0 )
    if chkhndl(vin{1},'uimenu')
       msg = menu_clb(vin{1});
       return
    elseif isequal(vin{1},'@GUI')
       [msg,file,fig] = make_gui(vin(2:end));
       return
    elseif isequal(vin{1},'@SCL')
       try
          [msg,file,fig,scl] = pointscale(Nout,vin{2:end});
       catch
           msg = sprintf('Error call POINTSCALE.\n%s',lasterr);
       end
       return
    end
end

%******************************************************
% Check for Figure

[fig,vin,nr] = hndlvarg(vin,[],'figure');

if isempty(fig)
   return
end

if ~( nr == 0 )
    ind(nr) = [];
    nin     = nin - 1;
end

%******************************************************
% Check for Axes, NewFontSize and NoScale

[axe,vin,nr] = hndlvarg(vin,[],'axes');

if ( nr == 0 )

    axe = [];

else

    for ii = nr+1 : min(nin,nr+2)

        v = vin{ii-1};

        ok = ( chkcstr(v,1) & isempty(nscl) );
        if ok
           nscl = v;
        else
           ok = ( isnumeric(v) & ( prod(size(v)) == 1 ) & isempty(nfs) );
           if ok
              ok = ( isfinite(v) & ( v > 0 ) );
           end
           if ok
              nfs = v;
           end
        end

        if ~ok
            break
        end

        nr = cat(2,nr,ii);

    end

    s2 = size(nr,2);
    if s2 > 1
       vin(nr(2:s2)-1) = [];
    end

    ind(nr) = [];
    nin     = nin - s2;

end

%******************************************************
% Check other

msg = cell(0,1);

for ii = 1 : nin
    
    v = vin{ii};

    s = size(v);
    p = prod(s);

    is_str = chkstr(v);
    is_num = isnumeric(v);

    is_fmt = 0;
    is_opt = 0;

    if is_str & ( s(2) >= 2 )
       is_fmt = (  any( v(1) == '.*' ) & isequal(fmt,-1) & ...
                  ~any( strcmp( v([1 2]) , { '..'  './'  '.\' } ) ) );
       is_opt =       ( v(1) == '-'  );
    end

    ok = ( is_str & ( isequal(file,-1) | is_fmt | is_opt ) );
    if ok
       if     is_opt
          opt = cat( 2 , opt , {v} );
       elseif is_fmt
          typ = v(1);
          fmt = v(2:s(2));
       else
          file = v;
       end
    elseif ( is_num & ( p == 1 ) )
       if ( real(v) == 0 )
          v = imag(v);
          ok = ( isfinite(v) & ( v > 0 ) & isempty(off) );
          if ok
             off = v;
          end
       else
          v = real(v);
          ok = ( isfinite(v) & ( v > 0 ) & isempty(scl) );
          if ok
             scl = v;
          end
       end
    end

    if ~ok
        str = sprintf('Invalid or unrecognized %.0f. Input.',ind(ii));
        msg = cat( 1 , msg , {str} );
    end

end

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    return
end

msg = '';

%**************************************************************
% Defaults

if isempty(scl)
   scl = dscl;
end

if isempty(off)
   off = doff;
end

if isequal(fmt,-1)
   typ = dfmt(1);
   fmt = dfmt(2:end);
end

fmt = lower(fmt);

if ~any(strcmp(fmt,no_file))
    [msg,file,fmt,qlt] = check(fig,file,fmt,dqlt);
end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,file,dev,opt] = check(fig,file,fmt,def)

% Check for FileName

msg = '';
opt = [];
dev = fmt;

%**************************************************************
% Check Format for JPEG | TIFF, get Quality

chk = { 'jpeg'  'jpg'  def
        'tiff'  'tif'  []  };


sf = size(fmt,2);

for ii = 1 : size(chk,1)

    for jj = [ 1  2 ]

        sc = size(chk{ii,jj},2);

        ok = ( sf >= sc );
        if ok
           ok = strcmp(fmt(1:sc),chk{ii,jj});
           if ok
              dev = chk{ii,1};
              if ~isempty(chk{ii,3})
                  ok = ( sf > sc );
                  if ok
                     opt = sprintf('[%s]',fmt(sc+1:sf));
                     opt = eval(opt,'[]');
                     ok  = ( isnumeric(opt) & ( prod(size(opt)) == 1 ) );
                     if ok
                        ok = ( ( mod(opt,1) == 0 ) & ( 0 < opt ) & ( opt <= 100 ) );
                     end
                  end
                  if ~ok
                      opt = chk{ii,3};
                  else
                      fmt = fmt(1:sc);
                      break
                  end
              end
           end
        end

    end

    if ok
       break
    end

end

ok = strcmp(fmt,chk(:,1));
if any(ok)
   fmt = chk{find(ok),2};
end

ext = fmt;
if     strcmp(ext(1:min(2,end)),'ps')
       ext = 'ps';
elseif strcmp(ext(1:min(3,end)),'eps')
       ext = 'eps';
end

%%% fmt = fmt( 1 : min(3,size(fmt,2)) );

%**************************************************************
% Check FileName for Directory

pfad = '';

if isempty(file)

   flt = sprintf('*.%s',ext);

   [file,pfad]= uiputfile( flt , 'Save Image ... ' );

   if isequal(file,0) | isequal(pfad,0) | isempty(file)
      file = '';
      return
   end

   file = fullfile(pfad,file);

end

if ~isequal(file,-1)

    [pfad,n,e] = fileparts(file);

    if ~isempty(pfad)
        ok = ( exist(pfad,'dir') == 7 );
        if ok
           ok = ~isempty(dir(pfad));
        end
        if ~ok
            msg = sprintf('Directory "%s" doesn''t exist.',pfad);
        end
    end

    if isempty(n) & isempty(e)
       file = -1;
    end

end
        
%**************************************************************
% Default FileName == FigureName

if isequal(file,-1)

   name = get(fig,'name');

   if isempty(name)

      file = sprintf('figure_%2.2d',round(fig));

   else 

      sp = '_';
      dp = double(sp);

      name = double(rmblank(name,2));

      % 0 .. 9  |  A .. Z   |  a .. z  |  .  |  _  | ~

      ok = ( ( ( 48 <= name )  &  ( name <= 57  ) ) | ...
             ( ( 65 <= name )  &  ( name <= 90  ) ) | ...
             ( ( 97 <= name )  &  ( name <= 122 ) ) | ...
               ( name == 46 )  |  ( name == 95  )   | ( name == 126 ) );

      if ~all(ok)
          name(find(~ok)) = dp;
          if size(name,2) > 2
             ok = find( name == dp );
             ok = ok( find( diff(ok) == 1 ) + 1 );
             name(ok) = [];
          end
          file = rmblank(char(name),2,sp);
      else
          file = char(name);
      end

   end

   file = cat(2,file,'.',ext);

   if ~isempty(pfad)
       file = fullfile(pfad,file);
   end

end

%**************************************************************
% Check FileName with Format


if ~any( ( file == filesep ) | ( file == '.' ) )
    file = cat(2,file,'.',ext);
else
   [p,n,e] = fileparts(file);
   if isequal(e,'.*')
      file = fullfile(p,cat(2,n,'.',ext));
   end
end

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function msg = menu_clb(h);

% MENU_CLB  CallBack of Menu

opt = get(h,'userdata');

t = '';
while ~strcmp(t,'figure')
      h = get(h,'parent');
      t = get(h,'type');
end

if ~iscell(opt)
    opt = {opt};
end


fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end

msg = feval(fcn,h,'',opt{:});

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,cnf,hm] = make_gui(vin)

% MAKE_GUI  Returns GUI_Structure / Creates Menu's
%
% [ParentHandle] , ['&Label'] , ['#FieldName'] , ['CallBackFcn'] , ...
%  ScaleVec , i*Margin
% 
% Label | FieldName   ==> SubStructure with Label
% ParentHandle        ==> Create in Parent
%
% Default CallBackFcn: SAVEIMG
%

msg = '';
cnf = [];
hm  = [];

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end

par = [];

clb = '';
lab = '';
fld = '';

scl = -1;

off = [];

nin = prod(size(vin));

ind = ( 1 : nin ) + 1;

%**************************************************
% Check for ParentHandle in Input

if ~( nin == 0 )

    for typ = { 'uimenu' 'uicontextmenu' 'figure' }
        [h,vin,nr] = hndlvarg(vin,[],typ{1});
        if ~( nr == 0 )
            break
        end
    end

    if ~( nr == 0 )
        par = h;
        nin = nin - 1;
        ind(nr) = [];
    end

end

%**************************************************
% Check for Function, Label, FieldName and Scale

msg = cell(0,1);

for ii = 1 : nin

    v = vin{ii};

    is_str =  chkstr(v,1);
    is_num = isnumeric(v);

    %---------------------------------------------------
    if is_str
    %---------------------------------------------------

       ok = ( size(v,2) > 1 );

       is_lab = ( any( v == '&' ) & ok );
       is_fld = ( ( v(1) == '#' ) & ok );

       is_clb = ~( is_lab | is_fld );
       
       is_clb = ( is_clb & isempty(clb) );
       is_lab = ( is_lab & isempty(lab) );
       is_fld = ( is_fld & isempty(fld) );


       if     is_clb
           if ~( exist(v,'file') == 2 )
               str = sprintf('CallBackFcn "%s" in %.0f. Input doen''t exist.',v,ind(ii));
               msg = cat( 1 , msg , {str} );
           else
              clb = v;
           end       
       elseif is_lab
              lab = v;
       elseif is_fld
              fld = v(2:end);
       else
           str = sprintf('Unrecognized %.0f. StringInput: %s',ind(ii),v);
           msg = cat( 1 , msg , {str} );
       end

    %---------------------------------------------------
    elseif is_num
    %---------------------------------------------------

       p = prod(size(v));

       ok = ( p == 1 );
       if ok
          ok = ( real(v) == 0 ); 
       end

       if ok

          if ~isempty(off)
              str = sprintf('Unrecognized %.0f. NumericInput.',ind(ii));
              msg = cat( 1 , msg , {str} );
          else
              v = imag(v);
              if ( isfinite(v) & ( v > 0 ) )
                 off = v;
              else
                 str = 'Value for Margin in %.0f. Input must be a finite positive.';
                 str = sprintf(str,ind(ii));
                 msg = cat(1,msg,{str});
              end
          end

       elseif ~isequal(scl,-1)

           str = sprintf('Unrecognized %.0f. NumericInput.',ind(ii));
           msg = cat( 1 , msg , {str} );

       else
           ok = isempty(v);
           if ~ok
               v = v(:)';
              ok = all( isfinite(v) & ( 0 < v ) & ( v <= 10 ) );
           end

           if ~ok
               str = 'Values for Scale in %.0f. Input must be finite positive, smaller 10.';
               str = sprintf(str,ind(ii));
               msg = cat(1,msg,{str});
           else
               scl = v;
           end 

       end

    %---------------------------------------------------
    else
    %---------------------------------------------------

        str = 'Inputs must be a ParentHandle, nonempty Strings or a Vector for Scale.';
        str = sprintf('Invalid %.0f. Input.\n%s',ind(ii),str);
        msg = cat( 1 , msg , {str} );

    %---------------------------------------------------
    end
    %---------------------------------------------------

end

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    return
end

msg = '';

%**************************************************

if isempty(clb)
   clb = fcn;
end

if     isempty(lab) & ~isempty(fld)
   lab = fld;
elseif isempty(fld) & ~isempty(lab)
   fld = upper(clb);
end

cnf = menu_gui(clb,lab,fld,scl,off);

if isempty(par)
   return
end


%---------------------------------------
% Check for first Seperator !!!

if ~isempty(get(par,'children'))
    fld = fieldnames(cnf);
    val =   getfield(cnf,fld{1});
    val{2} = 1;
    cnf =   setfield(cnf,fld{1},val);
end

%---------------------------------------
% Create Menus

hm = make_menu(par,'',cnf,get(par,'tag'));


%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function HM = make_menu(par,group,cnf,tag);

% MAKE_MENU  Creates UIMenus from Structure
%
% MAKE_MENU( ParentHandle ,  Group , Config , Tag );
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
    if isempty(tag)
       tag = group;
    else
       tag = cat( 2 , tag , '.' , group );
    end
end

for ii = 1 : n

  cc = getfield( cnf , field{ii} );

  CB    = '';
  usd   = [];

  if chkstr(cc{3},1)
     CB = sprintf('%s(gcbo,''%s'',''%s'',1);',cc{3},group,field{ii});
  end

  if prod(size(cc)) >= 4
     usd = cc{4};
  end

  tg = field{ii};
  if ~isempty(tag)
      tg = cat( 2 , tag , '.' , tg );
  end

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
       'tag'             , tg       , ...
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
       'tag'             , tg       , ...
       'Interruptible'   , 'off'    , ...
       'BusyAction'      , 'cancel'  );
 
     ch = cc{2};

   end 

  if isstruct(ch)
  % Children

     ch = make_menu(HC(ii),field{ii},ch,tag);

  else

     ch = [];

  end

  setappdata( HC(ii) , 'Children' , ch );

  HM = setfield( HM , field{ii} , HC(ii) );

end


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,lf,pf,fig] = pointscale(Nout,varargin);

% POINTSCALE(Scale,PointIntervall,MaxPoint,Plot)
%
% Returns PixelWidth for PointsUnits
%

msg = '';

lf = zeros(0,2);  % [ Linewidths  Pixels ]
pf = zeros(0,4);  % [ Pixel MeanLW MinLW MaxLW  ]  Unique Pixels

fig = [];

%-----------------------------------------------------------

frm = '.png';

tmp = cat(2,tempname,frm);

%-----------------------------------------------------------

Nin  = nargin - 1;

if Nin < 1, scl = 1.0; else, scl = varargin{1}; end
if Nin < 2, dlw = 0.1; else, dlw = varargin{2}; end
if Nin < 3, lwx = 5.0; else, lwx = varargin{3}; end
if Nin < 4, plt = [];  else, plt = varargin{4}; end

ok = ( isnumeric(scl) & ( prod(size(scl)) == 1 ) & ...
       isnumeric(dlw) & ( prod(size(dlw)) == 1 ) & ...
       isnumeric(lwx) & ( prod(size(lwx)) == 1 )       );
if ok
   ok = ( isfinite(scl) & ( scl > 0 ) & ...
          isfinite(dlw) & ( dlw > 0 ) & ...
          isfinite(lwx) & ( lwx > 0 )       );
end

if ~ok
    msg = 'Inputs must be positive finite numerics.';
    return
end

if isempty(plt)
   plt = ( Nout == 0 );
else
   plt = ~isequal(plt,0);
end

%-----------------------------------------------------------

ppi = 72;   % PointPerInch

spi = get(0,'screenpixelsperinch');

if dlw >= lwx
   lw = dlw;
else
   lw = ( dlw : dlw : lwx )';
end

nl = size(lw,1);

lp = lw / ppi * spi;

lp = ceil(lp); % estimated PixelWidth of Single Line on Screen

lp = lp + 2*min(max(lp,5),10);     % Space for Line

fp = [ 100  sum(lp) ];     % FigurePosition

lp = cumsum(lp,1) - lp/2;  % CenterPos of Line

%-----------------------------------------------------------
% Plot Lines, Save Figure into TempFile

fig = figure( 'units'    , 'pixels' , ...
              'position' , [ 10 10  fp ] , ...
              'visible'  , 'off'    , ...
              'color'    , [ 1  1  1 ]   );

axe = axes( 'parent' , fig , ...
             'units' , 'normalized' , ...
          'position' , [ 0 0 1 1 ]  , ...
              'xlim' , [ 0   1   ]  , ...
              'ylim' , [ 0 fp(2) ]  , ...
              'ydir' , 'reverse'    , ...
           'visible' , 'off'        , ...
          'nextplot' , 'add'        , ...
             'color' , 'none'             );

for ii = 1 : nl

    plot( [0 1] , lp([ii ii]) , 'k-' , ...
          'linewidth',lw(ii),'marker','none','parent',axe);

end

try
   msg = saveimg(fig,tmp,scl,frm);
catch
   msg = sprintf('Error call SAVEIMG( %s ).\n%s',tmp,lasterr);
end

delete(fig); fig = [];

if ~isempty(msg)
    return
end

%-----------------------------------------------------------
% Read Image, get PixelWidth's for Lines

try
   c = imread(tmp);
catch
   msg = sprintf('Error call IMREAD( %s ).\n%s',tmp,lasterr);
end

try, delete(tmp); end

if isempty(msg)
   if isempty(c)
      msg = sprintf('Empty Image written by SAVEIMG in "%s.',tmp);
   end
end

if isempty(msg)

   ic = floor(size(c,2)/2);

   c  = ( c(:,ic,1) == 0 );

   if ~any(c)

       msg = 'Can''t find Lines in Image';

   else

       c  = find(c);

       nc = size(c,1);

       c = cat( 1 , 1 , find( diff(c,1,1) > 1 )+1 , nc+1 );

       pw  = diff(c,1,1);

       if ~isequal(size(pw),size(lw))
           msg = 'Invalid Number of Lines in Image.';
       end

   end

end

if ~isempty(msg)
    return
end

[px,lm,nn,ld,l0,l1] = grpmean(pw,lw,1);


lf = [ lw pw ];
pf = [ px lm l0 l1 ];

if ~plt
    return
end

%-----------------------------------------------------------
% Show Result

fig = figure;

axe = axes( 'parent'   , fig , ...
            'ylim'     , [ 0  max(pw)+1  ] , ...
            'xlim'     , [ 0  lw(nl)+dlw ] , ...
            'box'      , 'on'              , ...
            'tickdir'  , 'out'             , ...
            'nextplot' , 'add'    );

fs = get(axe,'fontsize');

title( sprintf('Points <---> Pixels  at  Scale %.3g with %.0f px/in',scl,spi) , ...
       'fontsize' , fs+2 );

xlabel(sprintf('Points   Intervall: %.3g',dlw),'fontsize',fs+1)
ylabel('Pixels','fontsize',fs+1)

plot( lw , pw , 'b-' , 'parent' , axe );

plot( lm , px , 'b.' , 'parent' , axe );

for ii = 1 : size(px,1)
    text( l1(ii) , px(ii)+0.1 , ...
          sprintf('%.0f: %.2f  ',px(ii),lm(ii))     , ...
          'fontsize'   , max(8,fs-2) , ...
          'fontweight' , 'normal'   , ...
          'horizontal' , 'right'    , ...
          'vertical'   , 'bottom'   , ...
          'parent'     , axe );
end


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,cmp,wc,nn,cc] = strwcmp(str,cmp,wc,nn,cc)

% STRWCMP   Compare Strings, including WildCards
%
% OK = STRWCMP( STR , CMP , WildCard )
%
% STR: CharacterArray or CellStringArray to compare with Comp
%
% CMP: CharacterArray or CellStringArray with strings to compare with STR
%         strings can contains WildCards
%
% WildCard: specify WildCard to use, default: '*'
%
% OK : logical Array with same size of STR, contains 1 if 
%        any strings of CMP match string of STR
%
% Special Output:  [ok,cmp,wc,nn,cc] = STRWCMP( STR , CMP , [WildCard] )
%
%  to use in follwing statements: ok = STRWCMP( STR , cmp , wc , nn , cc );
%
%  which makes it a bit faster.
%
% see also: STRCMP, FINDSTR
%

Nin  = nargin;
Nout = nargout;

Msg = '';
 nl = char(10);

%***************************************************************
% Check Inputs

%---------------------------------------------------------------
% String

if Nin < 1
   str = cell(0,0);
end

if ~( iscell(str) & isempty(str) )
   [ok,str] = chkcstr(str,0);
   if ~ok
      Msg=[ Msg  nl(1:(end*(~isempty(Msg)))) ...
           'First Input must be a CharacterArray or CellStringArray.' ];
   end
end

%---------------------------------------------------------------
% CompareString

if Nin < 2
   cmp = cell(0,0);
end

if ~( iscell(cmp) & isempty(cmp) )
   [ok,cmp] = chkcstr(cmp,0);
   if ~ok
      Msg=[ Msg  nl(1:(end*(~isempty(Msg)))) ...
            'Second Input must be a CharacterArray or CellStringArray.' ];
   end
end

%---------------------------------------------------------------
% WildCard

if Nin < 3
   wc = '*';
elseif ~( ischar(wc) & ( prod(size(wc)) == 1 ) )
   Msg=[ Msg  nl(1:(end*(~isempty(Msg)))) ...
         'WildCard must be a Single Character.' ];
end
  
%---------------------------------------------------------------

if ~isempty(Msg)
    error(Msg)
end

%***************************************************************

si = size(str);

if ( isempty(str) | isempty(cmp) ) & ( Nout <= 1 ) 
   ok = zeros(si);
   return
end

cmp = cmp(:);

if ~isempty(cmp)
    if any(strcmp(cmp,wc)) & ( Nout <= 1 )  % Wildcard only
       ok = ones(si);
       return
    end
end

%***************************************************************
% Analyze CompareStrings

nc = size(cmp,1);

ok = ( Nin == 5 );

if ok
   ok = ( isequal(size(cc),[nc 1]) & iscell(cc) & ...
          isequal(size(nn),[nc 2]) & isnumeric(nn)    );
   if ok
      try
         ok = cat( 1 , cc{:} ); 
         ok = ( size(ok,2) == 3 );
      catch
         ok = 0;
      end
   end
end


%--------------------------------------------------
if ~ok
%--------------------------------------------------

  cc    = cell(nc,1);
  cc(:) = { zeros(0,3) };  % { [ Start End N  ] }

  nn = zeros(nc,2);        %   [ Ncmp  sum(N) ] 

  for ii = 1 : nc

    if ~isempty(cmp{ii})

       iwc = ( double(cmp{ii}) == double(wc) );

       if ~any( iwc )

          nn(ii,:) = size(cmp{ii},2);
          cc{ii}   = [ 1  nn(ii,:) ];

       else

         %--------------------------------------------
         % Remove Duplicate WildCards

         iwc = find( iwc );
         if ~( prod(size(iwc)) == 1 )
            jj = find( diff(iwc) == 1 );
            if ~isempty(jj)
               cmp{ii}(iwc(jj+1)) = [];
            end
         end

         %--------------------------------------------
         % Get Start End
     
         iwc = ( double(cmp{ii}) == double(wc) );
  
          n  = size(iwc,2);

          if ( n == 1 ) & ( iwc == 1 ) & ( Nout <= 1 ) % Wildcard only
             ok = ones(si);
             return
          end

          i0 = ~iwc(1);
          i1 = ~iwc(n);

         iwc = cat( 2 , ones(1,i0) , iwc , ones(1,i1) );

         iwc = find( iwc );

         iwc = iwc(:);

           n = size(iwc,1) - 1;

         cc{ii} = zeros(n,3);

         if n > 0      

            cc{ii}(:,[1 2]) = cat( 2 , iwc(1:n)+1 , iwc((1:n)+1)-1 ) - i0;

            cc{ii}(:,3) = cc{ii}(:,2) - cc{ii}(:,1) + 1;
 
         end

         nn(ii,:) = cat( 2 , size(cmp{ii},2) , sum(cc{ii}(:,3),1) );

       end

    end

  end

%--------------------------------------------------
end
%--------------------------------------------------

if ( Nout > 1 )

  if ( isempty(str) | isempty(cmp) )
     ok = zeros(si);
     return
  end

  if any(strcmp(cmp,wc))  % Wildcard only
     ok = ones(si);
     return
  end

end

%***************************************************************
% Compare

ok = zeros(si);

for ii = 1 : prod(si)
 
    s2 = size(str{ii},2);

    for jj = 1 : nc

        ok(ii) = ( ( s2 == 0 ) & ( nn(jj,1) == 0 ) );
       
        if ok(ii)
           break
        end
       
        if ( s2 >= nn(jj,2) ) & ~( nn(jj,1) == 0 )

           ok(ii) = compare(str{ii},cmp{jj},cc{jj});

           if ok(ii)
              break
           end
            
        end

    end

end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = compare( str , cmp , cc )

sc = size(cmp,2);

ok = 1;

for ii = 1 : size(cc,1) 

    s2 = size(str,2);

    ok = ( ok &  ( s2 >= cc(ii,3) ) );

    if ~ok
       break
    end

    ic  = ( cc(ii,1) : cc(ii,2) );

    i01 = ( cc(ii,[1 2]) == [ 1  sc ] );

    i0  = ( 1 : cc(ii,3) );
    i1  = s2 - cc(ii,3) + i0;

    i2  = cat( 2 , findstr( str , cmp(ic) ) , 0 );
    i2 = i2(1);
    i2 = ( i2 + cc(ii,3) - 1 ) * ~( i2 == 0 );

    i2  = s2       * (  i01(1) &  i01(2) & strcmp(str    ,cmp    ) ) + ...
          cc(ii,3) * (  i01(1) & ~i01(2) & strcmp(str(i0),cmp(ic)) ) + ...
          s2       * ( ~i01(1) &  i01(2) & strcmp(str(i1),cmp(ic)) ) + ...
          i2       * ( ~i01(1) & ~i01(2) );

    ok = ( ok & ~( i2 == 0 ) );

    if ~ok
       break
    end

    str = str( i2+1 : s2 );

end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [h,v,nr] = hndlvarg(v,n,typ,tag);

% HNDLVARG  Checks Input-varargin for Handle
%
%  [H,V,InputNr] = HANDLVARG( V , n , Type , Tag )
%
%  Search in V{n(1)} .. V{n(2)} for Handle of Type and Tag
%
%  If no Handle found and Type =  'figure' | 'axes', 
%    GCF | GCA is returned 
%
%  see also: CHKHNDL
%

Nin = nargin;

h  = [];
nr = 0;

if Nin == 0
   v = cell(1,0);
   return
end


if ~iscell(v)
   error('Input must be a CellArray.');
end

nv = prod(size(v));

if Nin < 2
   n = [];
end

if Nin < 3
   typ = 'figure';
end

chkin = { typ };

if Nin == 4
   chkin = cat( 2 , chkin , {tag} );
end

%------------------------------------------------------------

if isempty(n)
   n = [ 1  nv ];
else
   n = n(:)';
   n = n( 1 : max(1,size(n,2)) );
   n = min( n , nv );
   if size(n,2) == 1
      n = [ 1  n ];
   end
end


%-------------------------------------------------------------

ok = 0;

for nr = n(1) : n(2)

   ok = chkhndl(v{nr},chkin{:});

   if ok
      break
   end

end

%-------------------------------------------------------------

if ok

  h = v{nr};

  v(nr) = [];   % Reduce v

  return

end

nr = 0;

if ~any( strcmp( typ , { 'figure' 'axes' } ) )
   return
end

h = get( 0 , 'currentfigure' );
if isempty(h)
   return
end

if strcmp( typ , 'axes' )
   h = get( h , 'currentaxes' );
end

ok = chkhndl(h,chkin{:});

if ~ok
   h = [];
end


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [ok,msg] = chkhndl(h,typ,tag);

% CHKHNDL(H,Type,Tag)  Checks, if H is a Handle of specified Type and Tag
%
%  Tag   CharArray or CellStringArray which Tags, the Handle has to be.
%

Nin = nargin;

ok  = 0;
msg = [];


if Nin == 0
   return
end

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

is_typ = chkstr(typ,0);

if is_typ

  ok = isempty(typ);
  if ~ok
      ok = strcmp( get(h,'type') , typ  );
  end

  if ~ok | ( Nin < 3 )
     return
  end

elseif ( Nin == 3 )

  msg = 'Input Type must be a String.';
  ok  = 0;
  return

else
 
  tag = typ;
  
end


%-------------------------------------------------------------------------
% Check Tag

[ok,tag] = chkcstr(tag,0);

if ~ok
    msg = 'Input Tag must be a CharArray or CellStringArray.';
    return
end   


%-------------------------------------------------------------------------
% Check with Tag

 t = get(h,'tag');

nt = size(t,2);

ok = strcmp(t,tag);

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

function Config = menu_gui(fcn,lab,fld,scl,off)


qlt = [ ( 100 : -5 :  70 )  60  50 ];  % JPG-Quality

is_scl = ~isempty(scl);

if is_scl
   if isequal(scl,-1)
      scl = [ 50  75  90  110  125  150  200 ]; % default PrintScale
   else
      scl = scl * 100;  % PrintScale in Prozent !!!
   end
end

if isempty(off)
   off = cell(1,0);
else
   off = { off*i };
end

is_win = strcmp( upper(computer) , 'PCWIN' );

%------------------------------------------------------------------

ClipBoard = cell(1,0);

if is_win

   ClipBoard = struct( 'Meta' , { { '&MetaFile'  0  fcn  {'.meta'}   } } , ...
                       'Bmp'  , { { '&BitMap'    0  fcn  {'.bitmap'} } }       );

   ClipBoard = { 'ClipBoard'  , { { 'Clip&Board'     1   ClipBoard      } } };

end

%------------------------------------------------------------------
% PostScript

ps  = struct( 'RGB'  , { { '&RGB'  0  fcn  [{'.psc2'         } off ] } } , ...
              'CMYK' , { { '&CMYK' 0  fcn  [{'.psc2'  '-cmyk'} off ] } }        );

eps = struct( 'RGB'  , { { '&RGB'  0  fcn  [{'.epsc2'        } off ] } } , ...
              'CMYK' , { { '&CMYK' 0  fcn  [{'.epsc2' '-cmyk'} off ] } }        );

PostScript = struct(  'PS' , { {  '&PS'  0    ps } } , ...
                     'EPS' , { { '&EPS'  0   eps } }       );

%------------------------------------------------------------------
% Capture

 jpg = jpg_gui('*',qlt,fcn,{});

Capture = struct( 'JPG' , { { '&JPEG'  0   jpg           } } , ...
                  'BMP' , { { '&BMP'   0   fcn  {'*bmp'} } } , ...
                  'TIF' , { { '&TIFF'  0   fcn  {'*tif'} } } , ...
                  'PNG' , { { 'P&NG'   0   fcn  {'*png'} } } , ...
                  'PPM' , { { 'P&PM'   0   fcn  {'*ppm'} } } , ...
                  'XPM' , { { '&XPM'   0   fcn  {'*xpm'} } } , ...
                  'PCX' , { { 'P&CX'   0   fcn  {'*pcx'} } } , ...
                  'HDF' , { { '&HDF'   0   fcn  {'*hdf'} } } , ...
                  'XWD' , { { 'X&WD'   0   fcn  {'*xwd'} } }       );

%------------------------------------------------------------------
% Image

 jpg = jpg_gui('.',qlt,fcn,{});

Image = struct( 'JPG' , { { '&JPEG'  0   jpg           } } , ...
                'BMP' , { { '&BMP'   0   fcn  {'.bmp'} } } , ...
                'TIF' , { { '&TIFF'  0   fcn  {'.tif'} } } , ...
                'PNG' , { { 'P&NG'   0   fcn  {'.png'} } } , ...
                'PPM' , { { 'P&PM'   0   fcn  {'.ppm'} } }       );

%------------------------------------------------------------------
% Scaled

Scale = cell(1,0);

if is_scl

   jpg = jpg_gui('.',qlt,fcn,{});

   Scale = { 'Print'      , { { '&Print'         0   fcn   [{'.dlg'} off ] } } , ...
             ClipBoard{:} , ...
             'PostScript' , { { 'Post&Script'    1   PostScript     } } , ...
             'JPG'        , { { '&JPEG'          1   jpg            } } , ...
             'BMP'        , { { '&BMP'           0   fcn   {'.bmp'} } } , ...
             'TIF'        , { { '&TIFF'          0   fcn   {'.tif'} } } , ...
             'PNG'        , { { 'P&NG'           0   fcn   {'.png'} } } , ...
             'PPM'        , { { 'P&PM'           0   fcn   {'.ppm'} } }       };

   Scale = scl_gui( struct(Scale{:}) , scl );

   Scale = { 'Scale' , { { 'Sca&led'  1  Scale } } };

end


%------------------------------------------------------------------
% Print-Menu

 Config = { 'Print'      , { { '&Print'         0   fcn   [{'.dlg'} off] } } , ...
            ClipBoard{:} , ... 
            'PostScript' , { { 'Post&Script'    1   PostScript     } } , ...
            'Capture'    , { { 'S&creenShot'    1   Capture        } } , ...
            'Image'      , { { '&Image'         0   Image          } } , ...
            Scale{:}         };

 Config = struct(Config{:});

if ~isempty(lab)
    try
       Config = struct( fld        , { { lab  0  Config } } );
    catch
       Config = struct( upper(fcn) , { { lab  0  Config } } );
    end      
end
  
%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function cnf = jpg_gui(pre,qlt,fcn,opt)


n = prod(size(qlt));

cnf = cell(2,n);

for ii = 1 : n

    cnf{1,ii} = sprintf('Q%2.2d_%3.3d',ii,qlt(ii));

    lab = sprintf( 'Quality %.0f%%'  , qlt(ii) );
    usd = sprintf( '%sjpg%.0f' , pre , qlt(ii) );

    cnf{2,ii} = { { lab  0  fcn  [ {usd} opt ] } };

end

cnf = struct(cnf{:});

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function cnf = scl_gui(cfg,scl);

n = prod(size(scl));

cnf = cell(2,n);

for ii = 1 : n

    cnf{1,ii} = sprintf('S%2.2d_%3.3d',ii,scl(ii));

    lab = sprintf( '%.0f%%'  , scl(ii) );

    sep = ( ( ii == 2 ) & ( scl(1) == 100 ) );

    cnf{2,ii} = { { lab  sep  add_scl(cfg,scl(ii)/100) } };

end

cnf = struct(cnf{:});

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function cnf = add_scl(cnf,scl);

f = fieldnames(cnf);

for ff = f(:)'

    v = getfield(cnf,ff{1});

    if isstruct(v{3})

       v{3} = add_scl(v{3},scl);

    elseif iscell(v{4})

       v{4} = cat( 2 , v{4} , {scl} );

    end

    cnf = setfield(cnf,ff{1},v);

end

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

