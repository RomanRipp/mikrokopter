function [cmap,msg] = ferretc(name,N);

% FERRETC  returns Ferret-ColorMap
%
% [ ColorMap , Msg ] = FERRETC( Name , ColorNumber);
%
%  reads Ferret-Colormap from <Name>.spk
%
% A ColorMap of an SPK-File has 4 Columns:
%
%   [ Number Red  Green  Blue ]
%
% All Values must have a Range between 0 and 100 !
%
% Use a NonZero imaginary Part of ColorNumber to return
%  a GreyScaled ColorMap:  Grey = 0.299*R + 0.587*G + 0.114*B
%
%----------------------------------------------------
% Info:
%
% [ Info , Msg ] = FERRETC( 'list' , [ColorNumber] )
%
% Info = { Name FileName [ColorMap] }  % 2- or 3-Column-CellArray
%
%----------------------------------------------------
% Demo:
%
% [Fig,Msg] = FERRETC( 'test' , [ColorNumber] )
%
%  Shows a ExampleFigure with Ferret-Colormaps 
%   from Private-Directory
%
%

Nin  = nargin;
Nout = nargout;

if Nin < 1
   name = '';
end

if Nin < 2
   N = [];
end

FerretPath = '/usr/local/ferret/ppl/';
FerretExt  = '.spk';

cmap = zeros(0,3);

is_grey = 0;

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end

app = upper(fcn);

%******************************************************
% Check Inputs

msg = cell(0,1);

is_cb = isnumeric(name);

if is_cb

   ok = ( prod(size(name)) == 1 );

   if ok
      ok = ishandle(name);
      if ok
         tag = get(name,'tag');
         ok  = strcmp(tag(1:min(size(app,2),size(tag,2))),app);
      end
   end

   if ~ok
       msg = cat( 1 , msg , ...
             {sprintf('1. Numeric Input must be a %s-Handle.',app)} );
   end

   if ~chkstr(N,1)
       msg = cat( 1 , msg , ...
             {'2. Input after 1. Numeric Input must be a nonempty String.'} );
   end

elseif isempty(name)

   name = 'default';

elseif ~chkstr(name,1)

   msg = cat( 1 , msg , {'Name must be a String.'} );

end

if ~isempty(N) & ~is_cb

    ok = ( isnumeric(N) & ( prod(size(N)) == 1 ) );

    if ok
       is_grey = ~( imag(N) == 0 );
       N  = real(N);
       ok = ( ( mod(real(N),1) == 0 ) & ( real(N) > 0 ) );
    end 

    if ~ok
       msg = cat( 1 , msg , {'N must be a positive Integer.'} );
    end

end

%------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    msg = sprintf('Invalid Inputs.\n%s',msg);
    if Nout < 2
       error(msg);
    end
    return
end

msg = '';

%******************************************************
% Check for CallBack

if is_cb

   ferret_cb(name,N,app);

   return

end

%******************************************************
% Test or List !!!

if any( strcmp( name , { 'test' 'list' } ) )

   is_list = strcmp(name,'list');

   if isempty(N) & ~is_list
      N = 32;
   end

   [msg,cmap] = ferret_list(N+i*is_grey*is_list,FerretPath,FerretExt,name);

   if ~( isempty(msg) | is_list )
       cmap = [];
   end

   if     ~isempty(msg) & ( Nout < 2 )
           error(msg);
   elseif ~isempty(msg) | is_list
           return
   end
   
   cmap = ferret_test(N,cmap,is_grey,fcn,app);

   return

end

%******************************************************

if isempty(N) 

   fig = get( 0 , 'currentfigure' );

   if isempty(fig)
      c = get( 0 , 'defaultfigurecolormap' );
   else
      c = get( fig , 'colormap' );
   end
  
   N = size(c,1);

end

%******************************************************
% Load SPK-File

FerretExt  = FerretExt( 1+strcmp(name(end),'.') : end );
 
%------------------------------------------------------
% Try different Combinations for FileName

 file = { [            name FerretExt ]
          [            name           ]
          [ FerretPath name FerretExt ]
          [ FerretPath name           ] };


cc = zeros(0,3);
ok = 0;

cm = NaN;  % ColorMax  [ 100 | 255 ]

for ff = file(:)'

 ok = ( exist(ff{1},'file') == 2 );

 if ok
 
    try
       cc = load(ff{1});
    end

    ok = ( isnumeric(cc) & ( size(cc,1) >= 2 ) & ( size(cc,2) == 4 ) ); 

    if ok
       [h,si] = sort(cc(:,1));
          cc  = cc(si,:);
          ok  = all( diff(cc(:,1)) > 0 );
    end

    if ok

       cm = cc(:,[2 3 4]);

       cm = max(cm(:));

       ok = ( cm <= 255 );

       if ok

          cm = 100 + 155 * ( cm > 100 );

       end

    end

 end

 if ok
    break
 end

end

%------------------------------------------------------

if ~ok
 
    msg = sprintf('Can''t find valid FerretColormap for "%s".',name);

    if Nout < 2
       error(msg);
    end
    return

end

%******************************************************

ind  = linspace(cc(1,1),cc(end,1),N)';

cmap = interp1(cc(:,1),cc(:,2:4),ind)/cm;

if is_grey

   cmap = 0.299*cmap(:,1) + 0.587*cmap(:,2) + 0.114*cmap(:,3);
   cmap = cmap(:,[1 1 1]);

end

%******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,list] = ferret_list(nc,FerretPath,FerretExt,mode)

% Returns available SPK-Files from PRIVATE-Directroy 
%   or MatlabSearchPath
%
% LIST = { Name File [ColorMap] }
%

Nin  = nargin;
Nout = nargout;

if Nin <1
   nc = [];
end

msg = '';

is_cmap = ~isempty(nc);
is_test = strcmp(mode,'test');

list = cell(0,2+is_cmap);
 
%******************************************************
% Settings

ext = FerretExt( ( 1 + ( FerretExt(1) == '.' ) ) : end );

wcd = cat(2,'*',FerretExt);

%******************************************************
% Get SPK-Files from PRIVATE-Directory
%                 or FerretPath

src = fileparts(which(mfilename));
src = fullfile( src , 'private' );

pfad = { src  FerretPath  '' };  % Search in: Private FerretPath MatlabPath

np =  size(pfad,2);
nn = zeros(np,1+is_cmap);

for ii = 1 : np

    n = cell(0,1);
    f = cell(0,1);

    %--------------------------------------------------
    if isempty(pfad{ii}) 
    %--------------------------------------------------

       [n,f] = whichfile(FerretExt);    % MatlabPath

    %--------------------------------------------------
    else
    %--------------------------------------------------

        if ( exist(pfad{ii},'dir') == 7 )

           d = dir(fullfile(pfad{ii},wcd));

           if ~isempty(d)
               d = d(find(~cat(1,d.isdir)));
               if ~isempty(d)
                   n = {d.name};
                   n = n(:);
               end
           end

        end

        if ~isempty(n)
            f    = cell(size(n));
            f(:) = pfad(ii);
        end

    %--------------------------------------------------
    end
    %--------------------------------------------------

    nn(ii,1) = size(n,1);

    if ~( nn(ii,1) == 0 )

        for jj = 1 : nn(ii,1)
            f{jj} = fullfile(f{jj},n{jj});
        end

        if ~is_cmap

            l = [ n  f ];

        else

            [c,ok,n] = ferret_cmap(n,f,nc);
            nn(ii,2) = sum(ok);
            if     nn(ii,2) == 0
                   l = cell(0,3);
            elseif nn(ii,2) < nn(ii,1)
                   ok = find(ok);
                   l = [ n(ok) f(ok) c(ok) ];
            else
                   l = [ n  f  c ];
            end

        end

        list = cat(1,list,l);

        if is_test & ~( nn(ii,1+is_cmap) == 0 )
           break  % At Private, FerretPath or MatlabPath
        end

    end

end

%******************************************************

if ~any(nn(:,1))
    msg = sprintf('No %s-Files found.',upper(ext));
elseif is_cmap
    if ~any(nn(:,2))
        msg = sprintf('No valid %s-Files found.',upper(ext));
    end
end

%******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [cmap,ok,name] = ferret_cmap(name,file,nc);


cmap =  cell(0,1);
ok   = zeros(0,1);

if isempty(name)
   return
end

n = size(name,1);

cmap    = cell(n,1);
cmap(:) = { zeros(0,3) };

ok = zeros(n,1);

for ii = 1 : n

    [cmap{ii},msgf] = ferretc( file{ii} , nc );

    ok(ii) = isempty(msgf);

    if ok(ii)
       i1 = find( double(name{ii}) == double('.') );
       if ~isempty(i1)
           name{ii} = name{ii}( 1 : max(i1)-1 );
       end
    end

end


%******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function fig = ferret_test(nc,nfc,is_grey,fcn,app)

% FERRET_TEST    Shows Colormaps of PRIVATE-Directory
%
% NFC = { Name File ColorMap }
%

wdt = 30;   % Top/Bottom-Space and TitleWidth

fig = [];

n = size(nfc,1);

[h,si] = sort(nfc(:,1));

nfc = nfc(si,:);

%******************************************************
% FigureSettings

FontName   = 'helvetica';
FontUnits  = 'points';
FontSize   =  10 - 2*strcmp(computer,'PCWIN');
FontWeight = 'normal';
FontAngle  = 'normal';

xd = [ 1  nc ];
yd = [ 1  n  ];

cc = cat(3,nfc{:,3});
cc = permute(cc,[3 1 2]);

xl = xd + 0.5 * [ -1  1 ];
yl = yd + 0.5 * [ -1  1 ];

xt = nc/2 + 0.5;
yt = ( yl(1)+1 : 1 : yl(2)-1 );

if ~( mod(nc,2) == 0 )
   xt = [];
end

%******************************************************
% Create Figure and Axes

fig = figure('units'            , 'pixels'    , ...
             'position'         , [ 100 100 200 200 ] , ...
             'paperorientation' , 'portrait'  , ...
             'paperunits'       , 'inches'    , ...
             'color'            , [ 1  1  1 ] , ...
             'name'             , ['Demo ' app] , ...
             'numbertitle'      , 'off'       , ...
             'menubar'          , 'none'      , ...
             'toolbar'          , 'none'      , ...
             'tag'              , [app '_Figure']       );
 
papsi = get(fig,'papersize');

pappos = [0.8 1.5 papsi-2*[ 0.8  1.5 ] ];

set(fig,'paperposition',pappos)

wysiwyg
drawnow

%------------------------------------------------------
% Fit FigurePosition vertical

ssi = get(0,'screensize');
ssi = ssi(4);                   % ScreenHeight

figpos = get(fig,'position');

if figpos(2)+figpos(4)+wdt > ssi
   figpos(2) = ssi - wdt - figpos(4);
   dev = min( 0 , figpos(2)-wdt );
   figpos([2 4]) = figpos([2 4]) + [ -1  1 ] * dev;
   set(fig,'position',figpos);
end

%------------------------------------------------------

axe = axes( 'parent' , fig , ...
            'units'  , 'normalized' , ...
            'xlim'   , xl  , ...
            'ylim'   , yl  , ...
            'xdir'   , 'normal'  , ...
            'ydir'   , 'reverse' , ...
            'xgrid'  , 'on'  , ...
            'ygrid'  , 'on'  , ...
    'gridlinestyle'  , '-'  , ...
            'xtick'  , xt  , ...
            'ytick'  , yt  , ...
       'xticklabel'  , []  , ...
       'yticklabel'  , []  , ...
             'box'   , 'on'       , ...
           'color'   , 'none'     , ...
        'linestyle'  , '-'        , ...
        'linewidth'  , 0.1        , ...
            'layer'  , 'top'      , ...
         'nextplot'  , 'add'      , ...
              'tag'  , [app '_Axes']       );

ht = get(axe,'title');

set( ht , 'string'     , 'Ferret-Colors' , ...
          'fontname'   , FontName   , ...
          'fontunits'  , FontUnits  , ...
          'fontsize'   , FontSize+2 , ...
          'fontweight' , 'bold'     , ... 
          'fontangle'  , FontAngle  , ...
          'visible'    , 'on'       , ...
          'tag'        , [app '_Title']      );

%******************************************************
% Determine AxesPosition

h = text( 'parent'     , axe        , ...
          'units'      , 'inches'   , ...
          'fontname'   , FontName   , ...
          'fontunits'  , FontUnits  , ...
          'fontsize'   , FontSize   , ...
          'fontweight' , FontWeight , ... 
          'fontangle'  , FontAngle , ...
          'visible'    , 'off'           );

ext = zeros(2,4);

for ii = [ 1  2 ]
    str = char(nfc(ii:2:n,1));
    set( h , 'string' , str );
    ext(ii,:) = get( h , 'extent' );
    ext(ii,1) = ext(ii,3) / size(str,2);
    ext(ii,3) = ext(ii,3) * ( 1 + 2/size(str,2) );  % 2 Characters more
    ext(ii,4) = ext(ii,4) / size(str,1);            % Single Character
end

delete(h);

x0 = max(ext(:,1));   % XOffset for Text to Axes

axepos = zeros(1,4);

axepos(1) = ext(1,3);
axepos(3) = pappos(3) - ext(1,3) - ext(2,3);
axepos(2) = 1.5 * ext(1,4);
axepos(4) = pappos(4) - 3.5 * ext(1,4);

x0 = x0 / axepos(3) * ( xl(2) - xl(1) );  % XOffset in Units 'data'

pos = axepos ./ pappos([3 4 3 4]);

set( axe , 'position' , pos );

%******************************************************
% Show ColorMaps

hcc = uicontextmenu( 'parent' , fig , ...
                   'userdata' , nfc(:,1) , ...
                        'tag' , [app '_ColorContext'] );

ini = { ...
  'CNR'    0  '&Color   %.0f/%.0f  of  "%s" %s'
  'RGB'    1  ' &RGB: %4.2f %4.2f %4.2f | %3.3d %3.3d %3.3d | #%s '
  'HSV'    0  ' &HSV: %4.2f %4.2f %4.2f | %3.3d %3.3d %3.3d | #%s '
  'GRAY'   1  ' &Gray:  %4.2f  |  %3.3d  |  #%s '
  'EDIT'   1  '&Inquire ...' };

sets = { 'off'  'on' };

for ii = 1 : size(ini,1)

    hmc = uimenu( 'parent' , hcc , ...
                  'label'  , ''  , ...
                'userdata' , ini{ii,3} , ...
               'separator' , sets{1+ini{ii,2}} , ... 
                     'tag' , [ app '_' ini{ii,1} ]  );

   if ii == size(ini,1)

      set( hmc , 'label' , ini{ii,3} , ...
              'userdata' , [ 0  0  0 ] , ...
              'callback' , sprintf('%s(gcbo,''Inquire'');',fcn) );

   end

end

ButtonDown = sprintf('%s(gcbo,''ImageDown'');',fcn);

hi = image( 'parent' , axe , ...
            'xdata'  , xd  , ...
            'ydata'  , yd(:)  , ...
            'cdata'  , cc  , ...
         'userdata'  , cc  , ...      % Original
     'uicontextmenu' , hcc , ...
     'buttondownfcn' , ButtonDown , ...
               'tag' , [ app '_Image' ]  ); 

%******************************************************
% Labels with uicontextmenu

CallBack = sprintf('%s(gcbo,''Assign'');',fcn);

hcc = uicontextmenu( 'parent' , fig , ...
                        'tag' , [ app '_LabelContext' ]);

hmc = uimenu( 'parent' , hcc , ...
               'label' , '&Assign' , ...
            'callback' , CallBack , ...
                'tag' , [ app '_LabelContextMenu' ]);

ht = zeros(n,1);

hh = { 'right'  'left' };

ButtonDown = sprintf('%s(gcbo,''LabelDown'');',fcn);

for ii = 1 : n

    jj = mod( ii-1 , 2 ) + 1;

    ff = 2 * jj - 3;          % Sign for XOffset

    % Set CellString, taking care before 'default'

   ht(ii) = text( 'parent'    , axe        , ...
                  'string'    , nfc(ii,1) , ...
                 'units'      , 'data'     , ...
                 'position'   , [ xl(jj)+ff*x0 ii 0 ] , ...
        'horizontalalignment' , hh{jj}     , ...
          'verticalalignment' , 'middle'   , ...
                 'fontname'   , FontName   , ...
                 'fontunits'  , FontUnits  , ...
                 'fontsize'   , FontSize   , ...
                 'fontweight' , FontWeight , ... 
                 'fontangle'  , FontAngle  , ...
                 'clipping'   , 'off'      , ...
                 'visible'    , 'on'       , ...
                'interpreter' , 'none'     , ...
             'buttondownfcn'  , ButtonDown , ...
             'uicontextmenu'  , hcc        , ...
                'userdata'    , { ii  get(hi,'tag') }  , ...
                'tag'         , [app '_Label']         );
end

%******************************************************
% Slider

ButtonDown = sprintf('%s(gcbo,''SliderDown'');',fcn);

CallBack = sprintf('%s(gcbo,''SliderCallback'');',fcn);

ToolTip = [ 'Click Left to Rotate Colors, '  ...
            'Click Right to Invert Colors.'  char(10) ...
            'Click Right on Label to Assign ColorMap. ' ...
            'Click Right on Image to see ColorValue.'       ];

pos([2 4]) = axepos(2) * [ 1/3  2/3 ] ./ pappos([4 4]);

hc = uicontrol('parent' , fig      , ...
               'units'  , 'normalized' , ...
             'position' , pos      , ...
               'style'  , 'slider' , ...
               'min'    , 0 , ...
               'max'    , 1 , ...
               'value'  , 0 , ...
          'sliderstep'  , [ 5/360 30/360 ] , ...
           'userdata'   , []          , ...        % Set below !!!
           'callback'   , CallBack    , ...
      'buttondownfcn'   , ButtonDown  , ...
      'tooltipstring'   , ToolTip     , ...
        'interruptible' , 'off'       , ...
                'tag'  , [app '_Slider']  );

%------------------------------------------------------
% Label for SliderValue

    jj = mod( (n+1)-1 , 2 ) + 1;

    ff = 2 * jj - 3;          % Sign for XOffset

    frm = [' %.0f' char(176)  ' '];

    p2 = yl(2) + diff(yl) * ( axepos(2)/pappos(4) - pos(2) );

       hl = text( 'parent'    , axe        , ...
                  'string'    , ''         , ...
                 'units'      , 'data'     , ...
                 'position'   , [ xl(jj)+ff*x0 p2 0 ] , ...
        'horizontalalignment' , hh{jj}     , ...
          'verticalalignment' , 'baseline' , ...
                 'fontname'   , FontName   , ...
                 'fontunits'  , FontUnits  , ...
                 'fontsize'   , FontSize   , ...
                 'fontweight' , 'bold'     , ... 
                 'fontangle'  , FontAngle  , ...
                 'clipping'   , 'off'      , ...
                 'visible'    , 'on'       , ...
                'interpreter' , 'none'     , ...
             'buttondownfcn'  , ''         , ...
             'uicontextmenu'  , []         , ...
                'userdata'    , frm        , ...
                'tag'         , [app '_SliderLabel']       );


set( hc , 'userdata' , { get(hi,'tag')  get(hl,'tag')  is_grey  0 } );

%------------------------------------------------------------------
% Initialisation of Label

cb = get(hc,'callback');

cb = strrep(cb,'gcbo',epsstr(hc));

eval(cb);   % Initialisation of Label

%------------------------------------------------------------------
% FigureUserdata: XRGB-ColorNames

  [c,g] = rgb_colors;

  cg = struct( 'c' , { struct( 'c' , { cat(1,c{:,2}) } , ...
                               'm' , { size(c,1) }     , ...
                               'n' , { c(:,1) }            ) } , ...
               'g' , { struct( 'c' , { cat(1,g{:,2}) } , ...
                               'm' , { size(g,1) }     , ...
                               'n' , { g(:,1) }            ) } );


  cg.c.c = cg.c.c / 255;
  cg.g.c = cg.g.c(:,1) / 255;

  set(fig,'userdata',cg);

%******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function fig = ferret_cb(h,action,app)

ud = get(h,'userdata');

switch upper(action)

%***********************************************
case 'SLIDERDOWN'
%***********************************************

f = get(h,'parent');

hi = findobj(f,'type','image','tag',ud{1});

set( hi , 'cdata' , 1 - get(hi,'cdata')   , ...
       'userdata' , 1 - get(hi,'userdata') );

ud{4} = ~ud{4};

set( h , 'userdata' , ud );

%***********************************************
case 'SLIDERCALLBACK'
%***********************************************

f = get(h,'parent');

v = round( get(h,'value') * 360 ) / 360;

hi = findobj(f,'type','image','tag',ud{1});
hl = findobj(f,'type','text' ,'tag',ud{2});

cc = get(hi,'userdata');

cc = rgb2hsv(get(hi,'userdata'));

cc(:,:,1) = cc(:,:,1) + v;
cc(:,:,1) = cc(:,:,1) - floor(cc(:,:,1));

cc = min(max(hsv2rgb(cc),0),1);

if ud{3}
   cc = 0.299*cc(:,:,1) + 0.587*cc(:,:,2) + 0.114*cc(:,:,3);
   cc = cc(:,:,[1 1 1]);
end

set( hi , 'cdata'  , cc);
set( hl , 'string' , sprintf(get(hl,'userdata'),360*v) );
set( h  , 'value' , v );

%***********************************************
case 'LABELDOWN'
%***********************************************

c = get(h,'uicontextmenu');
m = get(c,'children');       % AssignMenu

f = get(c,'parent');

n =  char(get(h,'string'));  % Name

set( c , 'userdata' , n );

set( m , 'userdata' , ud , ...
            'label' , sprintf('&Assign "%s" %s',n,gri(f,app)) );

%***********************************************
case 'ASSIGN'
%***********************************************

c = get(h,'parent');

f = get(c,'parent');

hi = findobj(f,'type','image','tag',ud{2});

cc = get( hi , 'cdata' );

cc = cc(ud{1},:,:);

cc = permute(cc,[2 3 1]);

assignin('base',get(c,'userdata'),cc);

%***********************************************
case 'IMAGEDOWN'
%***********************************************

cc = get( h , 'cdata' );

si = size(cc(:,:,1));

axe = get( h ,'parent');

cp = get( axe , 'currentpoint' );

cp = min( floor(cp(1,[2 1])-0.5)+1 , si );

cc = permute(cc(cp(1),cp(2),:),[1 3 2]);

c = get(h,'uicontextmenu');

%---------------------------------------------
% Get ColorName

rad = 0.03;

f = get(c,'parent');     % Figure

cg = get(f,'userdata');  % ColorNames

ig = all( abs(cc-mean(cc)) <= 1/255 );

name = '';

if ~ig
    d = ( cg.c.c - cc(ones(1,cg.c.m),:) );
    d = sum( d.^2 , 2 );
    ok = ( d <= rad^2 );
    if any(ok)
       nk =  sum(ok);
       ok = find(ok);
       if nk > 1
          [d,ii] = min(d(ok));
          ok = ok(ii);
       end
       name = cg.c.n{ok};
    end
end

if isempty(name)
   if ig
      d = ( cg.g.c - cc(1) );
   else
      d = ( cg.g.c(:,[1 1 1]) - cc(ones(1,cg.g.m),:) );
   end
    d = sum( d.^2 , 2 );
    ok = ( d <= rad^2 );
    if any(ok)
       nk =  sum(ok);
       ok = find(ok);
       if nk > 1
          [d,ii] = min(d(ok));
          ok = ok(ii);
       end
       name = cg.g.n{ok};
    end
end

%---------------------------------------------

n = get(c,'userdata');   % Names

m = get(c,'children');   % SubMenu's

t = strrep( get(m,'tag') , [app '_'] , '' );

%  'CNR'    0  '&Color   %.0f/%.0f  of  "%s" '
%  'RGB'    1  ' &RGB: %4.2f %4.2f %4.2f / %3.3d %3.3d %3.3d / #%s '
%  'HSV'    0  ' &HSV: %4.2f %4.2f %4.2f / %3.3d %3.3d %3.3d / #%s '
%  'GRAY'   1  ' &Gray:  %4.2f  /  %3.3d  /  #%s '
%  'EDIT'   1  '&Inquire ...' };

ii = find(strcmp(t,'CNR'));
if ~isempty(ii)
    frm = get(m(ii),'userdata');
    set( m(ii) , 'label' , sprintf(frm,cp(2),si(2),n{cp(1)},gri(f,app)) );
end

%-----------------------------------------------------------
% Set ColorLabels etc.

for id = { 'RGB' 'HSV' 'GRAY' 'EDIT' }

    ii = strcmp(id{1},t);

    if sum(ii) == 1

        ii = find(ii);
        
        if strcmp(id,'HSV')
           ci = rgb2hsv(cc);
        elseif strcmp(id,'GRAY')
           ci = 0.299*cc(1) + 0.587*cc(2) + 0.114*cc(3);
        elseif strcmp(id,'EDIT')
           frm = '&Inquire %s...';
           if ~isempty(name)
               name = sprintf('"%s" ',name);
           end
           lab = sprintf(frm,name);
           set( m(ii) , 'label' , lab , 'userdata' , cc );
        else
           ci = cc;
        end

        if ~strcmp(id,'EDIT')

            cj = round( 255 * ci );

           frm = get(m(ii),'userdata');

      if strcmp(id,'GRAY')
           ch = sprintf(' %s',dec2hex(cj,2));
      else
           ch = sprintf( ' %s %s %s' , dec2hex(cj(1),2) , ...
                                       dec2hex(cj(2),2) , ...
                                       dec2hex(cj(3),2)         );
      end

           lab = sprintf(frm,ci,cj,ch);

           set( m(ii) , 'label' , lab );

        end

    end

end

%***********************************************
case 'INQUIRE'
%***********************************************

c = get(h,'parent');

m = findobj(c,'type','uimenu','tag',[app '_CNR']);

if isempty(m)
   t = 'Inquire';
else
   t = get(m,'label');
   t = strrep(t,'&','');
end

uisetcolor(get(h,'userdata'),t);

%***********************************************

end

%******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = gri(fig,app);

% Get Grey, Rotation, Invert

sld = findobj(fig,'type','uicontrol','style','slider','tag',[app '_Slider']);

rot = '';
inv = '';
gry = '';

frm = [' %.0f' char(176)  ' '];

if ~isempty(sld)
    rot = sprintf(frm,get(sld,'value')*360);
    sud = get(sld,'userdata');
    if sud{4}
       inv = ' inverted ';
    end
    if sud{3}
        gry = ' grey ';
    end
end

str = sprintf('%s%s%s',gry,rot,inv);

%******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [ h , form ] = epsstr( h );

% EPSSTR  Transform Number exact into String,
%          
% using Matlab's floating point relative accuracy.
%
%  Form = EPSSTR;   
%    returns Format for using with SPRINTF
%
%  [ String , Form ] = EPSSTR( Number ); 
%    returns exact String for Number
%
%  Form = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 )
%


form = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 );

if nargin < 1

  h = form;

  return

end


if ~isnumeric(h)  |  ( prod(size(h)) > 1 )
 error('Handle must be a single Numeric.');
end


  h = sprintf(form,h);

%******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function wysiwyg(fig)

% WYSIWYG  WhatYouSeeIsWhatYouGet
%
% WYGIWYS( FigureHandle )
%
%WYSIWYG -- changes the size of the figure on the screen to equal
%       the size of the figure that would be printed, 
%       according to the papersize attribute.  Use this function
%       to give a more accurate picture of what will be 
%       printed.
%       Dan(K) Braithwaite, Dept. of Hydrology U.of.A  11/93
 
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

unis = get(fig,'units');
ppos = get(fig,'paperposition');

org = get(fig,'position');

set(fig,'units',get(fig,'paperunits'));

pos = get(fig,'position');
pos(3:4) = ppos(3:4);

set(fig,'position',pos);
set(fig,'units',unis);

new = get(fig,'position');

if ~all( new([1 2]) == org([1 2]) )
    new([1 2]) = org([1 2]);
    set(fig,'position',new);
end


%******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [f,p] = whichfile(varargin)

% WHICHFILE   List Files in Matlab's SearchPath
%
% [ FileName , DirName ] = WHICHFILE( FileEnd1 , FileEnd2 , ... )
%
% List Files in Matlab's SearchPath with the Syntax: *FileEnd#
%
% FileEnd# are Strings, FileName and DirName are CellArray's of Strings.
%

Nout = nargout;

f = cell(0,1);
p = cell(0,1);

if nargin == 0
   return
end

[ok,v] = chkcstr(varargin);

if ~ok
   error('Inputs must be Strings');
end

nv = prod(size(v));

pfad = sepname(matlabpath,NaN,pathsep);

fs = filesep;

for pp = pfad

    for ii = 1 : nv

        d = dir(cat(2,pp{1},fs,'*',v{ii}));

        if ~isempty(d)
            ok = find(~cat(1,d.isdir));
            if ~isempty(ok)
                f = cat( 1 , f , {d(ok).name}' );
                p = cat( 1 , p , pp(ones(length(ok),1)) );
            end
        end
 
    end

end   

if Nout == 0

   s = { filesep };
   s = s(ones(size(f)));

   s = strhcat(permute(cat(2,p,s,f),[2 1]),'',3);

   fprintf(1,'\n%s\n\n',s);

   clear f

end

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = sepname(name,dep,sep,mode);

% SEPNAME  separates Name
%
% String = SEPNAME( Name , Depth , Seperator ,Mode )
%
% Depth gives the number of Recursions
% 
% Use Depth = NaN  to get all parts of Name in String.
%  In this case String is a CellStringArray.
%
% Mode =  1  Recursion from Start
%        -1  Recursion from End
%         0  Recursion from End, 
%             special Handling for Class- and Private-Directories
%
% Defaults:  Depth     = 0
%            Seperator = FILESEP
%            Mode      = 0    ( 1  if Depth = NaN )
%

Nin = nargin;

if Nin < 1
   name = '';
end

if Nin < 2
   dep = 0;
end

dep_nan = isnan(dep);

if Nin < 3
   sep = filesep;
end

if Nin < 4
   mode = dep_nan;
end

%********************************************

if dep_nan
   str = cell(1,0);
else
   str = char(zeros(1,0));
end

if isempty(name)
   return
end

if ~chkstr(name)
   error('Name must be a String.')
end

if ~( chkstr(sep,1) & ( prod(size(sep)) == 1 ) )
   error('Seperator must be a single Character.')
end

n = size(name,2);

%---------------------------------------------
% Find Seperator in Name

is = ( double(name) == double(sep) );

if all(is)
   if dep_nan
      str    = cell(1,n-1);
      str(:) = { '' };
   end      
   return
end

%---------------------------------------------

i0 = ~is(1);
i1 = ~is(n);

is = cat( 2 , ones(1,i0) , is , ones(1,i1) );

is = find( is );

is = is(:);

ni = size(is,1) - 1;

if ni == 0 
   return
end
     
%---------------------------------------------
% [ Start  End ]

ind = ( 1 : ni ) ;

is  = cat( 2 , is(ind)+1 , is(ind+1)-1 ) - i0;

%---------------------------------------------
% Take care for duplicate Seperators

if ~dep_nan | ( double(sep) == char(32) )

   is( find( is(:,1) > is(:,2) ) , : ) = [];

end

%---------------------------------------------

ni = size(is,1);

if ni == 0
   return
end

%---------------------------------------------
if dep_nan
   
   ind = [ 1  ni ];

  flip = ~( mode == 1 );

   ind = ind( [ 1  2 ] + [ 1 -1 ]*flip );

   ind = ( ind(1) : 1-2*flip : ind(2) );

   is = is(ind,:);

   str = cell(1,ni);
  
   for ii = 1 : ni

       str{ii} = name( is(ii,1) : is(ii,2) );

   end

   return

end

%---------------------------------------------

ii = ni - 1 * ( ni > 1 );

nn = name( is(ii,1) : is(ii,2) );

ic = strcmp( nn(1) , '@'       );
ip = strcmp( nn    , 'private' );

id = 1 * ic + 2 * ip;

dep = dep + 1 + id * ( mode == 0 ) * ( ni > 1 );

dep = min( dep , ni );

%---------------------------------------------

is(1,1) = is(1,1) - 1 * ( is(1,1) > 1 );  % Start incl. Seperator

ind = ( 1 : ni-1 );

is(ind,2) = is(ind,2) + 1;                % End incl. Seperator

is(:,2) = is(:,2) - is(:,1) + 1;          % Length

%---------------------------------------------

ind = ( 1 : dep ) + ( ni - dep ) * (~( mode == 1 ));

is  = is(ind,:);

ind = grp2ind(is(:,1),is(:,2));

str = name(ind);

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ii = grp2ind(i0,l);

% GRP2IND  Built IndexVector from StartIndex and Length
%
% Index = GRP2IND( StartIndex , GroupLength )
%

if isempty(i0);
   ii = [];
   return
end

si = size(i0);

if ( sum( si > 1 ) > 1 )
   error('StartIndex must be a Vector.');
end

i0 = i0(:);
l  =  l(:);

if ~isequal(size(i0,1),size(l,1))
   error('Size of StartIndex and GroupLenght must be the same.');
end

n = size(l,1);

ii = ones(sum(l),1);
jj = cumsum( cat(1,1,l) , 1 );

ii(jj(1:n)) = i0;

if n > 1
   ii(jj(2:n)) = ii(jj(2:n))-(i0(1:n-1,1)+l(1:n-1)-1);
end

ii = cumsum(ii,1);

if n == 1
   return
end

jj = find( si > 1 );

if jj == 1
   return
end
  
perm = cat( 2 , (1:jj-1)+1 , 1 , (jj+1:size(si,2)) );

ii = permute(ii,perm);

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

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

function  str = strhcat(str,del,n,nl)

% STRHCAT  Concatenates Strings into ONE
%
% STRHCAT( StringArray , Delimiter )
%   Forms one long String from the Strings in the
%   StringArray, delimited with the delimiter.
%   The EndDelimiter will be removed.
%
% STRHCAT( StringArray , Delimiter , N , NewLine )
%   Build a  NewLine after each N-th String.
%   default: N = 10;  NewLine = char(10);
%
% Example:  
%         >> strhcat({'apples' 'pies' 'prunes'},', ')
%    
%         ans =
%
%         apples, pies, prunes
%
%         >> strhcat({'apples';'pies';'prunes'},', ',2)
%    
%         ans =
%
%         apples, pies
%         prunes
%



Nin = nargin;

if Nin < 4
 nl = char(10);
end
if Nin < 3
 n = [];
end
if Nin < 2
 del = char((32*ones(1,3)));
end


if isempty(str)
 str = '';
 return
end


if ischar(str)
  str = cellstr(str);
end

str = str(:);

if isempty(n)
   n = size(str,1) + 1;
end

str( : , 2 ) = { del };

str(n:n:size(str,1),2) = { nl };

str(    size(str,1),2) = { '' };


str = permute( str , [ 2  1 ] );

str = cat(2,str{:});


%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [c,g] = rgb_colors

c = { ...

  'alice blue'               [ 240  248  255 ]
  'antique white'            [ 250  235  215 ]
  'AntiqueWhite1'            [ 255  239  219 ]
  'AntiqueWhite2'            [ 238  223  204 ]
  'AntiqueWhite3'            [ 205  192  176 ]
  'AntiqueWhite4'            [ 139  131  120 ]
  'aquamarine1'              [ 127  255  212 ]
  'aquamarine2'              [ 118  238  198 ]
  'aquamarine3'              [ 102  205  170 ]
  'aquamarine4'              [ 069  139  116 ]
  'azure1'                   [ 240  255  255 ]
  'azure2'                   [ 224  238  238 ]
  'azure3'                   [ 193  205  205 ]
  'azure4'                   [ 131  139  139 ]
  'beige'                    [ 245  245  220 ]
  'bisque1'                  [ 255  228  196 ]
  'bisque2'                  [ 238  213  183 ]
  'bisque3'                  [ 205  183  158 ]
  'bisque4'                  [ 139  125  107 ]
  'blanched almond'          [ 255  235  205 ]
  'blue1'                    [ 000  000  255 ]
  'blue2'                    [ 000  000  238 ]
  'blue3'                    [ 000  000  205 ]
  'blue4'                    [ 000  000  139 ]
  'blue violet'              [ 138  043  226 ]
  'brown'                    [ 165  042  042 ]
  'brown1'                   [ 255  064  064 ]
  'brown2'                   [ 238  059  059 ]
  'brown3'                   [ 205  051  051 ]
  'brown4'                   [ 139  035  035 ]
  'burlywood'                [ 222  184  135 ]
  'burlywood1'               [ 255  211  155 ]
  'burlywood2'               [ 238  197  145 ]
  'burlywood3'               [ 205  170  125 ]
  'burlywood4'               [ 139  115  085 ]
  'cadet blue'               [ 095  158  160 ]
  'CadetBlue1'               [ 152  245  255 ]
  'CadetBlue2'               [ 142  229  238 ]
  'CadetBlue3'               [ 122  197  205 ]
  'CadetBlue4'               [ 083  134  139 ]
  'chartreuse1'              [ 127  255  000 ]
  'chartreuse2'              [ 118  238  000 ]
  'chartreuse3'              [ 102  205  000 ]
  'chartreuse4'              [ 069  139  000 ]
  'chocolate'                [ 210  105  030 ]
  'chocolate1'               [ 255  127  036 ]
  'chocolate2'               [ 238  118  033 ]
  'chocolate3'               [ 205  102  029 ]
  'chocolate4'               [ 139  069  019 ]
  'coral'                    [ 255  127  080 ]
  'coral1'                   [ 255  114  086 ]
  'coral2'                   [ 238  106  080 ]
  'coral3'                   [ 205  091  069 ]
  'coral4'                   [ 139  062  047 ]
  'cornflower blue'          [ 100  149  237 ]
  'cornsilk1'                [ 255  248  220 ]
  'cornsilk2'                [ 238  232  205 ]
  'cornsilk3'                [ 205  200  177 ]
  'cornsilk4'                [ 139  136  120 ]
  'cyan1'                    [ 000  255  255 ]
  'cyan2'                    [ 000  238  238 ]
  'cyan3'                    [ 000  205  205 ]
  'cyan4'                    [ 000  139  139 ]
  'dark goldenrod'           [ 184  134  011 ]
  'DarkGoldenrod1'           [ 255  185  015 ]
  'DarkGoldenrod2'           [ 238  173  014 ]
  'DarkGoldenrod3'           [ 205  149  012 ]
  'DarkGoldenrod4'           [ 139  101  008 ]
  'dark green'               [ 000  100  000 ]
  'dark khaki'               [ 189  183  107 ]
  'dark olive green'         [ 085  107  047 ]
  'DarkOliveGreen1'          [ 202  255  112 ]
  'DarkOliveGreen2'          [ 188  238  104 ]
  'DarkOliveGreen3'          [ 162  205  090 ]
  'DarkOliveGreen4'          [ 110  139  061 ]
  'dark orange'              [ 255  140  000 ]
  'DarkOrange1'              [ 255  127  000 ]
  'DarkOrange2'              [ 238  118  000 ]
  'DarkOrange3'              [ 205  102  000 ]
  'DarkOrange4'              [ 139  069  000 ]
  'dark orchid'              [ 153  050  204 ]
  'DarkOrchid1'              [ 191  062  255 ]
  'DarkOrchid2'              [ 178  058  238 ]
  'DarkOrchid3'              [ 154  050  205 ]
  'DarkOrchid4'              [ 104  034  139 ]
  'dark salmon'              [ 233  150  122 ]
  'dark sea green'           [ 143  188  143 ]
  'DarkSeaGreen1'            [ 193  255  193 ]
  'DarkSeaGreen2'            [ 180  238  180 ]
  'DarkSeaGreen3'            [ 155  205  155 ]
  'DarkSeaGreen4'            [ 105  139  105 ]
  'dark slate blue'          [ 072  061  139 ]
  'dark slate gray'          [ 047  079  079 ]
  'DarkSlateGray1'           [ 151  255  255 ]
  'DarkSlateGray2'           [ 141  238  238 ]
  'DarkSlateGray3'           [ 121  205  205 ]
  'DarkSlateGray4'           [ 082  139  139 ]
  'dark turquoise'           [ 000  206  209 ]
  'dark violet'              [ 148  000  211 ]
  'deep pink'                [ 255  020  147 ]
  'DeepPink2'                [ 238  018  137 ]
  'DeepPink3'                [ 205  016  118 ]
  'DeepPink4'                [ 139  010  080 ]
  'deep sky blue'            [ 000  191  255 ]
  'DeepSkyBlue2'             [ 000  178  238 ]
  'DeepSkyBlue3'             [ 000  154  205 ]
  'DeepSkyBlue4'             [ 000  104  139 ]
  'dodger blue'              [ 030  144  255 ]
  'DodgerBlue2'              [ 028  134  238 ]
  'DodgerBlue3'              [ 024  116  205 ]
  'DodgerBlue4'              [ 016  078  139 ]
  'firebrick'                [ 178  034  034 ]
  'firebrick1'               [ 255  048  048 ]
  'firebrick2'               [ 238  044  044 ]
  'firebrick3'               [ 205  038  038 ]
  'firebrick4'               [ 139  026  026 ]
  'floral white'             [ 255  250  240 ]
  'forest green'             [ 034  139  034 ]
  'ghost white'              [ 248  248  255 ]
  'gold1'                    [ 255  215  000 ]
  'gold2'                    [ 238  201  000 ]
  'gold3'                    [ 205  173  000 ]
  'gold4'                    [ 139  117  000 ]
  'goldenrod'                [ 218  165  032 ]
  'goldenrod1'               [ 255  193  037 ]
  'goldenrod2'               [ 238  180  034 ]
  'goldenrod3'               [ 205  155  029 ]
  'goldenrod4'               [ 139  105  020 ]
  'green1'                   [ 000  255  000 ]
  'green2'                   [ 000  238  000 ]
  'green3'                   [ 000  205  000 ]
  'green4'                   [ 000  139  000 ]
  'green yellow'             [ 173  255  047 ]
  'honeydew1'                [ 240  255  240 ]
  'honeydew2'                [ 224  238  224 ]
  'honeydew3'                [ 193  205  193 ]
  'honeydew4'                [ 131  139  131 ]
  'hot pink'                 [ 255  105  180 ]
  'HotPink1'                 [ 255  110  180 ]
  'HotPink2'                 [ 238  106  167 ]
  'HotPink3'                 [ 205  096  144 ]
  'HotPink4'                 [ 139  058  098 ]
  'indian red'               [ 205  092  092 ]
  'IndianRed1'               [ 255  106  106 ]
  'IndianRed2'               [ 238  099  099 ]
  'IndianRed3'               [ 205  085  085 ]
  'IndianRed4'               [ 139  058  058 ]
  'ivory1'                   [ 255  255  240 ]
  'ivory2'                   [ 238  238  224 ]
  'ivory3'                   [ 205  205  193 ]
  'ivory4'                   [ 139  139  131 ]
  'khaki'                    [ 240  230  140 ]
  'khaki1'                   [ 255  246  143 ]
  'khaki2'                   [ 238  230  133 ]
  'khaki3'                   [ 205  198  115 ]
  'khaki4'                   [ 139  134  078 ]
  'lavender'                 [ 230  230  250 ]
  'lavender blush'           [ 255  240  245 ]
  'LavenderBlush2'           [ 238  224  229 ]
  'LavenderBlush3'           [ 205  193  197 ]
  'LavenderBlush4'           [ 139  131  134 ]
  'lawn green'               [ 124  252  000 ]
  'lemon chiffon'            [ 255  250  205 ]
  'LemonChiffon2'            [ 238  233  191 ]
  'LemonChiffon3'            [ 205  201  165 ]
  'LemonChiffon4'            [ 139  137  112 ]
  'light blue'               [ 173  216  230 ]
  'LightBlue1'               [ 191  239  255 ]
  'LightBlue2'               [ 178  223  238 ]
  'LightBlue3'               [ 154  192  205 ]
  'LightBlue4'               [ 104  131  139 ]
  'light coral'              [ 240  128  128 ]
  'light cyan'               [ 224  255  255 ]
  'LightCyan2'               [ 209  238  238 ]
  'LightCyan3'               [ 180  205  205 ]
  'LightCyan4'               [ 122  139  139 ]
  'light goldenrod'          [ 238  221  130 ]
  'LightGoldenrod1'          [ 255  236  139 ]
  'LightGoldenrod2'          [ 238  220  130 ]
  'LightGoldenrod3'          [ 205  190  112 ]
  'LightGoldenrod4'          [ 139  129  076 ]
  'light goldenrod yellow'   [ 250  250  210 ]
  'light green'              [ 144  238  144 ]
  'light pink'               [ 255  182  193 ]
  'LightPink1'               [ 255  174  185 ]
  'LightPink2'               [ 238  162  173 ]
  'LightPink3'               [ 205  140  149 ]
  'LightPink4'               [ 139  095  101 ]
  'light salmon'             [ 255  160  122 ]
  'LightSalmon2'             [ 238  149  114 ]
  'LightSalmon3'             [ 205  129  098 ]
  'LightSalmon4'             [ 139  087  066 ]
  'light sea green'          [ 032  178  170 ]
  'light sky blue'           [ 135  206  250 ]
  'LightSkyBlue1'            [ 176  226  255 ]
  'LightSkyBlue2'            [ 164  211  238 ]
  'LightSkyBlue3'            [ 141  182  205 ]
  'LightSkyBlue4'            [ 096  123  139 ]
  'light slate blue'         [ 132  112  255 ]
  'light slate gray'         [ 119  136  153 ]
  'light steel blue'         [ 176  196  222 ]
  'LightSteelBlue1'          [ 202  225  255 ]
  'LightSteelBlue2'          [ 188  210  238 ]
  'LightSteelBlue3'          [ 162  181  205 ]
  'LightSteelBlue4'          [ 110  123  139 ]
  'light yellow'             [ 255  255  224 ]
  'LightYellow2'             [ 238  238  209 ]
  'LightYellow3'             [ 205  205  180 ]
  'LightYellow4'             [ 139  139  122 ]
  'lime green'               [ 050  205  050 ]
  'linen'                    [ 250  240  230 ]
  'magenta1'                 [ 255  000  255 ]
  'magenta2'                 [ 238  000  238 ]
  'magenta3'                 [ 205  000  205 ]
  'magenta4'                 [ 139  000  139 ]
  'maroon'                   [ 176  048  096 ]
  'maroon1'                  [ 255  052  179 ]
  'maroon2'                  [ 238  048  167 ]
  'maroon3'                  [ 205  041  144 ]
  'maroon4'                  [ 139  028  098 ]
  'medium orchid'            [ 186  085  211 ]
  'MediumOrchid1'            [ 224  102  255 ]
  'MediumOrchid2'            [ 209  095  238 ]
  'MediumOrchid3'            [ 180  082  205 ]
  'MediumOrchid4'            [ 122  055  139 ]
  'medium purple'            [ 147  112  219 ]
  'MediumPurple1'            [ 171  130  255 ]
  'MediumPurple2'            [ 159  121  238 ]
  'MediumPurple3'            [ 137  104  205 ]
  'MediumPurple4'            [ 093  071  139 ]
  'medium sea green'         [ 060  179  113 ]
  'medium slate blue'        [ 123  104  238 ]
  'medium spring green'      [ 000  250  154 ]
  'medium turquoise'         [ 072  209  204 ]
  'medium violet red'        [ 199  021  133 ]
  'midnight blue'            [ 025  025  112 ]
  'mint cream'               [ 245  255  250 ]
  'misty rose'               [ 255  228  225 ]
  'MistyRose2'               [ 238  213  210 ]
  'MistyRose3'               [ 205  183  181 ]
  'MistyRose4'               [ 139  125  123 ]
  'moccasin'                 [ 255  228  181 ]
  'navajo white'             [ 255  222  173 ]
  'NavajoWhite2'             [ 238  207  161 ]
  'NavajoWhite3'             [ 205  179  139 ]
  'NavajoWhite4'             [ 139  121  094 ]
  'navy'                     [ 000  000  128 ]
  'old lace'                 [ 253  245  230 ]
  'olive drab'               [ 107  142  035 ]
  'OliveDrab1'               [ 192  255  062 ]
  'OliveDrab2'               [ 179  238  058 ]
  'OliveDrab4'               [ 105  139  034 ]
  'orange1'                  [ 255  165  000 ]
  'orange2'                  [ 238  154  000 ]
  'orange3'                  [ 205  133  000 ]
  'orange4'                  [ 139  090  000 ]
  'orange red'               [ 255  069  000 ]
  'OrangeRed2'               [ 238  064  000 ]
  'OrangeRed3'               [ 205  055  000 ]
  'OrangeRed4'               [ 139  037  000 ]
  'orchid'                   [ 218  112  214 ]
  'orchid1'                  [ 255  131  250 ]
  'orchid2'                  [ 238  122  233 ]
  'orchid3'                  [ 205  105  201 ]
  'orchid4'                  [ 139  071  137 ]
  'pale goldenrod'           [ 238  232  170 ]
  'pale green'               [ 152  251  152 ]
  'PaleGreen1'               [ 154  255  154 ]
  'PaleGreen3'               [ 124  205  124 ]
  'PaleGreen4'               [ 084  139  084 ]
  'pale turquoise'           [ 175  238  238 ]
  'PaleTurquoise1'           [ 187  255  255 ]
  'PaleTurquoise2'           [ 174  238  238 ]
  'PaleTurquoise3'           [ 150  205  205 ]
  'PaleTurquoise4'           [ 102  139  139 ]
  'pale violet red'          [ 219  112  147 ]
  'PaleVioletRed1'           [ 255  130  171 ]
  'PaleVioletRed2'           [ 238  121  159 ]
  'PaleVioletRed3'           [ 205  104  137 ]
  'PaleVioletRed4'           [ 139  071  093 ]
  'papaya whip'              [ 255  239  213 ]
  'peach puff'               [ 255  218  185 ]
  'PeachPuff2'               [ 238  203  173 ]
  'PeachPuff3'               [ 205  175  149 ]
  'PeachPuff4'               [ 139  119  101 ]
  'peru'                     [ 205  133  063 ]
  'pink'                     [ 255  192  203 ]
  'pink1'                    [ 255  181  197 ]
  'pink2'                    [ 238  169  184 ]
  'pink3'                    [ 205  145  158 ]
  'pink4'                    [ 139  099  108 ]
  'plum'                     [ 221  160  221 ]
  'plum1'                    [ 255  187  255 ]
  'plum2'                    [ 238  174  238 ]
  'plum3'                    [ 205  150  205 ]
  'plum4'                    [ 139  102  139 ]
  'powder blue'              [ 176  224  230 ]
  'purple'                   [ 160  032  240 ]
  'purple1'                  [ 155  048  255 ]
  'purple2'                  [ 145  044  238 ]
  'purple3'                  [ 125  038  205 ]
  'purple4'                  [ 085  026  139 ]
  'red1'                     [ 255  000  000 ]
  'red2'                     [ 238  000  000 ]
  'red3'                     [ 205  000  000 ]
  'red4'                     [ 139  000  000 ]
  'rosy brown'               [ 188  143  143 ]
  'RosyBrown1'               [ 255  193  193 ]
  'RosyBrown2'               [ 238  180  180 ]
  'RosyBrown3'               [ 205  155  155 ]
  'RosyBrown4'               [ 139  105  105 ]
  'royal blue'               [ 065  105  225 ]
  'RoyalBlue1'               [ 072  118  255 ]
  'RoyalBlue2'               [ 067  110  238 ]
  'RoyalBlue3'               [ 058  095  205 ]
  'RoyalBlue4'               [ 039  064  139 ]
  'salmon'                   [ 250  128  114 ]
  'salmon1'                  [ 255  140  105 ]
  'salmon2'                  [ 238  130  098 ]
  'salmon3'                  [ 205  112  084 ]
  'salmon4'                  [ 139  076  057 ]
  'sandy brown'              [ 244  164  096 ]
  'sea green'                [ 046  139  087 ]
  'SeaGreen1'                [ 084  255  159 ]
  'SeaGreen2'                [ 078  238  148 ]
  'SeaGreen3'                [ 067  205  128 ]
  'seashell1'                [ 255  245  238 ]
  'seashell2'                [ 238  229  222 ]
  'seashell3'                [ 205  197  191 ]
  'seashell4'                [ 139  134  130 ]
  'sienna'                   [ 160  082  045 ]
  'sienna1'                  [ 255  130  071 ]
  'sienna2'                  [ 238  121  066 ]
  'sienna3'                  [ 205  104  057 ]
  'sienna4'                  [ 139  071  038 ]
  'sky blue'                 [ 135  206  235 ]
  'SkyBlue1'                 [ 135  206  255 ]
  'SkyBlue2'                 [ 126  192  238 ]
  'SkyBlue3'                 [ 108  166  205 ]
  'SkyBlue4'                 [ 074  112  139 ]
  'slate blue'               [ 106  090  205 ]
  'SlateBlue1'               [ 131  111  255 ]
  'SlateBlue2'               [ 122  103  238 ]
  'SlateBlue3'               [ 105  089  205 ]
  'SlateBlue4'               [ 071  060  139 ]
  'slate gray'               [ 112  128  144 ]
  'SlateGray1'               [ 198  226  255 ]
  'SlateGray2'               [ 185  211  238 ]
  'SlateGray3'               [ 159  182  205 ]
  'SlateGray4'               [ 108  123  139 ]
  'snow1'                    [ 255  250  250 ]
  'snow2'                    [ 238  233  233 ]
  'snow3'                    [ 205  201  201 ]
  'snow4'                    [ 139  137  137 ]
  'spring green'             [ 000  255  127 ]
  'SpringGreen2'             [ 000  238  118 ]
  'SpringGreen3'             [ 000  205  102 ]
  'SpringGreen4'             [ 000  139  069 ]
  'steel blue'               [ 070  130  180 ]
  'SteelBlue1'               [ 099  184  255 ]
  'SteelBlue2'               [ 092  172  238 ]
  'SteelBlue3'               [ 079  148  205 ]
  'SteelBlue4'               [ 054  100  139 ]
  'tan'                      [ 210  180  140 ]
  'tan1'                     [ 255  165  079 ]
  'tan2'                     [ 238  154  073 ]
  'tan4'                     [ 139  090  043 ]
  'thistle'                  [ 216  191  216 ]
  'thistle1'                 [ 255  225  255 ]
  'thistle2'                 [ 238  210  238 ]
  'thistle3'                 [ 205  181  205 ]
  'thistle4'                 [ 139  123  139 ]
  'tomato1'                  [ 255  099  071 ]
  'tomato2'                  [ 238  092  066 ]
  'tomato3'                  [ 205  079  057 ]
  'tomato4'                  [ 139  054  038 ]
  'turquoise'                [ 064  224  208 ]
  'turquoise1'               [ 000  245  255 ]
  'turquoise2'               [ 000  229  238 ]
  'turquoise3'               [ 000  197  205 ]
  'turquoise4'               [ 000  134  139 ]
  'violet'                   [ 238  130  238 ]
  'violet red'               [ 208  032  144 ]
  'VioletRed1'               [ 255  062  150 ]
  'VioletRed2'               [ 238  058  140 ]
  'VioletRed3'               [ 205  050  120 ]
  'VioletRed4'               [ 139  034  082 ]
  'wheat'                    [ 245  222  179 ]
  'wheat1'                   [ 255  231  186 ]
  'wheat2'                   [ 238  216  174 ]
  'wheat3'                   [ 205  186  150 ]
  'wheat4'                   [ 139  126  102 ]
  'yellow1'                  [ 255  255  000 ]
  'yellow2'                  [ 238  238  000 ]
  'yellow3'                  [ 205  205  000 ]
  'yellow4'                  [ 139  139  000 ]
  'yellow green'             [ 154  205  050 ]
  'black'                    [ 000  000  000 ]

};


g = { ...

  'gray1'                    [ 003  003  003 ]
  'gray2'                    [ 005  005  005 ]
  'gray3'                    [ 008  008  008 ]
  'gray4'                    [ 010  010  010 ]
  'gray5'                    [ 013  013  013 ]
  'gray6'                    [ 015  015  015 ]
  'gray7'                    [ 018  018  018 ]
  'gray8'                    [ 020  020  020 ]
  'gray9'                    [ 023  023  023 ]
  'gray10'                   [ 026  026  026 ]
  'gray11'                   [ 028  028  028 ]
  'gray12'                   [ 031  031  031 ]
  'gray13'                   [ 033  033  033 ]
  'gray14'                   [ 036  036  036 ]
  'gray15'                   [ 038  038  038 ]
  'gray16'                   [ 041  041  041 ]
  'gray17'                   [ 043  043  043 ]
  'gray18'                   [ 046  046  046 ]
  'gray19'                   [ 048  048  048 ]
  'gray20'                   [ 051  051  051 ]
  'gray21'                   [ 054  054  054 ]
  'gray22'                   [ 056  056  056 ]
  'gray23'                   [ 059  059  059 ]
  'gray24'                   [ 061  061  061 ]
  'gray25'                   [ 064  064  064 ]
  'gray26'                   [ 066  066  066 ]
  'gray27'                   [ 069  069  069 ]
  'gray28'                   [ 071  071  071 ]
  'gray29'                   [ 074  074  074 ]
  'gray30'                   [ 077  077  077 ]
  'gray31'                   [ 079  079  079 ]
  'gray32'                   [ 082  082  082 ]
  'gray33'                   [ 084  084  084 ]
  'gray34'                   [ 087  087  087 ]
  'gray35'                   [ 089  089  089 ]
  'gray36'                   [ 092  092  092 ]
  'gray37'                   [ 094  094  094 ]
  'gray38'                   [ 097  097  097 ]
  'gray39'                   [ 099  099  099 ]
  'gray40'                   [ 102  102  102 ]
  'gray41'                   [ 105  105  105 ]
  'gray42'                   [ 107  107  107 ]
  'gray43'                   [ 110  110  110 ]
  'gray44'                   [ 112  112  112 ]
  'gray45'                   [ 115  115  115 ]
  'gray46'                   [ 117  117  117 ]
  'gray47'                   [ 120  120  120 ]
  'gray48'                   [ 122  122  122 ]
  'gray49'                   [ 125  125  125 ]
  'gray50'                   [ 127  127  127 ]
  'gray51'                   [ 130  130  130 ]
  'gray52'                   [ 133  133  133 ]
  'gray53'                   [ 135  135  135 ]
  'gray54'                   [ 138  138  138 ]
  'gray55'                   [ 140  140  140 ]
  'gray56'                   [ 143  143  143 ]
  'gray57'                   [ 145  145  145 ]
  'gray58'                   [ 148  148  148 ]
  'gray59'                   [ 150  150  150 ]
  'gray60'                   [ 153  153  153 ]
  'gray61'                   [ 156  156  156 ]
  'gray62'                   [ 158  158  158 ]
  'gray63'                   [ 161  161  161 ]
  'gray64'                   [ 163  163  163 ]
  'gray65'                   [ 166  166  166 ]
  'gray66'                   [ 168  168  168 ]
  'dark gray'                [ 169  169  169 ]
  'gray67'                   [ 171  171  171 ]
  'gray68'                   [ 173  173  173 ]
  'gray69'                   [ 176  176  176 ]
  'gray70'                   [ 179  179  179 ]
  'gray71'                   [ 181  181  181 ]
  'gray72'                   [ 184  184  184 ]
  'gray73'                   [ 186  186  186 ]
  'gray74'                   [ 189  189  189 ]
  'gray'                     [ 190  190  190 ]
  'gray75'                   [ 191  191  191 ]
  'gray76'                   [ 194  194  194 ]
  'gray77'                   [ 196  196  196 ]
  'gray78'                   [ 199  199  199 ]
  'gray79'                   [ 201  201  201 ]
  'gray80'                   [ 204  204  204 ]
  'gray81'                   [ 207  207  207 ]
  'gray82'                   [ 209  209  209 ]
  'light gray'               [ 211  211  211 ]
  'gray83'                   [ 212  212  212 ]
  'gray84'                   [ 214  214  214 ]
  'gray85'                   [ 217  217  217 ]
  'gray86'                   [ 219  219  219 ]
  'gainsboro'                [ 220  220  220 ]
  'gray87'                   [ 222  222  222 ]
  'gray88'                   [ 224  224  224 ]
  'gray89'                   [ 227  227  227 ]
  'gray90'                   [ 229  229  229 ]
  'gray91'                   [ 232  232  232 ]
  'gray92'                   [ 235  235  235 ]
  'gray93'                   [ 237  237  237 ]
  'gray94'                   [ 240  240  240 ]
  'gray95'                   [ 242  242  242 ]
  'gray96'                   [ 245  245  245 ]
  'gray97'                   [ 247  247  247 ]
  'gray98'                   [ 250  250  250 ]
  'gray99'                   [ 252  252  252 ]
  'white'                    [ 255  255  255 ]

};
