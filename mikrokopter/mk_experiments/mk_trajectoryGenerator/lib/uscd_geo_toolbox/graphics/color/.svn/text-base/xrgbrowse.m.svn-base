function fig = xrgbrowse(fig,action,varargin)

% XRGBROWSE  Browser for XRGB-Colors
%
% Creates a GUI to explore the Colors by XRGB (397)
%
% see also: XRGB (required), COLSPEC, FERRETC
%

Nin  = nargin;
Nout = nargout;

msg = '';

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
app = upper(fcn);   % ApplicationData of Figure and ColorAxes

%**********************************************************
% Check Inputs
%----------------------------------------------------------
if Nin == 0
%----------------------------------------------------------

   fig = [];

   try
      fig = newfig(fcn,app);
   catch
      msg = lasterr;
   end

   if ischar(fig)
      msg = fig;
      fig = [];
   end

   if ~isempty(msg)
       msg = sprintf('Can''t build XRGB-GUI.\n%s',msg);
   end

%----------------------------------------------------------
elseif Nin < 2
%----------------------------------------------------------

   msg = 'Not enough InputArguments';

%----------------------------------------------------------
else
%----------------------------------------------------------

   hh = fig;

   ok = ( isnumeric(fig) & ( prod(size(fig)) == 1 ) );
   if ok
      ok = ishandle(fig);
      if ok
         typ = get(fig,'type');
         ok  = strcmp(typ,'figure');
         if ~ok
             ok = any(strcmp(typ,{'uicontrol' 'uicontextmenu' 'uimenu'}));
             if ok
                switch typ
                  case {'uicontrol' 'uicontextmenu'}
                     fig = get(fig,'parent');
                  case  'uimenu'
                     while ~strcmp(typ,'figure')
                        fig = get(fig,'parent');
                        typ = get(fig,'typ');
                     end
                end
             end
         end
         if ok  
            ok = isappdata(fig,app);
         end
      end
   end

   if ~ok
       msg = 'First Input must be a FigureHandle or UI-Handle.';
   end

   if ~( ischar(action) & ~isempty(action) & ...
        ( prod(size(action)) == size(action,2) )   )
       m = 'Action must be a nonempty String.';
       if isempty(msg)
          msg = m;
       else
          msg = sprintf('%s\n%s',msg,m);
       end
   end

%----------------------------------------------------------
end
%----------------------------------------------------------

if ~isempty(msg)
    if Nout == 0
          error(msg);
    end
    fig = msg;
    return
end

if Nin == 0
   if Nout == 0
      clear fig
   end
   return
end

%**********************************************************

apd = getappdata(fig,app);

switch upper(action)

%----------------------------------------------------------
case 'BUTTONDOWN'
%----------------------------------------------------------

      nr = varargin{1};

      tag_frame(apd.Frame,'Activate',nr,1);

%----------------------------------------------------------
case 'ACTIVATE'
%----------------------------------------------------------

     nr = varargin{2};

     ht = apd.Label;
     hc = apd.Color;

     ini = apd.ini;
    name = apd.name;

     rgb = apd.rgb;
     hsv = apd.hsv;
     gry = apd.gry;

    urgb = apd.urgb;
    uhsv = apd.uhsv;
    ugry = apd.ugry;

    hrgb = apd.hrgb;
    hhsv = apd.hhsv;
    hgry = apd.hgry;

     fld = fieldnames(ini);

     set(ht,'visible','off');
     set(hc,'visible','off');
     
     if nr == 0
        set( apd.Intro , 'visible' , 'on' );
        return
     end

     set( apd.Intro , 'visible' , 'off' );

     %-----------------------------------------------------------------
     % Set Labels and Colors

     frm = [ '\n %s\n\n' ...
             '  RGB: %4.2f %4.2f %4.2f | %3.3d %3.3d %3.3d | # %s %s %s \n' ...
             '  HSV: %4.2f %4.2f %4.2f | %3.3d %3.3d %3.3d | # %s %s %s \n' ...
             '  Gray:  %4.2f  |  %3.3d  |  # %s \n\n Inquire by click Right ...\n' ];

     cc = getfield( ini , fld{nr} );

     for ii = 1 : size(cc,2)

         set( ht(ii) , 'string' , [ cc{1,ii} ' ' ], ...
                   'fontweight' , 'bold'   , ...
                   'fontangle'  , 'normal' , ...
                      'visible' , 'on' );

         if strcmp(cc{1,ii},'single')
            set( ht(ii) , 'string' , 'additional ' , ...
                      'fontweight' , 'normal' , ...
                      'fontangle'  , 'italic'        );
         end

         ind = cc{3,ii};

         for jj = 1 : prod(size(ind))

             kk = ind(jj);

             tip = sprintf( frm , name{kk} , ...
                            rgb(kk,:) , urgb(kk,:) , hrgb{kk,:} , ...
                            hsv(kk,:) , uhsv(kk,:) , hhsv{kk,:} , ...
                            gry(kk)   , ugry(kk)   , hgry{kk}         );
                            
             set( hc(ii,jj) , 'backgroundcolor' , rgb(kk,:) , ...
                                'tooltipstring' , tip       , ...
                                     'userdata' , name{kk}  , ...
                                      'visible' , 'on' );

         end

     end

end

%**********************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function fig = newfig(fcn,app);

%**********************************************************

ini = struct( ...
...
'red'     , { { 'pale violet red' 'violet red' 'orange red' 'indian red' 'red' 'fire brick' ...
                'tomato'  'coral' 'light salmon' 'salmon' 'rosy brown' 'peach puff' 'misty rose'  } } , ...
...
'brown'   , { { 'peru' 'tan' 'chocolate' 'brown' 'burly wood' 'navajo white'                ...
                'bisque' 'wheat' 'rosy brown'                                                     } } , ...
...
'gold'    , { { 'light golden rod ' 'dark golden rod' 'golden rod' 'gold'                   ...
                'sienna' 'orange' 'dark orange' 'orange red' 'peach puff'                         } } , ...
...
'lemon'   , { { 'lemon chiffon' 'light yellow' 'khaki' 'yellow'                                   } } , ...
...
'green'   , { { 'olive drab' 'dark olive' 'chartreuse' 'green' 'spring green'               ...
                'pale green' 'aquamarine' 'sea green' 'dark sea green'                            } } , ...
...
'cyan'    , { { 'aquamarine' 'cadet blue' 'dark slate gray' 'pale turquoise' 'turquoise'    ...
                'cyan' 'light cyan' 'azure'                                                       } } , ...
...
'blue'    , { { 'cadet blue' 'light sky blue' 'sky blue' 'deep sky blue'                    ...
                'dark slate gray' 'slate gray' 'light blue' 'light steel blue' 'steel blue' ...
                'dodger blue' 'royal blue' 'blue' 'slate blue'                                    } } , ...
...
'violet'  , { { 'medium purple' 'purple' 'dark orchid' 'medium orchid' 'orchid'             ...
                'pale violet red' 'violet red' 'violet'                                           } } , ...
...
'pink'    , { { 'magenta' 'maroon' 'deep pink' 'hot pink' 'light pink' 'pink' 'plum'        ...
                'thistle' 'lavender' 'rosy brown'                                                 } } , ...
...
'gray'    , { { 'white' 'snow' 'seashell' 'cornsilk' 'ivory' 'honey dew'                    ...
                'light cyan' 'azure' 'misty rose' 'antique white' 'lemon chiffon ' 'light yellow' } } , ...
...
'single'  , { { 'linen' 'blanched almond' 'papaya whip' 'moccasin' 'old lace' 'beige' 'mint cream' } } );

col = struct( ...
'red'     , { [ 1   0   0   ] } , ...
'brown'   , { [ 0.8 0.4 0.1 ] } , ...
'gold'    , { [ 1   0.7 0.1 ] } , ...
'lemon'   , { [ 1   1   0   ] } , ...
'green'   , { [ 0   1   0   ] } , ...
'cyan'    , { [ 0   1   1   ] } , ...
'blue'    , { [ 0   0   1   ] } , ...
'violet'  , { [ 0.5 0   1   ] } , ...
'pink'    , { [ 1   0   1   ] } , ...
'gray'    , { [ 0.5 0.5 0.5 ] } , ...
'single'  , { [ 1   1   1   ] }       );

%**********************************************************

su = get(0,'units');      set(0,'units','pixels')
si = get(0,'ScreenSize'); set(0,'units', su );

fig = figure( 'units'            , 'pixels'    , ...
              'name'             , 'XRGB-Browser' , ...
              'numbertitle'      , 'off'       , ...
              'menubar'          , 'none'      , ...
              'toolbar'          , 'none'      , ...
              'color'            , 'w'         , ...
              'visible'          , 'off'       , ...
              'tag'              ,  app              );

%-----------------------------------------------------------------

is_win = strcmp(computer,'PCWIN');

h0 = uicontrol( 'parent' , fig , ...
                 'units' , 'pixels' , ...
                 'style' , 'text'   , ...
            'FontWeight' , 'bold'   , ...
       'ForeGroundColor' , 'k'      , ...
       'BackGroundColor' , 'w'             );

hh = get(h0,'position');
hh = hh(4);

hi = ( 1 + 0.5*is_win ) * hh;

pos = [ NaN NaN NaN 5*hi ];

pos(3) = ceil( 1.5 * pos(4) );

pos(1) = si(3)/2 - pos(3)/2;
pos(2) = si(4) - pos(4) - 60;

set( h0 , 'position' , [ hi 3*hi pos(3)-2*hi hi ] , ...
          'string'   , 'Initialize XRGB-Colors'   , ...
          'horizontalalignment' , 'center' );

h1 = uicontrol( 'parent' , fig , ...
                 'units' , 'pixels' , ...
              'position' , [ hi 1*hi pos(3)-2*hi hi ] , ...
                 'style' , 'text'   , ...
              'string'   , ''       , ...
            'FontWeight' , 'bold'   , ...
       'ForeGroundColor' , 'k'      , ...
       'BackGroundColor' , 'w'      , ...
   'horizontalalignment' , 'center'       );

if is_win
   set( [ h0 h1 ] , 'fontsize' , get(h1,'fontsize')-2 );
end

set( fig , 'units' , 'pixels' , 'position' , pos , 'visible' , 'on' );

%-----------------------------------------------------------------

[ini,nn,name,rgb] = init(ini,h1);

fld = fieldnames(ini);

vv = ini.single;
nv = cat(2,vv{2,:});
if all( nv == 1 )
   nv = sum(nv);
   vv = { 'single' ; nv ; cat(1,vv{3,:}) };
   ini.gray = cat( 2 , ini.gray , vv );
   ii = find( strcmp( fld , 'gray' ) );
   nn(ii,:) = nn(ii,:) + [ 1  max(0,nv-nn(ii,2)) nv ];
   ii = find( strcmp( fld , 'single' ) );
   fld(ii)  = [];
   nn(ii,:) = [];
   ini = rmfield(ini,'single');
end

%-----------------------------------------------------------------

set( h1 , 'string' , 'Transform Colors' );

hsv = rgb2hsv(rgb);
gry = 0.299*rgb(:,1) + 0.587*rgb(:,2) + 0.114*rgb(:,3);

urgb = round( 255 * rgb );
uhsv = round( 255 * hsv );
ugry = round( 255 * gry );

hrgb = cat( 2 , cellstr(dec2hex(urgb(:,1),2)) , ...
                cellstr(dec2hex(urgb(:,2),2)) , ...
                cellstr(dec2hex(urgb(:,3),2))       );

hhsv = cat( 2 , cellstr(dec2hex(uhsv(:,1),2)) , ...
                cellstr(dec2hex(uhsv(:,2),2)) , ...
                cellstr(dec2hex(uhsv(:,3),2))       );

hgry = cellstr(dec2hex(ugry,2));

pause(0.1);

%-----------------------------------------------------------------

set(h1,'string','Build GUI');

nb = size(fld,1);

fs = 12 - 2 * is_win;

[msg,hf,hb] = tag_frame( fig , 'new'         , ...
                      'TagNumber' ,  nb      , ...
                     'Fontsize'   ,  fs      , ...
                     'FontWeight' , 'bold'   , ...
                'ForeGroundColor' , 'k'      , ...
                'BackGroundColor' , 'w'      , ...
                          'CBFcn' , { fcn fig } , ...
                  'Interruptible' , 'off'    , ...
                     'BusyAction' , 'cancel'       );

if ~isempty(msg)
    msg = sprintf('Error call TAG_FRAME.\n%s',msg);
    try, delete(fig), end
    fig = msg;
    return
end

set( hf , 'deletefcn' , '' );

tag_frame(hf,'visible','off');

drawnow

frm = '%s(gcbo,''ButtonDown'',%.0f);';

btd = sprintf(frm,fcn,0);

for ii = 1 : nb
    cb = sprintf(frm,fcn,ii);
    cl = getfield(col,fld{ii});
    cl = rgb2hsv(cl);
    cl(2) = 0.5 * cl(2);
    cl(3) = cl(3) + ( 1 - 0.5 ) * ( 1 - cl(3) );
    cl = hsv2rgb(cl);
    set( hb(ii) , 'string' , fld{ii} , ...
                'callback' , cb      , ...
           'buttondownfcn' , btd     , ...
         'backgroundcolor' , cl            );
end

%-----------------------------------------------------------------

fud = get(hf,'userdata');
fud = fud.TAG_FRAME;

hh = max( hh , fud.PixelHight );

dh = ceil( hh / 2 );

nn = max(nn,[],1);

nt = 0;
for ii = 1 : nb
    vv = getfield(ini,fld{ii});
    nt = max( nt , size(char(vv(1,:)),2) );
end

wt = nt * hh * 1/2;           % LabelWidth

vp = ( 1 : nn(1) ) - 1;
vp = vp * ( hh + dh ) + hh;   % VertPos

hp = ( 1 : nn(2) ) - 1;
hp = hp * hh + 2*dh + wt;     % HorizPos of ColorFields

ht = zeros(nn(1),1);          % Label
hc = zeros(nn(1),nn(2));      % ColorFields

cb = 'uisetcolor(get(gcbo,''backgroundcolor''),get(gcbo,''userdata''));';
 
for ii = 1 : nn(1)

    [m,ht(ii)] = tag_frame( hf , 'add' , 1 , 'uicontrol'   , ...
                        'units' , 'pixels'             , ...
                     'position' , [ dh -vp(ii)-0.2*hh wt hh*0.8 ] , ... 
                        'style' , 'text'   , ...
                       'string' ,  ''      , ...
                   'Fontsize'   ,  fs      , ...
                   'FontWeight' , 'bold'   , ...
              'ForeGroundColor' , 'k'      , ...
              'BackGroundColor' , 'w'      , ...
          'horizontalalignment' , 'right'        );

     for jj = 1 : nn(2)

        [m, hc(ii,jj)] = tag_frame( hf , 'add' , 1 , 'uicontrol'   , ...
                        'units' , 'pixels'             , ...
                     'position' , [ hp(jj) -vp(ii) hh hh ] , ... 
                        'style' , 'frame'   , ...
              'ForeGroundColor' , 'k'       , ...
              'BackGroundColor' , 'w'       , ...
                     'UserData' , ''        , ... 
                'ButtonDownFcn' , cb     );


     end

end

%-----------------------------------------------------------------
% Resize

fp = get(hf,'position');

np = [ hp(nn(2))  vp(nn(1)) ] + hh + dh + 2*fud.BorderWidth;

np(1) = max( np(1) , 3 * hh * nb );

dp = np - fp([3 4]);

pos    = get( fig , 'position' );
pos(1) = pos(1) - dp(1)/2;
pos(3) = pos(3) + dp(1);
pos(2) = pos(2) - dp(2);
pos(4) = pos(4) + dp(2);

delete([h0 h1]);

set( fig , 'position' , pos , 'resize' , 'off' );

tag_frame(hf,'resize','absolut');

%-----------------------------------------------------------------
% Colored Background

brd = 3;                  % PixelBorder

p =  get(hf,'position');

p = p([3 4]) - 2*brd;

ok = ( ( gry > 0.4 ) & ( gry < 0.8 ) & ( hsv(:,2) > 0.6 ) );

nc = sum(ok);

ok = find(ok);

cc = hsv(ok,:);

cc(:,2) = 0.1 * cc(:,2);
cc(:,3) = cc(:,3) + ( 1 - 0.1 ) * ( 1 - cc(:,3) );

cc = hsv2rgb(cc);

nn = sqrt( nc * ( (p(1)/p(2)) .^ [ 1  -1 ] ) );  % [ Nx  Ny ]

nn = ceil(nn);

ix = min( nn(1) , floor(linspace(1,nn(1)+1,p(1))) );
iy = min( nn(2) , floor(linspace(1,nn(2)+1,p(2))) )';

ic = iy(:,ones(1,p(1))) + ( ix(ones(1,p(2)),:) - 1 ) * nn(2);

nn = prod(nn);

[h,si] = sort(rand(nc,1));

nn = ceil( nc * rand(nn-nc,1) );
nn = min(max(nn,1),nc);

si = cat( 1 , si , nn );

ic = si(ic);

cc = cc(cat(3,ic,ic+nc,ic+2*nc));


hb = uicontrol(  'parent' , fig , ...
                  'units' , 'pixels' , ...
               'position' ,  get(hf,'position')       , ...
                  'style' , 'pushbutton'   , ...
                 'string' ,  ''      , ...
                 'cdata'  ,  cc      , ...
        'ForeGroundColor' , 'k'      , ...
        'BackGroundColor' , 'w'      , ...
    'horizontalalignment' , 'center'        );

%-----------------------------------------------------------------
% Intro

intro = { ''
  sprintf('Explore %.0f Colors by XRGB',size(rgb,1));
          ''
          'Select by ColorType' 
          '' 
          'Move the Mouse over ColorFields to get the Color'
          ''  };

cpr = { [ char(169) ' 2005  Christian Begler' ]
        ''
        'IFM-GEOMAR Kiel'
        ''
        'cbegler@ifm-geomar.de' };

ni = prod(size(intro));
nc = prod(size(cpr));

p = get(hf,'position');

p([1 3]) = p([1 3]) + hh * [ 1 -2 ];
p(2)     = floor(pos(4)/2) - hh * ( ni + nc + 1 ) / 2 + hh * ( nc + 1 );
p(4)     = hh * ni;

hi = zeros(1,2);

hi(1) = uicontrol(  'parent' , fig , ...
                     'units' , 'pixels' , ...
                  'position' ,  p       , ...
                     'style' , 'text'   , ...
                    'string' ,  intro   , ...
                'Fontsize'   ,  fs+2    , ...
                'FontWeight' , 'bold'   , ...
           'ForeGroundColor' , 'k'      , ...
           'BackGroundColor' , 'w'      , ...
       'horizontalalignment' , 'center'        );

p(2) = p(2) - hh * ( nc + 1 );
p(4) = hh * nc;

hi(2) = uicontrol(  'parent' , fig , ...
                     'units' , 'pixels' , ...
                  'position' ,  p       , ...
                     'style' , 'text'   , ...
                    'string' ,  cpr     , ...
                'Fontsize'   ,  fs      , ...
                'FontWeight' , 'bold'   , ...
           'ForeGroundColor' , 'k'      , ...
           'BackGroundColor' , 'w'      , ...
       'horizontalalignment' , 'center'        );


tag_frame(hf,'visible','on');

%------------------------------------------------------------------------------

p = get(hf,'position');

a = axes( 'units' , 'pixels' , 'position' , p , 'visible' , 'off' );

figure(fig), drawnow

c = getframe(a); delete(a);

p = p([3 4]) - 2*brd;

c = c.cdata;

sc = size(c);

if all( sc([2 1]) >= p([1 2])+brd )

   c = c( brd+(1:p(2)) , brd+(1:p(1)) , : );
   c = double(c) / 255;

   cc = cc .* ~( c == 0 );

   set( hb , 'cdata' , cc );

   delete(hi);

   hi = hb;

else
    
   delete(hb);
   
end

%------------------------------------------------------------------------------

% Activate not by TAG_FRAME
set( fud.HideHandle(1) , 'userdata' , get(fud.HideHandle(2),'userdata') );

%------------------------------------------------------------------------------

apd = struct( 'Frame'  , { hf } , ...
              'Button' , { hb } , ...
              'Intro'  , { hi } , ...
              'Label'  , { ht } , ...
              'Color'  , { hc } , ...
              'ini'    , { ini  } , ...
              'name'   , { name } , ...
              'rgb'    , {  rgb } , ...
              'hsv'    , {  hsv } , ...
              'gry'    , {  gry } , ...
             'urgb'    , { urgb } , ...
             'uhsv'    , { uhsv } , ...
             'ugry'    , { ugry } , ...
             'hrgb'    , { hrgb } , ...
             'hhsv'    , { hhsv } , ...
             'hgry'    , { hgry }       );

setappdata(fig,app,apd);

%**********************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ini,nn,n,c] = init(ini,ht)


set( ht , 'string' , 'Call XRGB' );

[n,c] = xrgb(-1); pause(0.1);

m = strrep(lower(n),' ','');

[m,si] = sort(m);

n = n(si);
c = c(si,:);

m = sprintf('%s\n',m{:});

in = find( m == 10 );

nc = size(n,1);

fi = fieldnames(ini);

ni = size(fi,1);

nn =  zeros(ni,3);  % [ nv  max  sum ]

ok = zeros(nc,1);

for ii = 1 : ni

    set( ht , 'string' , fi{ii} ), drawnow

    v = getfield(ini,fi{ii});

    w = strrep(v,' ','');

     bl     = sum( char(w) == 32 , 2 );
    [bl,si] = sort(bl);
    
    nv = size(v,2);

    v      = v(ones(1,3),:);
    v(2,:) = {00};             % Number
    v(3,:) = {[]};             % Index

    ff = ones(nc,1);

    ind{ii} = cell(1,nv);

    for jj = si(:)'

        kk = findstr(m,w{jj});

        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if strcmp(w{jj},'royalblue')
           kk = cat( 2 , kk , findstr(m,'navy') );
        end
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if ~isempty(kk)

            for ik = 1 : prod(size(kk))
                kk(ik) = sum( in < kk(ik) ) + 1;
                kk(ik) = kk(ik) * ff(kk(ik));
            end

            ik = ( kk == 0 );
            if     all(ik)
                   kk = [];
            elseif any(ik)
                   ik = find(~ik);
                   kk = kk(ik);
            end

        end

        if ~isempty(kk)

            ff(kk) = 0;

            v{2,jj} = prod(size(kk));

            v{3,jj} = kk;

        end

    end

    ini = setfield( ini , fi{ii} , v );

    vv = cat( 2 , v{2,:} );

    nn(ii,:) = [ nv max(vv) sum(vv) ];

    ok = ( ok | ~ff );

    pause(0.1);

end

if any(~ok)
   ii = find(~ok);
   w = sprintf('%s\n',n{ii});
   w = sprintf('Unmatched XRGB-Colors:\n%s',w);
   ww = warnstat;
   if ~strcmp(ww,'off')
       warning('on');
       warning(w);
       warning(ww);
   end
end

%**********************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [Msg,varargout] = tag_frame(h,action,varargin);

% TAG_FRAME   Creates a TagFrame
%
% [ Msg, ... ] = TAG_FRAME( Handle , Action , ... );
%
%  Note:  The first Output contains ErrorMessages for invalid Inputs.
%         If all is ok, Msg is empty.
%
%-------------------------------------------------------------------------
% 1. Create a new TagFrame in a Figure, Action 'NEW'
%
%   
% [ Msg, FrameHandle, ButtonHandle, PixelPosition, ButtonHight ] = ...
%    TAG_FRAME( FigureHandle, 'New', Property1, PropertyValue1, ... );
%
% Creates a new TagFrame in the Figure, defined by 
%   FigureHandle, and returns the Handle of the FrameBox, TabButtons, and 
%   the Position in Pixels. The Units for the Position of the ListBox 
%   in the Figure are set to 'pixels'.
% 
% Additional Inputs are Properties and their Values:
%
%  'TagNumber'   , Number of Tab's                 , default: 3
%  'Position'    , [  Left  Bottom  -Right -Top  ] , default: [ 3  3 -3 -3]
%                  [  Left  Bottom   Width  High ] , ... and all Combinations ...
%                  [ -Right -Top     Width  High ]
%  'BorderWidth' , PixelWidth of UIControl-Border  , default: 3
%  'HighPadding' , Padding of TabButton in Hight   , default: 1.5
%  'CBFcn'       , CallBackFcn, called if a Tag is activated:
%                  feval( CBFcn{:} , 'Activate' , 1 , Number )
%
% More Inputs are FontProperties and their Values, by default 
%  the DefaultUIControlFontProperties of the Root will used.
% 
%
% Please:  don't change the Property 'Tag' of the TagFrameHandle,
%           this Property will used to identify the Handle as valid 
%           TagFrame.
%
%
% The UserData of the TagFrameHandle contains all Properties in a StructArray.
%
%-------------------------------------------------------------------------
% 2. Add an Object, Action 'ADD'
%
%  [ Msg , Handles ] = TAG_FRAME( FrameHandle , 'Add' , nr , type , varargin )
%  
%   nr = [ 0 .. ud.TagNumber ]  |  FrameHandle  |  ButtonHandle
%
% type = 'msg_logo' |  'tag_frame'  |  'tab_list'  |  'sel_list'  |  'uicontrol'
%
% Following OutPuts: 
%
%   'msg_logo' :  MessageListHandle, LogoFrameHandle, LogoTextHandle
%   'tag_frame':  TagFrameHandle,  TagButtonHandle, PixelPosition, ButtonHight
%   'sel_list' :  SelFrameHandle,  TextHandle, PixelPosition
%
% UIControls can be positioned equal to the Positioning of the TagFrame, see under NEW
%
%-------------------------------------------------------------------------
% 3. Activate a Tag, Action 'ACTIVATE'
%
%  [ Msg , ActiveNumber ] = TAG_FRAME( FrameHandle , 'Activate' , nr )
%
%  nr = [ 0 .. ud.TagNumber ]  |  FrameHandle  |  ButtonHandle
%
%  nr = NaN sets TagFrame visible off
%
%-------------------------------------------------------------------------
% 4. Set Enability of Tag's, Action 'ENABLE'
%
%  Msg  = TAG_FRAME( FrameHandle , 'Enable' , sets , nr )
%
%  sets = 'off'  |  'on'
%
%   nr = [ 0 .. ud.TagNumber ]  |  FrameHandle  |  ButtonHandle
%  
%-------------------------------------------------------------------------
% 5. Set Visibility of TagFrame, Action 'VISIBLE'
%
%  Msg  = TAG_FRAME( FrameHandle , 'Visible' , sets )
%
%  sets = 'off'  |  'on'
%
%
%-----------------------------------------------------------------------------
% 6. Resize the TagFrame after the Figure was resized, Action 'RESIZE' 
%
%  [ Msg, PixelPosition, ButtonHight ] = TAG_FRAME( FrameHandle, 'Resize', Mode )
%
%  Resize the TagFrame in the Figure,  depending on Mode
%   call this from the ResizeFcnProperty of the Figure.  
%
%  Mode:  'normalized'  ( short:  'n' )
%
%   The TagFramePosition will set normalized to the Figure, with the normalized
%   Position when the TagFrame was created. 
%
%  Mode:  'absolut'     ( short:  'a' )
%
%   The TagFramePosition will set absolut, 
%    defined by the original Position.
%
%-----------------------------------------------------------------------------
% 7. Delete the TagFrame, Action 'DELETE'
%
%  MSG = TAG_FRAME( FrameHandle , 'Delete' )
%
%-----------------------------------------------------------------------------
%
% see also: MAKE_GUI, TAB_LIST, SEL_LIST, MSG_LIST, MSG_LOGO
%
%

Nout = nargout - 1;

Nout = Nout * ( Nout > 0 );

varargout = cell(1,Nout);



Msg = '';

nl = char(10);


Msg0 = 'TAG_FRAME: ';

Nin = nargin;

if Nin < 2
  Msg = [ Msg0  'Inputs H and Action are undefined.' ];
  return
end


ok = ( isnumeric(h) & ( prod(size(h)) == 1 ) );
if ok
 ok = ishandle(h);
end

if ~ok
 Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
         ' First Input must be a Handle.' ];
end


if ~ischar(action)
 Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
          'Action must be a String,'  ];
end

if ~isempty(Msg)
  Msg = [ Msg0   Msg ];
  return
end


action = upper(action);

Msg0 = [ 'TAG_FRAME( ' action ' ): ' ];

%------------------------------------------------------------
% Check Handle

if strcmp(action,'NEW')

 if ~strcmp(get(h,'type'),'figure')
   Msg = [ Msg0  'Handle must be a Figure.' ];
   return
 end

 fig = h;
   h = [];

%------------------------------------------------------------
else

 ok = strcmp(get(h,'type'),'uicontrol');
 if ok
   ok = ( strcmp(get(h,'style'),'frame')      &  ...
          strcmp(get(h,'tag'),'TAG_FRAME')          );
 end

 if ~ok
   Msg = [ Msg0  'Handle must be a ' ...
           'valid TagFrame-UIControl.' ];
   return
 end

 fig = get(h,'parent');

 ud0 = get(h,'userdata');

 ud  = ud0.TAG_FRAME;
  
end
%------------------------------------------------------------


  VarIn = varargin;
  VarIn = VarIn(:);
 

switch action

%*********************************************************
case 'NEW'

  
  if ( mod(size(VarIn,1),2) ~= 0 )
    Msg = [ 'Additional Inputs must contain UIControl-Property-Value-Pairs.' ];
  else 
    VarIn = reshape(VarIn,2,size(VarIn,1)/2)';
    if ~iscellstr(VarIn(:,1))
     Msg = [ 'Additional Inputs must contain UIControl-Property-Value-Pairs.' nl ...
             'Properties must be Strings.' ];
    end
  end

  if ~isempty(Msg)
     Msg = [ Msg0  Msg ];
     return
  end


  % FigurePosition in Pixels

  figuni = get(fig,'units');
           set(fig,'units','pixels');
  figpos = get(fig,'position');
           set(fig,'units',figuni);


  % DefaultUserData

  e1 = zeros(0,1);

  ch = struct( 'Handle'   , { e1        }  , ...
               'Position' , { cell(0,1) }        );

  ud = struct( 'TagNumber'     , {  3  }         , ... % Nr of Tab's
               'Position'      , {[ 3  3 -3 -3 ]} , ... % Position 
               'BorderWidth'   , {  3 }          , ... % BorderPixelWidth of UIControl
               'HighPadding'   , { 1.5 }         , ... % Position per ListboxLine
               'CBFcn'         , { {} }          , ... % CallBackFcn
               'PixelHight'    , {  0 }          , ... % PixelHight of Font
               'FrameHandle'   , { NaN }         , ... % Handle of Frame arround
               'ButtonHandle'  , { e1  }         , ... % Handle of Buttons
               'HideHandle'    , { e1  }         , ... % Handle of TextObject, using for width
               'ActiveNr'      , {  0  }         , ... % Nr of Active Button
               'FrameChildren' , { ch  }         , ... % Children of Frame
               'Visibility'    , { 'on' }        , ... % Visibility of TagFrame
               'PixelPosition' , { zeros(1,4) }  , ... % Actual PixelPosition [ Left Bottom Width High ]
               'NormPosition'  , { zeros(1,4) }  , ... % Default normalized Position
               'FigPos0'       , { figpos     }  );   % Default FigurePosition


  %-------------------------------------------------------------------
  % Look, if the first 5 Properties are given in Inputs

  fields = fieldnames(ud);
  fields = fields([1 2 3 4 5]);

  nv = size(VarIn,1);

  is_ud = zeros(nv,1);

  for ii = 1 : nv

     jj = find( strcmp( lower(VarIn{ii,1}) , lower(fields) ) );

     if ~isempty(jj)

       val = VarIn{ii,2};
       msg = '';

       if strcmp( fields{jj} , 'CBFcn' )

           if ischar(val)
              val = cellstr(val);
           end

           ok = iscell(val);
           if ok & ~isempty(val)
              ok = ( ischar(val{1}) & ~isempty(val{1}) & ...
                    ( prod(size(val{1})) == size(val{1},2) )  );
           end
 
           if ~ok
               msg = [  msg  nl(1:(end*(~isempty(msg)))) ...
                        'Value for '  fields{jj} ' must be a CharacterArray' ...
                        ' or CellArray with a String in the 1. Element.' ];

           end

       elseif ~isnumeric(val) | isempty(val)

            msg = [ 'Value for '  fields{jj} ' must be numeric and not empty.' ];

       elseif strcmp(fields{jj},'Position')

            ok =  all(isfinite(val));
            if ~ok
                msg = [ 'Value for '  fields{jj} ' must contain finite numerics.' ];
            end
            si0 = size(getfield(ud,fields{jj}));
            if ~isequal(size(val),si0)
               str = sprintf(' %.0f by',si0);
               str = [ '[' str(1:end-2) ']' ];
               msg = [  msg  nl(1:(end*(~isempty(msg)))) ...
                        'Value for '  fields{jj} ' must be of Size ' str '.' ];
            elseif ok
               if any( ( val([1 2]) < 0 )  &  ( val([3 4]) < 0 ) )

                  msg = [ 'Invalid Values for '  fields{jj} ...
                          ', type "help tag_frame" for more Informations.' ];

               end
            end

       else
            val_min = 1 - strcmp(fields{jj},'TagNumber');

            ok1 =  ( all( val >= val_min )  &  all(isfinite(val)) );
            ok2 =  ( all( mod(val,1) == 0 )  |  strcmp(fields{jj},'HighPadding') );
            if ~( ok1 & ok2 )
                msg = [ 'Value for '  fields{jj} ' must contain Integers >= ' ...
                         sprintf('%.0f',val_min) '.' ];
            end
            si0 = size(getfield(ud,fields{jj}));
            if ~isequal(size(val),si0)
               str = sprintf(' %.0f by',si0);
               str = [ '[' str(1:end-2) ']' ];
               msg = [  msg  nl(1:(end*(~isempty(msg)))) ...
                        'Value for '  fields{jj} ' must be of Size ' str '.' ];
            end
       end

        if isempty(msg)
          is_ud(ii) = 1;
             ud     = setfield(ud,fields{jj},val);
        else
             Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) msg ];
        end

     end
     % ~isempty(jj)
  end
  % ii

  if ~isempty(Msg)
    Msg = [ Msg0  Msg ];
    return
  end

  
  VarIn(find(is_ud),:) = [];

  %-------------------------------------------------------------------
  % Check UIControlProperties in the other Inputs

  VarIn(:,3) = { [] };  % DefaultValues


  % Set DefaultProperties 

  for ii = 1 : size(VarIn,1)

    msg = '';

    try

      VarIn{ii,3} = get(0,['DefaultUIControl' VarIn{ii,1}]);

    catch

       msg = lasterr;

    end

    if isempty(msg)
       set(0,['DefaultUIControl' VarIn{ii,1}],VarIn{ii,2});
    else
       Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
               'Invalid UIControlProperty: ' VarIn{ii,1} ];
       break
    end

  end

  if ~isempty(Msg)
     % Set DefaultProperties back
    for jj = 1 : ii-1
      set(0,['DefaultUIControl' VarIn{jj,1}],VarIn{jj,3});
    end

    Msg = [ Msg0  Msg ];
    return

  end

  %--------------------------------------------------------
  % Define DummyAxes, containing TextObject to get PixelHigh

  axe = axes('parent'  , fig , ...
                      'units'   , 'pixels' , ... 
                      'position',[ 0 0 1 1 ] , ...
                      'xlim'    , [ 0  1 ] , ...
                      'ylim'    , [ 0  1 ] , ...
                      'xtick'   , [] , ...
                      'ytick'   , [] , ...
                      'visible' , 'off' , ...
                      'nextplot' , 'add' , ... 
                      'handlevisibility' , 'callback'   );

  ht = text('parent' , axe , ...
          'units' , 'pixels' , ...
          'position',[ 1 1 0 ] , ...
          'string'  , 'H' , ... 
          'interpreter' , 'none' , ...
          'fontunits'  , get(0,'DefaultUIControlFontUnits') , ... 
          'fontsize'   , get(0,'DefaultUIControlFontSize')  , ... 
          'fontname'   , get(0,'DefaultUIControlFontName')  , ... 
          'fontangle'  , get(0,'DefaultUIControlFontAngle') , ... 
          'fontweight' , get(0,'DefaultUIControlFontWeight') , ...
          'visible'    , 'off'       );

   % Determine PixelHight

    e = get( ht , 'extent' );

   delete(ht);
   delete(axe);

   ud.PixelHight = e(4);

   %--------------------------------------------------------
 

   ud.ButtonHandle = NaN*ones(ud.TagNumber,1);
   ud.HideHandle   = NaN*ones(ud.TagNumber,1);


   for ii = 1 : ud.TagNumber

     ud.ButtonHandle(ii) = uicontrol('parent'   , fig     , ...
                                     'visible'  , 'off'   , ...
                                     'selected' , 'off'   , ...
                                     'units'    , 'pixels', ...
                                     'position' ,  [ 0  0  1  1 ] , ...
                                     'style'    , 'pushbutton' , ...
                                     'string'   , ''           , ...
                                     'enable'   , 'on'         , ...
                                     'visible'  , 'on'         , ...
                                     'tag'      , 'TAG_BUTTON'  );
 
   end
   
   ud.FrameHandle = uicontrol('parent'   , fig     , ...
                              'visible'  , 'off'   , ...
                              'selected' , 'off'   , ...
                              'units'    , 'pixels', ...
                              'position' ,  [ 0  0  1  1 ] , ...
                              'style'    , 'frame' , ...
                              'enable'   , 'on'         , ...
                              'visible'  , 'on'         , ...
                              'tag'      , 'TAG_FRAME'  );

                      

   form = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 );

     HF = sprintf(form,ud.FrameHandle);


   DelCB = [ 'tag_frame('  HF  ',''Delete'',-1);' ];

   set(ud.FrameHandle,'DeleteFcn',DelCB);


   % ButtonUserData

   udb = struct( 'Root'     , { fig           }  , ...
                 'Parent'   , { ud.FrameHandle } , ...
                 'Children' , { [] }                     );

   for ii = 1 : ud.TagNumber

     CB = [ 'tag_frame('  HF  ',''Activate'','  sprintf('%.0f',ii)  ',1);'  ];

     set( ud.ButtonHandle(ii) , 'callback' , CB , ...
                                'userdata' , udb          );
    
     ud.HideHandle(ii) = uicontrol(  'parent'   , fig     , ...
                                     'visible'  , 'off'   , ...
                                     'selected' , 'off'   , ...
                                     'units'    , 'pixels', ...
                                     'position' ,  [ 0  0  1  1 ] , ...
                                     'style'    , 'text'     , ...
                                     'string'   , { '' }     , ... 
                                     'enable'   , 'on'         , ...
                                     'userdata' ,   ch       , ...
                                     'tag'      , 'TAG_HIDE'  );

    
   end


   % Set DefaultProperties back
   for ii = 1 : size(VarIn,1)
     set(0,['DefaultUIControl' VarIn{ii,1}],VarIn{ii,3});
   end


   ud0 = struct( 'TAG_FRAME' , { ud  } , ...
                 'Root'      , {  0  } , ...
                 'Parent'    , { fig } , ...
                 'Children'  , { []  }         );

   set( ud.FrameHandle , 'userdata' , ud0 );

   %-------------------------------------------
   % Set correct Size of UIControls

   try
     [Msg,pos,ButtonHigh] = tag_frame(ud.FrameHandle,'Resize','absolut');
   catch
      Msg = lasterr;
   end

   if ~isempty(Msg)

     Msg = [ Msg0 'Error call TAG_FRAME( Resize ).' nl Msg ];

     set(ud.FrameHandle,'DeleteFcn','');

     delete(ud.FrameHandle);
     delete(ud.ButtonHandle);
     delete(ud.HideHandle);

     return

   end
 
   ud0.TAG_FRAME.PixelPosition = pos;
   ud0.TAG_FRAME.NormPosition  = pos ./ figpos([ 3  4  3  4 ]);

   set( ud.FrameHandle , 'userdata' , ud0 );

   %-------------------------------------------
   
   try
      Msg = tag_frame(ud.FrameHandle,'Activate');
   catch
      Msg = lasterr;
   end

   if ~isempty(Msg)
     Msg = [ Msg0 'Error call TAG_FRAME( Activate ).' nl Msg ];
     delete(ud.FrameHandle);
     delete(ud.ButtonHandle);
     delete(ud.HideHandle);
     return
   end

   %-------------------------------------------

   out = { ud.FrameHandle  ud.ButtonHandle  pos  ButtonHigh }; 

   n = min(Nout,size(out,2));

   varargout(1:n) = out(1:n);


%*********************************************************
case 'RESIZE'

  if isempty(VarIn)
    mode = 'a';   % Absolut
  else
    mode = VarIn{1};

    ok = ( ischar(mode)  &  ~isempty(mode) );
    if ok
        mode = lower(mode(1,1));
          ok = any(strcmp(mode,{'a'  'n'}));
    end
    if ~ok
      Msg = [  Msg0 'Input must be a String: ''normalized'' or ''absolut''.' ];
      return
    end
  end


  % FigurePosition in Pixels

  figuni = get(fig,'units');
           set(fig,'units','pixels');
  figpos = get(fig,'position');
           set(fig,'units',figuni);


  ButtonHigh = ceil(ud.PixelHight*ud.HighPadding) + 2*ud.BorderWidth;

  ButtonHigh = ButtonHigh * ( ud.TagNumber > 0 );

  %--------------------------------------------------------------------------
  if strcmp(mode,'n')
  % Normalized

     pos = ud.NormPosition .* figpos([3 4 3 4]);

     pos([1 3]) =  ceil( pos([1 3]) );
     pos([2 4]) = floor( pos([2 4]) );
 

  %--------------------------------------------------------------------------
  else
  % Absolut

    %  'Position'  , [  Left  Bottom  -Right -Top  ] , default: [ 5  5 -5 -5]
    %                [  Left  Bottom   Width  High ]
    %                [ -Right -Top     Width  High ]


     po = ud.Position;

     pos = ( [ 1  1  0  0 ] + po ) .* ( po >= 0 );

     pos([3 4]) = pos([3 4]) + ...
                   ( figpos([3 4]) - pos([1 2]) + 1 + po([3 4]) ) .* ...
                   ( po([1 2]) >= 0 ) .* ( po([3 4]) <= 0 );

     pos([1 2]) = pos([1 2]) + ...
                   ( figpos([3 4]) - pos([3 4]) + 1 + po([1 2]) ) .* ...
                   ( po([1 2]) < 0 ) .* ( po([3 4]) > 0 );

     pos([1 3]) =  ceil( pos([1 3]) );
     pos([2 4]) = floor( pos([2 4]) );
 
             
  end
  %--------------------------------------------------------------------------

     pos(4) = pos(4);

     ud.PixelPosition = pos;

  % Frame

     posf = pos;
     posf(4) = posf(4) - ButtonHigh + ud.BorderWidth * ( ButtonHigh > 0 );  % !!!

     
  % Buttons & Hide

  n = ud.TagNumber;

  if n > 0

     % Buttons

     pos(1) = pos(1)-0;
     pos(3) = pos(3)+0;

     posb = [ pos(1)          pos(2)+pos(4)-ButtonHigh  ...
             floor(pos(3)/n)  ButtonHigh                    ];

     posb = posb(ones(1,n),:);

     dw =  pos(3) - n * posb(1,3);

     ii = ( n-dw+1 : ud.TagNumber );

     posb(ii,3) = posb(ii,3) + 1;

     posb(2:n,1) = posb(2:n,1)+cumsum(posb(1:n-1,3));

  
     % HideHandle

     posh = posb;
     posh(:,1) = posh(:,1) + ud.BorderWidth; 
     posh(:,2) = posf(2)+posf(4)-ud.BorderWidth;
     posh(:,3) = posh(:,3) - 2 * ud.BorderWidth; 
     posh(:,4) = 1 * ud.BorderWidth;

  end
  % n > 0

  ud0.TAG_FRAME = ud;

     posf([3 4]) = max( posf([3 4]) , 1 );

     set(ud.FrameHandle , 'units'     , 'pixels' , ...
                          'position'  , posf     , ...
                          'userdata'  , ud0            )

   for ii = 1 : n

     posb(ii,[3 4]) = max( posb(ii,[3 4]) , 1 );

     set(ud.ButtonHandle(ii) , 'units'    , 'pixels' , ...
                               'position' , posb(ii,:)    )

     posh(ii,[3 4]) = max( posh(ii,[3 4]) , 1 );

     set(ud.HideHandle(ii) , 'units'    , 'pixels' , ...
                             'position' , posh(ii,:)       )

   end


  % Resize Children

    for ii = 1 : n

      ch = get(ud.HideHandle(ii),'userdata');

      for jj = 1 : size(ch.Handle,1)

         tag_frame(ud.FrameHandle,'ResizeChildren',ch.Handle(jj),ch.Position{jj},mode);

      end

    end

    for jj = 1 : size(ud.FrameChildren.Handle,1)

       tag_frame(ud.FrameHandle,'ResizeChildren', ...
            ud.FrameChildren.Handle(jj),ud.FrameChildren.Position{jj},mode);

    end

   %-------------------------------------------

   out = { ud.PixelPosition  ButtonHigh-ud.BorderWidth };  % !!!

   n = min(Nout,size(out,2));

   varargout(1:n) = out(1:n);


%*********************************************************
case 'ACTIVATE'

  clb = 0;   % CallBack

  if ~isempty(VarIn)

    nr = VarIn{1};

    ok = ( isnumeric(nr) &  ( prod(size(nr)) == 1 ) );
    if ok
       ok = ( any( nr == ( 0 : ud.TagNumber ) )  |  isnan(nr) );
       if ok
          if nr == ud.ActiveNr
             return
          end
       else
          hh = [ ud.FrameHandle ; ud.ButtonHandle ];
          ok = any( nr == hh );
          if ok
             nr = find( nr == hh );
             nr = nr(1)-1;
          end
       end
    end

    if ~ok

      Msg = [ Msg0  'Input must be a ButtonNumber or ZERO ' nl ...
                    ' or a ButtonHandle or TagFrameHandle.' ];
      return

    end

    if prod(size(VarIn)) > 1
       clb = VarIn{2};
    end

  else
     
    nr = ud.ActiveNr;

  end


    if isnan(nr)  

       sets = { 'off'  'off' };

    else

       ud.ActiveNr = nr;
 
       sets = { 'off'  ud.Visibility };

    end

    ud0.TAG_FRAME = ud;

    set(ud.FrameHandle  , 'visible'  , sets{2} , ...
                          'userdata' ,  ud0       );
    set(ud.ButtonHandle , 'visible'  , sets{2}       );


    %-----------------------------------------

    for ii = 1 : ud.TagNumber

      set( ud.HideHandle(ii) , 'visible' , sets{1+(ii==nr)} );

      ch = get(ud.HideHandle(ii),'userdata');

      ch = ch.Handle;

      set( ch , 'visible' , sets{1+(ii==nr)} );


      cht = get(ch,'tag');

        %---------------------------------------
        % Check for SelList

        tt = find( strcmp(cht,'SEL_FRAME') );

        for jj = tt(:)'

          sel_list(ch(jj),'visible',sets{1+(ii==nr)});
  
        end

        %---------------------------------------
        % Check for Tabular

        tt = find( strcmp(cht,'TAB_FRAME') );

        for jj = tt(:)'

          tab_list(ch(jj),'visible',sets{1+(ii==nr)});
  
        end

        %---------------------------------------
        % Check for TagFrame

        tt = find( strcmp(cht,'TAG_FRAME') );

        for jj = tt(:)'
    
          tag_frame(ch(jj),'visible',sets{1+(ii==nr)});
  
        end

        %---------------------------------------
        % Check for MsgLogo

        tt = find( strcmp(cht,'LOGO_FRAME') );

        for jj = tt(:)'
    
          msg_logo(ch(jj),'visible',sets{1+(ii==nr)});
  
        end

        %---------------------------------------
        % Check for MsgList

        tt = find( strcmp(cht,'MESSAGE_FRAME') );

        for jj = tt(:)'
    
          msg_list(ch(jj),'visible',sets{1+(ii==nr)});
  
        end


    end


    %-----------------------------------------
    % Check for activated Frame

    ch = ud.FrameChildren.Handle;

    set( ch , 'visible' , sets{1+(0==nr)} );


    cht = get(ch,'tag');


      %---------------------------------------
      % Check for Tabular

      tt = find( strcmp(cht,'TabFrame') );

      for jj = tt(:)'

         tab_list(ch(jj),'visible',sets{1+(0==nr)});

      end

      %---------------------------------------
      % Check for TagFrame

      tt = find( strcmp(cht,'TAG_FRAME') );

      for jj = tt(:)'
   
         tag_frame(ch(jj),'visible',sets{1+(0==nr)});
  
      end



    %-----------------------------------------
    % CallBack

     if isequal(clb,1)  &  ~isempty(ud.CBFcn)
        try
           feval( ud.CBFcn{:} , 'Activate' , 1 , nr );
        catch
           fprintf([nl Msg0 'Error call CBFcn.' nl lasterr nl ]);
        end
     end

    %-----------------------------------------

    out = { nr };

    n = min(Nout,size(out,2));

    varargout(1:n) = out(1:n);


%*********************************************************
case 'ADD'

  if ( prod(size(VarIn)) < 2 )

      Msg = [ Msg0  'Not enough Input Arguments.' ];

      return

  end

  %----------------------------------------------------------
  % Check 1. Input

    nr = VarIn{1};

    ok = ( isnumeric(nr) &  ( prod(size(nr)) == 1 ) );
    if ok
       ok = any( nr == ( 0 : ud.TagNumber ) );
       if ~ok
          hh = [ ud.FrameHandle ; ud.ButtonHandle ];
          ok = any( nr == hh );
          if ok
             nr = find( nr == hh );
             nr = nr(1)-1;
          end
       end
    end

    if ~ok

      Msg = [ Msg0  '1. Input must be a ButtonNumber or ZERO ' nl ...
                    ' or a ButtonHandle or TagFrameHandle.' ];
      return

    end

  %----------------------------------------------------------
  % Check 2. Input

    type = VarIn{2};

    ok = ( ischar(type) & ~isempty(type) & ...
           ( prod(size(type)) == size(type,2) ) );
    if ok
       type = lower(type);
       typs = {'uicontrol' 'sel_list'  'tab_list' 'tag_frame' 'msg_list' 'msg_logo'};
       ok = any(strcmp(type,typs));
    end

 
    if ~ok

      Msg = [ Msg0  '2. Input must be ' ...
                    '''uicontrol'', ''tab_list'', ''tag_frame'', ' ...
                    '''sel_list'' , msg_list'' or ''msg_logo'' .'];
      return

    end

  %----------------------------------------------------------
  % Check other Input

   VarIn = VarIn(3:end);

   if ~isempty(VarIn)

      VarIn = VarIn(:);

      is_logo = strcmp(type,'msg_logo');

      ok = is_logo;

      if ~ok

        ok = ( mod(size(VarIn,1),2) == 0 );

        if ok
  
           VarIn = reshape(VarIn,2,size(VarIn,1)/2);

              ok = iscellstr(VarIn(1,:));

        end

      
        if ~ok
  
          Msg = [ Msg0  'Following Inputs must be ' ...
                        'Property-Value-Pairs, Properties must be Strings.'];
          return

        end
 
      end

   end

  %----------------------------------------------------------
  % Build Control

   switch type

     %-------------------------------------------------------------
     case 'uicontrol'

       % Search for "position"
       if ~isempty(VarIn)

         VarIn0 = VarIn;

         ii = find( strcmp( lower(VarIn(1,:)) , 'position' ) );

         if ~isempty(ii)
           for jj = ii(:)'
             pos         = VarIn{2,jj};
             pos([3 4])  = abs(pos([3 4]));
             pos([3 4])  = pos([3 4]) + ( pos([3 4]) == 0 );
             VarIn{2,jj} = pos;
           end
         end
         
       end

       try
          h = uicontrol('parent',fig,VarIn{:});
       catch
          Msg = [ Msg0 'Error call UICONTROL.' nl lasterr ];
          return
       end

       if ~isempty(VarIn)

         VarIn = VarIn0;

       end

     %-------------------------------------------------------------
     case 'sel_list'

       try
         [Msg,h,ht] = sel_list(fig,'New',VarIn{:});
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)

          Msg = [ Msg0 'Error call SEL_LIST.' nl Msg ];
          return

       end

     %-------------------------------------------------------------
     case 'tab_list'

       try
         [Msg,h] = tab_list(fig,VarIn{:});
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)

          Msg = [ Msg0 'Error call TAB_LIST.' nl Msg ];
          return

       end


     %-------------------------------------------------------------
     case 'tag_frame'

       try
         [Msg,h,hb,p,bh] = tag_frame(fig,'New',VarIn{:});
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)
          Msg = [ Msg0 'Error call TAG_FRAME( New ).' nl Msg ];
          return

       end

       % Add Root to Buttons

         if nr == 0
            hr = ud.FrameHandle;
         else
            hr = ud.ButtonHandle(nr);
         end

         for hh = hb(:)'
           bud = get( hh , 'userdata' );
           bud.Root = hr;
           set( hh , 'userdata' , bud );
         end
  
       % Parent to FrameHandle

         hud = get( h , 'userdata' );
         hud.Parent = hr;
         hud.Root   = get( hr , 'Parent' );
       
         set( h , 'userdata' , hud );

     %-------------------------------------------------------------
     case 'msg_list'

       try
         [Msg,h] = msg_list(fig,'New',VarIn{:});
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)

          Msg = [ Msg0 'Error call MSG_LIST( New ).' nl Msg ];
          return

       end

     %-------------------------------------------------------------
     case 'msg_logo'

       try
         [Msg,hl,h,ht] = msg_logo(fig,'New',VarIn{:});
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)

          Msg = [ Msg0 'Error call MSG_LOGO( New ).' nl Msg ];
          return

       end

   end


  %----------------------------------------------------------
  % Get Position

   switch type

     %-------------------------------------------------------------
     case 'uicontrol'

       if isempty(VarIn)

         pos = get(h,'position');

       else

         ii = find( strcmp( lower(VarIn(1,:)) , 'position' ) );

         if isempty(ii)
           pos = get(h,'position');
         else
           pos = VarIn{2,ii(end)};
         end
         
       end

     %-------------------------------------------------------------
     case 'sel_list'
   
       tud = get(h,'userdata');
    
       pos = tud.Position;

     %-------------------------------------------------------------
     case 'tab_list'
   
       tud = get(h,'userdata');
    
       pos = tud.position;

     %-------------------------------------------------------------
     case 'tag_frame'
   
       tud = get(h,'userdata');
    
       pos = tud.TAG_FRAME.Position;
   
     %-------------------------------------------------------------
     case 'msg_list'
   
       tud = get(h,'userdata');
    
       pos = tud.PixelOffset;
   
     %-------------------------------------------------------------
     case 'msg_logo'
   
       tud = get(h,'userdata');
    
       pos = tud.Offset;

   end


  %-----------------------------------------
  % Resize

    Msg = tag_frame( ud.FrameHandle,'ResizeChildren',h,pos,'absolut',nr);

    if ~isempty(Msg)

        delete(h);

        return

    end

  %-----------------------------------------
  % Add new Children
 
    if nr == 0

      ud.FrameChildren.Handle   = cat( 1 , ud.FrameChildren.Handle , h );
      ud.FrameChildren.Position = cat( 1 , ud.FrameChildren.Position , {pos} );

      ud0.TAG_FRAME = ud;

      set( ud.FrameHandle , 'userdata' , ud0 );

    else

      ch = get(ud.HideHandle(nr),'userdata');
 
      ch.Handle = cat(1,ch.Handle,h);

      ch.Position = cat( 1 , ch.Position , {pos} );

      set( ud.HideHandle(nr) , 'userdata' , ch );

    end


  %-----------------------------------------

    if strcmp(type,'msg_logo');
       out = { hl h ht };        % MessageListHandle, LogoFrameHandle, LogoTextHandle
    elseif strcmp(type,'msg_logo');
       out = { h ht };           % FrameHandle TextHandle
    elseif strcmp(type,'tag_frame')
       tud = get(h,'userdata');
       out = { h  hb tud.TAG_FRAME.PixelPosition bh }; % TagFrameHandle TagButtonHandle
    else
       out = { h };
    end

    n = min(Nout,size(out,2));

    varargout(1:n) = out(1:n);


%*********************************************************
case 'RESIZECHILDREN'

  if prod(size(VarIn)) < 2

    Msg = [ Msg0 'ChildrenHandle and Position required.' ];
    return
 
  end


  h = VarIn{1};

%  ok = ( isnumeric(h)  &  ( prod(size(h)) == 1 ) );
%  if ok 
%     ok = ishandle(h);
%     if ok
%        ok = strcmp( get(h,'type') , 'uicontrol' );
%     end
%  end
%  if ~ok
%      Msg = [ Msg0  'Input must be a ChildrenHandle.' ];
%      return
%  end


  pos = VarIn{2};

  mode = lower(VarIn{3}(1));

  if prod(size(VarIn)) < 4

     nr = ud.ActiveNr;

  else

    nr = VarIn{4};

%    ok = ( isnumeric(nr) &  ( prod(size(nr)) == 1 ) );
%    if ok
%       ok = any( nr == ( 0 : ud.TagNumber ) );
%       if ~ok
%          hh = [ ud.FrameHandle ; ud.ButtonHandle ];
%          ok = any( nr == hh );
%          if ok
%             nr = find( nr == hh );
%             nr = nr(1)-1;
%          end
%       end
%    end
%    if ~ok
%      Msg = [ Msg0  '2. Input must be a ButtonNumber or ZERO ' nl ...
%                    ' or a ButtonHandle or TagFrameHandle.' ];
%      return
%    end

  end

 
  %-----------------------------------------------------
  % FramePosition

    % FigurePosition in Pixels

     figuni = get(fig,'units');
              set(fig,'units','pixels');
     figpos = get(fig,'position');
              set(fig,'units',figuni);

    fpos        = get(ud.FrameHandle,'position');

    opos        = zeros(1,4);  % [  Left  Bottom  -Right -Top  ]
    opos([1 2]) = fpos([1 2]) - 1 + ud.BorderWidth;
    opos([3 4]) = - ( figpos([3 4]) - ( fpos([1 2]) + fpos([3 4]) - 1 ) + ...
                      ud.BorderWidth );
 
             ff = figpos([3 4])./ud.FigPos0([3 4]);
 

  switch(get(h,'tag'));

     %-------------------------------------------------------------
     case 'TAG_FRAME'

    %  'Position'  , [  Left  Bottom  -Right -Top  ] , default: [ 5  5 -5 -5]
    %                [  Left  Bottom   Width  High ]
    %                [ -Right -Top     Width  High ]


       if strcmp(mode,'n');

          pos([1 3]) = pos([1 3]) * ff(1);
          pos([2 4]) = pos([2 4]) * ff(2);

       end
  
       pos([1 2]) = pos([1 2]) + opos([1 2]) .* ( pos([1 2]) >= 0 ) + ...
                                 opos([3 4]) .* ( pos([1 2]) <  0 );
       pos([3 4]) = pos([3 4]) + opos([3 4]) .* ( pos([3 4]) <= 0 );

       tud = get(h,'userdata');

       tud.TAG_FRAME.Position = pos;

       set(h,'userdata',tud);
 
       try
          Msg = tag_frame(h,'Resize','absolut');
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)

          Msg = [ Msg0 'Error call TAG_FRAME( Resize ).' nl Msg ];
        
       end


       if nr ~= ud.ActiveNr
 
          tag_frame( h , 'visible' , 'off' );

       end

     %-------------------------------------------------------------
     case 'SEL_FRAME'

       if strcmp(mode,'n');

          pos = pos .* ff;
 
       end

       pos([1 2]) = pos([1 2]) + opos([1 2]) .* ( pos([1 2]) >= 0 ) + ...
                                 opos([3 4]) .* ( pos([1 2]) <  0 );

       tud = get(h,'userdata');

       tud.Position = pos;

       set(h,'userdata',tud);

       try
          Msg = sel_list(h,'Resize','absolut');
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)

          Msg = [ Msg0 'Error call SEL_LIST( RESIZE ).' nl Msg ];

       end

       if nr ~= ud.ActiveNr
 
          sel_list( h , 'visible' , 'off' );

       end

     %-------------------------------------------------------------
     case 'TAB_FRAME'

       if strcmp(mode,'n');

          pos = pos .* ff;
 
       end

       pos([1 2]) = pos([1 2]) + opos([1 2]) .* ( pos([1 2]) >= 0 ) + ...
                                 opos([3 4]) .* ( pos([1 2]) <  0 );

       try
          Msg = tab_list(h,'position',pos);
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)

          Msg = [ Msg0 'Error call TAB_LIST( Position ).' nl Msg ];

       end

       if nr ~= ud.ActiveNr
 
          tab_list( h , 'visible' , 'off' );

       end

     %-------------------------------------------------------------
     case 'MESSAGE_FRAME'

        % pos = [  Left  Right  Top     ] 
        %       [  Left  Right -Bottom  ]
 
       if strcmp(mode,'n');

          pos([1 2]) = pos([1 2]) * ff(1);
          pos([ 3 ]) = pos([ 3 ]) * ff(2);

       end

       pos([1 2]) = pos([1 2]) + opos([1 3]).*[1 -1];
       pos([ 3 ]) = ( pos(3) - opos(4) ) * ( pos(3) >= 0  ) + ...
                    ( pos(3) - opos(2) ) * ( pos(3) <  0  );

       tud = get(h,'userdata');

       tud.PixelOffset = pos;

       set(h,'userdata',tud);

       try
          Msg = msg_list(h,'resize','absolut');
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)

          Msg = [ Msg0 'Error call MSG_LIST( Resize ).' nl Msg ];

       end

       if nr ~= ud.ActiveNr
 
          msg_list( h , 'visible' , 'off' );

       end


     %-------------------------------------------------------------
     case 'LOGO_FRAME'

       %  pos.MsgOffset  =  [ Left Right Top ]
       %  pos.LogoOffset =  [ Left  Right ]

       if strcmp(mode,'n');

          pos.Msg([1 2]) = pos.Msg([1 2]) * ff(1);
          pos.Msg([ 3 ]) = pos.Msg([ 3 ]) * ff(2);

          pos.Logo = pos.Logo * ff(1);

       end

       pos.Logo(1) = pos.Logo(1) + opos(1);
       pos.Msg(2)  = pos.Msg(2) - opos(3);
       pos.Msg(3)  = ( pos.Msg(3) - opos(4) ) * ( pos.Msg(3) >= 0  ) + ...
                        ( pos.Msg(3) - opos(2) ) * ( pos.Msg(3) <  0  );

       tud = get(h,'userdata');

       tud.Offset = pos;

       set(h,'userdata',tud);

       lud = get(tud.Handle.List,'userdata');
       lud.PixelOffset = pos.Msg;
       set(tud.Handle.List,'userdata',lud);

       try
          Msg = msg_logo(h,'resize','absolut');
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)

          Msg = [ Msg0 'Error call MSG_LOGO( Resize ).' nl Msg ];

       end

       if nr ~= ud.ActiveNr
 
          msg_logo( h , 'visible' , 'off' );

       end



     %-------------------------------------------------------------
     otherwise
     % uicontrol

    %  'Position'  , [  Left  Bottom  -Right -Top  ] , default: [ 5  5 -5 -5]
    %                [  Left  Bottom   Width  High ]
    %                [ -Right -Top     Width  High ]

       if strcmp(mode,'n');

          pos([1 3]) = pos([1 3]) * ff(1);
          pos([2 4]) = pos([2 4]) * ff(2);

       end

       poo  = pos;

       pos([1 2]) = pos([1 2]) + opos([1 2]) .* ( pos([1 2]) >= 0 ) + ...
                                 opos([3 4]) .* ( pos([1 2]) <  0 );
       pos([3 4]) = pos([3 4]) + opos([3 4]) .* ( pos([3 4]) <= 0 );


       po = pos;

       pos = ( [ 1  1  0  0 ] + po ) .* ( po >= 0 );

       pos([3 4]) = pos([3 4]) + ...
                     ( figpos([3 4]) - pos([1 2]) + 1 + po([3 4]) ) .* ...
                     ( po([1 2]) >= 0 ) .* ( po([3 4]) <= 0 );

       pos([1 2]) = pos([1 2]) + ...
                     ( figpos([3 4]) - pos([3 4]) + 1 + po([1 2]) ) .* ...
                     ( po([1 2]) < 0 ) .* ( po([3 4]) > 0 );

       uni = get(h,'units');

       pos([3 4]) = max( pos([3 4]) , 1 );

       set( h , 'units'    , 'pixels' , ...
                'position' , pos            );

       set(h,'units',uni);

       if nr ~= ud.ActiveNr
 
         set(h,'visible','off');

       end

 end


%***********************************************************************
case { 'VISIBLE'  'ENABLE' }

   
  ok = ~isempty(VarIn);

  if ok

    sets = VarIn{1};

    ok = ( ischar(sets) & ~isempty(sets) & ...
           ( prod(size(sets)) == size(sets,2) ) );
    if ok
       ok = any(strcmp(sets,{'on'  'off'}));
    end

  end

  if ~ok

    Msg = [ Msg0  'Input for VISIBLE must be ''on'' or ''off''.'];
    return

  end

  %-------------------------------------------------- 
  if strcmp(action,'VISIBLE')

    ud.Visibility = sets;

    ud0.TAG_FRAME = ud;

    set ( ud.FrameHandle , 'userdata' , ud0 );

    if strcmp(sets,'off')

       tag_frame( ud.FrameHandle , 'Activate' , NaN );
  
    else

       tag_frame( ud.FrameHandle , 'Activate' );

    end
     
    return

  end

  %--------------------------------------------------
  if size(VarIn,1) > 1

    nr = VarIn{2};

    ok = ( isnumeric(nr) &  ( prod(size(nr)) == 1 ) );
    if ok
       ok = ( any( nr == ( 0 : ud.TagNumber ) )  |  isnan(nr) );
       if ~ok
          hh = [ ud.FrameHandle ; ud.ButtonHandle ];
          ok = any( nr == hh );
          if ok
             nr = find( nr == hh );
             nr = nr(1)-1;
          end
       end
    end

    if ~ok

      Msg = [ Msg0  'Input must be a ButtonNumber or ZERO ' nl ...
                    ' or a ButtonHandle or TagFrameHandle.' ];
      return

    end

  else
     
    nr = ud.ActiveNr;

  end


  %---------------------------------------------------------
  if strcmp(sets,'off')  &  any( nr == [ ud.ActiveNr  0 ] )

         tag_frame( ud.FrameHandle , 'Activate' , 0 );

  end


  hh = ud.ButtonHandle;

  if nr > 0

     hh = hh(nr);

  end

  set( hh , 'enable' , sets );


%***********************************************************************
case 'DELETE'

   
  % Delete Children

    for ii = 1 : ud.TagNumber

      if ishandle(ud.HideHandle(ii))

        ch = get(ud.HideHandle(ii),'userdata');

        ch = ch.Handle;

        for h = ch(:)'

           try
             delete(h);
           end

        end

      end

    end


    hh = [ ud.ButtonHandle ; ud.HideHandle ; ud.FrameChildren.Handle ];

    for h = hh'

         try
           delete(h);
         end

    end
    
    ok = isempty(VarIn);
    if ~ok
       ok = ~isequal(VarIn{1},-1);
    end

    if ok
       set(ud.FrameHandle,'deletefcn','');
       delete(ud.FrameHandle)
    end

end


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

