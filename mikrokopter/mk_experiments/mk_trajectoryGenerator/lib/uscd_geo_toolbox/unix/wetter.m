function msg = wetter(h,grp,act,clb)

% WETTER  a UNIX-Browser for actual WeatherGrafics of IfM Kiel
%
% "http://www.ifm.uni-kiel.de/fb/fb1/me/kieldaten/kieldata-d.html"
%
% Configuration at the end of this M-File:
%
% URL:         String for URL
% EXT:         String for ImageExtension in URL 
%
% Location:  { Label  Prefix }        [ N by 2 ]
% Time:      { Label SubDir Suffix }  [ M by 3 ]
% Type:      { Label Nr1 Nr2 ...   }  [ K by N ]
%
% Scale:       CharArray z.B. BeaufortScale
%
% Logo:      { Logo }                 [ 3 by 1 ]
% Size:       [ Width  High ]
%
% The Adress is build like: 
%
%       fullfile( URL , SubDir , [Prefix Nr Suffix EXT] )
%
%---------------------------------------------------------------------------
%
% requires: MKGUI, STDGUI, Unix: wget, convert
%

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
app = upper(fcn);

name = 'Wetter Kiel';

Nin  = nargin;
Nout = nargout;

%**********************************************************
% New Figure

if Nin == 0
   msg = newfig(fcn,app,name);
   if ~isempty(msg) & ( Nout == 0 )
       error(msg)
   end
   return
elseif Nin < 3
   return
end

%**********************************************************
% Check Handle

ok = ( isnumeric(h) & ( prod(size(h)) == 1 ) );
if ok
   ok = ishandle(h);
   if ok
      [msg,fig] = recpar(h,'figure');
      ok = ( isempty(msg) & ~isempty(fig) );
      if ok 
         fig = fig(1);
         ok = isappdata(fig,app);
      end
   end
end

if ~ok
    msg = 'Invalid Handle.'
    if ( Nout == 0 )
       error(msg)
    end
    return
end

%**********************************************************
% Check Action

if ~chkstr(act,1)
    msg = 'Action must be a String.';
    if ( Nout == 0 )
       error(msg)
    end
    return
end

%**********************************************************

ud = get(fig,'userdata');

ch = get(ud.Children.Frame,'userdata');
hl = ch.Children.Logo.Message;

ch = ch.Children.Control.Main;

%**********************************************************

switch upper(act)

      case { 'LOCATION'  'TIME'  'TYPE' }

         v = getappdata(fig,app);

         % URL
         % Location: { Label  Prefix }
         % Time:     { Label SubDir Suffix }
         % Type:     { Label NrLoc1 NrLoc2 }
%

         loc  = get(ch.Location,'value');
         time = get(ch.Time,'value');
         typ  = get(ch.Type,'value');

         nr = cat(2,v.Type{typ,2:end});
         if isnan(nr(loc))
            if all(isnan(nr))
               return
            end
            loc = min(find(~isnan(nr)));
            nr = nr(loc);
         else
            nr = nr(loc);
         end

         file = sprintf('%s%.0f%s%s',v.Location{loc,2},nr,v.Time{time,3},v.EXT);

         url = fullfile(v.URL,v.Time{time,2},file);

         setimage(url,ch.Image,hl);

      case 'IMAGE'

         url = get(ch.Image,'userdata');
         if ~chkstr(url,1)
             return
         end

         setimage(url,ch.Image,hl);

      case 'SCALE'

         set(ch.Scale,'value',[]);

end

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function setimage(url,hc,hl);


[p,n,e] = fileparts(url);

tmp = sprintf('%s_%s.png',tempname,n);

command = sprintf('wget -O - %s | convert - %s',url,tmp);

msg_list(hl,'Message',sprintf('Call UNIX: %s',command),'wait');

%%% msg = lib_kim(0);

[s,w] = unix(command);

%%% msg = lib_kim(1);

ok = ( s == 0 );
if ok
   ok = ( exist(tmp,'file') == 2 );
   if ok
      d = dir(tmp);
      ok = ( prod(size(d)) == 1 );
      if ok
         ok = ( ~( d.bytes == 0 ) & ~d.isdir );
      end
   end
end

if ~ok
    msg_list(hl,'Message','error','append');
    if ~isempty(w)
        [m,w] = char2cell(rmblank(w,2));
        w = strhcat(w,' ',2);
        msg_list(hl,'Message',w,'append',char(10));
    end
    if ( exist(tmp,'file') == 2 )
       try, delete(tmp); end
    end
    return
end

msg_list(hl,'Message','ok','append');

%***********************************************************
% Read Image

msg_list(hl,'Message',sprintf('Call IMREAD: %s',tmp),'wait');

msg = '';
try
   [c,m] = imread(tmp);
   if isempty(c)
      msg = 'Empty Image.';
   end
catch
   msg = lasterr;
end

try, delete(tmp); end

if ~isempty(msg)
    msg_list(hl,'Message',sprintf('error\n%s',msg),'append');
    return
end

msg_list(hl,'Message','ok','append');

%***********************************************************
% Transform Data

cl = class(c);

if size(c,3) == 3
   if strcmp(cl,'uint8')
      c = double(c) / 255;
   end
else
   if isempty(m)
      if strcmp(cl,'uint8')
         c = double(c) / 255;
      end
      c = c(:,:,ones(1,3));
   else
      if strcmp(cl,'uint8')
         c = double(c) + 1;
      end
      n = size(m,1);
      c = m(cat(3,c,c+n,c+2*n));
   end
end

set(hc,'cdata',c,'userdata',url);

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function msg = newfig(fcn,app,name);

v = config;

bg = [ 0.90 0.95 1.0 ];

[msg,fig] = mkgui(gui_config(fcn,v),name,bg);

ud = get(fig,'userdata');

ch = get(ud.Children.Frame,'userdata');
ch = ch.Children.Control.Main;

set(ch.Image,'style','pushbutton','horizontalalignment','center'),

set(ch.Scale,'value',[]);

setappdata(fig,app,v);

set(fig,'tag',app,'visible','on');


%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function msg = lib_kim(mode)

% LIB_KIM  Renames Links of Matlab-LIB-files
%
% Nessecary to open an Window under X using Matlab5 on SuSE9
%
% Works on HOST kim, VERSION 5, BINDIR /opt/matlab53/bin/lnx86
%
% Changes: libXt.so.6  <---> libXt.so.6.tmp
%          libX11.so.6 <---> libX11.so.6.tmp
%
% LIB_KIM( Mode )
%
% Mode:    0   --->
%          1  <---
%
% 

host = 'kim';
vers = '5';

bin = '/opt/matlab53/bin/lnx86';

lnk = { 'libX11.so.6' 
        'libXt.so.6'  };

suf = '.tmp';

n = size(lnk,1);

%**************************************************************

if nargin == 0
   error('Mode required.');
end

mode = isequal(mode,1);

%**************************************************************
% Check for HOST, VERSION, BINDIR, FILES

msg = cell(0,1);

if ~strcmp(getenv('HOSTNAME'),host)
    msg = cat(1,msg,{sprintf('HOSTNAME must be "%s"',host)});
end

v = version;

if ~( v(1) == vers )
    msg = cat(1,msg,{sprintf('VERSION must be "%s"',vers)});
end

if ~( exist(bin,'dir') == 7 )
    msg = cat(1,msg,{sprintf('Directory doesn''t exist: %s',bin)});
else
   for ii = 1 : n
       lnk{ii} = fullfile(bin,lnk{ii});
       f = cat( 2 , lnk{ii} , suf( 1 : end*mode ) );
       if ~( exist(f,'file') == 2 )
           msg = cat(1,msg,{sprintf('File doesn''t exist: %s',f)});
       end
   end
end

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    if nargout == 0
       warning(msg);
       clear msg
    end
    return
end

%**************************************************************

for ii = 1 : n

    src = cat( 2 , lnk{ii} , suf( 1 : end*mode ) );
    dst = cat( 2 , lnk{ii} , suf( 1 : end*(~mode) ) );

    cmd = sprintf('mv %s %s',src,dst);

    [s,w] = unix(cmd);

    if ~( s == 0 )
        msg = cat(1,msg,{sprintf('Error call UNIX: %s\n%s',cmd,w)});
    end

end

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    if nargout == 0
       warning(msg);
    end
else
    msg = '';
end

if nargout == 0
   clear msg
end


%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function Config = gui_config(fcn,v)

cl = [ 0.0  0.0   2/3 ];  % LogoBack

Config = struct( 'Logo'    , { { v.Logo  1-cl cl } } , ...
                 'Control' , {  gui_control(fcn,v) }        );

%**************************************************************
function Config = gui_control(fcn,v)

ww = 0;

for ff = { 'Location'  'Time'  'Type'  'Scale' }
    f  = getfield(v,ff{1});
    if iscell(f)
       f = char(f(:,1));
    end
    ww = max( ww , size(f,2) );
end

ww = ww - 2;

str = v.Location(:,1);
str = cellstr(cat(2,char(32*ones(size(str))),char(str)));

Location = stdgui('Type','listbox','Width',ww+size(str,1)*i,'String',str,...
                  'Text','Mess-Station','CBFcn',fcn);

str = v.Time(:,1);
str = cellstr(cat(2,char(32*ones(size(str))),char(str)));

Time = stdgui('Type','listbox','Width',ww+size(str,1)*i,'String',str,...
              'Text','Zeit','CBFcn',fcn);

str = v.Type(:,1);
str = cellstr(cat(2,char(32*ones(size(str))),char(str)));

Typ = stdgui('Type','listbox','Width',ww+size(str,1)*i,'String',str,...
             'Text','Messwert','CBFcn',fcn);

str  = v.Scale;
Scale = stdgui('Type','listbox','Width',ww+size(str,1)*i,'String',str,...
              'Text','Beaufort-Skala','CBFcn',fcn,'Value',[1 0 2]);
 
str = sprintf('%s %s','Click on a Selection to load the Image.', ...
                      'Click on the Image to reload it.' );

Img = stdgui('Type','listbox','Width',-v.Size(1)-v.Size(2)*i,'CBFcn',fcn, ...
             'String',str);


Sep = stdgui('Type' , 'separator' , 'Width' , ww , ... 
             'Color' , [0 0 0] , ...
             'String' , '' , 'Style' , 'List' );
                 
Config = struct( 'Location' , { { 1+0*i Location } } , ...
                 'Time'     , { { 1+1*i Time   } } , ...
                 'Type'     , { { 1+2*i Typ    } } , ...
                'Seperator' , { { 1+3*i Sep    } } , ...
                 'Scale'    , { { 1+4*i Scale  } } , ...
                 'Image'    , { { 2+0*i Img } }       );
                 
Config = struct( 'Main' , { Config } );

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function v = config

%
% URL:         String for URL
% EXT:         String for ImageExtension in URL 
%
% Location:  { Label  Prefix }        [ N by 2 ]
% Time:      { Label SubDir Suffix }  [ M by 3 ]
% Type:      { Label Nr1 Nr2 ...   }  [ K by N ]
%
% Scale:       CharArray z.B. BeaufortScale
%
% Logo:      { Logo }                 [ 3 by 1 ]
% Size:       [ Width  High ]
%
% The Adress is build like: 
%
%       fullfile( URL , SubDir , [Prefix Nr Suffix EXT] )
%

v = struct(  'URL'  , { '' } , ...
             'EXT'  , { '.gif' } , ...
        'Location'  , { {} } , ...
            'Time'  , { {} } , ...
            'Type'  , { {} } , ...
            'Scale' , { '' } , ...
            'Logo'  , { {} } , ...
            'Size'  , { [ 600 600 ] } );


v.URL = 'http://www.ifm.uni-kiel.de/fb/fb1/me/images/kieldaten/';

v.Logo = { 'Wetter'
           'Aktuelle Daten'
           'www.ifm.uni-kiel.de' };

%                Label                    Prefix
v.Location = { 'Institut & Leuchtturm'  'see-ifm'  
               'Institut'               'ifm'      };


%           Label               SubDir        Suffix
v.Time = { 'letzte Woche'       '7days'       '-woche'
           'gestern & heute'    '2days'       '-2tage'
           'gestern'            'SavePlots'   ''
           'heute'              ''            ''        };

%           Label                  NrLoc1  NrLoc2
v.Type = { 'Windgeschwindigkeit'     2       2
           'Windrichtung'            1       1
           'Luftdruck'              NaN      7
           'Pegel'                  NaN      5
           'Lufttemperatur'          3       3
           'Wassertemperatur'        4       4
           'Solare Einstrahlung'     5       6
           'Relative Luftfeuchte'    8       8   };

v.Scale = [ ...
            ...
' bft |    m/s |  km/h |   kn '
'-----------------------------'
'  0  | <  0.2 | <   1 | <  1 '
'  1  | -  1.5 | -   5 | -  3 '
'  2  | -  3.3 | -  11 | -  6 '
'  3  | -  5.4 | -  19 | - 10 '
'  4  | -  7.9 | -  28 | - 16 '
'  5  | - 10.7 | -  38 | - 21 '
'  6  | - 13.8 | -  49 | - 27 '
'  7  | - 17.1 | -  61 | - 33 '
'  8  | - 20.7 | -  74 | - 40 '
'  9  | - 24.4 | -  88 | - 47 '
' 10  | - 28.4 | - 102 | - 55 '
' 11  | - 32.6 | - 117 | - 63 '
' 12  | > 32.7 | > 118 | > 64 ' ];

