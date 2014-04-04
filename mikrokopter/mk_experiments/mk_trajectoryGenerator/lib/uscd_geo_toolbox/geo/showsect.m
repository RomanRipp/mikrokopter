function [msg,dat,fig] = showsect(cnf,varargin)

% SHOWSECT  Extract and Display a bathymetric Section 
%
% SHOWSECT requires GET_TOPO to extract and TOPOSURF to display
%  the bathymetric Section.
%
% [Msg,Section,Fig] = SHOWSECT( CNF , Pos1  ,  Pos2 )
% [Msg,Section,Fig] = SHOWSECT( CNF , LON , LAT )
%
% Shows the Bathymetry between the Locations Pos1 and Pos2,
%  or at the section defined by the Vectors LAT and LON.
%
%  CNF = Configuration for Bathymetric Database, see GET_TOPO
%
% Pos1 = [ Lon1 Lat1 ], Pos2 = [ Lon2 Lat2 ]
%
%     NetCDF-Config: { CNF_File  [DataVar]  [NC_FileName] }
%
%        MAT-Config: { MAT_File  [XVar]  [YVar]  [ZVar] }
%
% Section = [ Lon Lat Depth Distance ]
%
% SHOWSECT use GET_TOPO to extract the Section and TOPOSURF for display.
%
% SHOWSECT( CNF , ... , Stride , SmoothWindow , Step )
%
% defines additional Parameter for GET_TOPO.
%
% A Title for the plot can be given by a String as last Input:
%
% SHOWSECT( CNF , ... , Name )
%
% See also: GET_TOPO, TOPOSURF, LOAD_CDF, CNF_CDF
% 

Nout = nargout;

msg = '';
fig = [];
dat = [];

if nargin < 2
   error('Not enough InnputArguments.');
end

%----------------------------------------------------
% Check for AxesType and Name

vin = varargin;

name = 0;
btax = 0;
cnt  = 0;

while ~isempty(vin) & ( cnt < 2 )
    cnt = cnt + 1;
    v = vin{end};
    isc = ( iscell(v) & isempty(v) );
    if ~( ischar(v) | iscellstr(v) | isc )
        break
    end
    ok = ( chkstr(vin{end},1) & isequal(btax,0) );
    if ok
       ok = ( size(v,2) <= 2 );
       if ok
          btax = v;
       end
    end
    if ~ok & isequal(name,0)
        ok = isc;
        if ok
           name = v;
        else
           [ok,v] = chkcstr(v,0);
           if ok
              name = v;
           end
        end
    end
    if ~ok
        warning('Not recognized String-Input.')
    end
    vin = vin(1:(end-1));
end

if isequal(btax,0)
   btax = 'd';
end

if isequal(name,0)
   name = {'Bathymetric Section'};
end

%----------------------------------------------------
% Get Section 

try
   dat = get_topo(cnf,vin{:});
catch
   msg = lasterr;
end

if ~isempty(msg)
    msg = sprintf('Error call GET_TOPO.\n%s',msg);
    if Nout == 0
       error(msg)
    end
    return
end

dat(:,4) = dat(:,4) / 1000;  % [m] --> [km]

%****************************************************
% Prepare Data

%-----------------------------------------------------
% Check for Y-Achses

yl = [ min(dat(:,3))  max(dat(:,3)) ];

yl(2) = yl(2) * ( yl(1) > 0 );
yl(1) = yl(1) * ( yl(1) < 0 );

yl(2) = yl(2) * ( yl(2) > 0 );

ng = ( yl(1) < 0 );

lb = { 'Elevation'  'Depth' };
yd = { 'normal'  'reverse' };

lb = lb{ 1 + ng };
yd = yd{ 1 + ng };

sg = 1 - 2 * ng;

yl = yl * sg;

ni = [  2   5   10 ];

ym = max(yl);
pl = floor( log(ym) / log(10) );
nl = ceil( ym / 10^pl );
nl = sum( ni < nl ) + 1;
nl = min( nl , size(ni,2) );
dl = ni(nl) * 10^(pl-1);

yl = [ 0   dl*(ceil(ym/dl)+1) ];

y0 = yl(1+ng);

%-----------------------------------------------------
% Check for X-Achses

rev = ( btax(1) == upper(btax(1)) );

btax = lower(btax);

xd = { 'normal'  'reverse' };

xd = xd{ 1 + rev };

xlab = { 'Distance   [km]'   'Longitude'   'Latitude' }; % XLabels
gmod = { ''                  'x'           'X'        }; % GEOAXIS

xb = { 'none'  'lon'  'lat' };  % XBase for TOPOSURF
xc = [   4       1      2   ];  % Column-Index in Section-Data DAT

ix = 1 + 1 * ( btax(1) == 'x' ) + 2 * ( btax(1) == 'y' );

xb  = xb{ix};
xi  = xc(ix);
xlb = xlab{ix};

xl  = dat([1 end],xi)';

if xl(1) == xl(2)
   xl = xl + [ -1  1 ] * 1e-6;
   dat(:,xi) = linspace(xl(1),xl(2),size(dat,1))';
elseif xl(1) > xl(2)
   xl = xl([2 1]);
end

% Check for TopAxis

if size(btax,2) == 1
   xlab = '';
   gmod = '';
   xlm = [];
else 
   ix = 1 + 1 * ( btax(2) == 'x' ) + 2 * ( btax(2) == 'y' );
   xlab = xlab{ix};
   gmod = gmod{ix};
   xlm  = dat([1 end],xc(ix))';

   if xlm(1) == xlm(2)
      xlm = xlm + [ -1  1 ] * 1e-6;
      dat(:,xi) = linspace(xlm(1),xlm(2),size(dat,1))';
   elseif xlm(1) > xlm(2)
      xlm = xlm([2 1]);
   end
end

iszoom = { 'on'  'off' };

iszoom = iszoom{ 1 + ~isempty(xlab) };

%-----------------------------------------------------
% Check Title

isc = ( iscell(name) & isempty(name) );

if ~isc

    if isempty(name{1})
       name{1} = ' ';    % Use NonEmpty Title in TOPOSURF
    end

    if ~isempty(xlab)    % If TopAxis !!!
        name = cat(1,name(:),{' ';' '});
    end

    name = sprintf('%s\n',name{:});
    name = name(1:(end-1));

elseif ~isempty(xlab)
    
    name = sprintf('%s\n%s',' ',' ');       % NonEmpty  Title if TopAxis !!!

end


%****************************************************

[fig,apd] = toposurf( []                 , ...
          'Verbose'  , 'off'             , ...
          'Zoom'     , iszoom            , ... 
          'Area'     , [ xl  yl ]        , ...
          'Mercator' , 'none'            , ...
          'Geometry' , [ 800 400 ]       , ...
          'FrameOffset' , [ 0 20 0 0 ]   , ...
          'XBase'    , xb                , ...
          'XLabel'   , xlb               , ...
          'YLabel'   , [lb '   [m]']     , ...
          'FontSize' , 12  , ...
          'YDir'     , yd , ...
          'XDir'     , xd , ...
          'Title'    , name );

h = toposurf( fig , 'patch' , dat([1 (1:end) end],xi) , ...
                        cat(1,y0,sg*dat(:,3),y0) , ...
              'edgecolor','k','facecolor','k','tag','BathyPatch');

set( h , 'buttondownfcn' , get(get(h,'parent'),'buttondownfcn') );

if ~isempty(xlab)

    drawnow

    set(0,'showhiddenhandles','on');

    set(apd.axe.Data,'box','off');

    a1 = copyobj(apd.axe.Data,fig);
    delete(get(a1,'children'));
    set(a1,'xtick',[],'yticklabel',{},'yaxislocation','right');

    set(findobj(a1,'type','text'),'string','');

    ax = copyobj(apd.axe.Data,fig);
    delete(get(ax,'children'));
    set(ax,'ytick',[], ...
        'xtickmode','auto','xticklabelmode','auto', ...
        'xlim',xlm,'xaxislocation','top');

    set(findobj(ax,'type','text'),'string','');

    set(get(ax,'xlabel'),'string',xlab, ...
        'fontname'  ,get(apd.obj.XLabel,'fontname'), ...
        'fontsize'  ,get(apd.obj.XLabel,'fontsize'), ...
        'fontweight',get(apd.obj.XLabel,'fontweight'), ...
        'fontangle' ,get(apd.obj.XLabel,'fontangle')      );

    if ~isempty(gmod)
        geoaxis(ax,gmod);
    end

end

if isc
   return
end

ht = { 'left'  'right' };
pp = [ 0  1 ] + [ 1  -1 ] * 0.005;

ps = dat([1 end],[2 1]);  % [ Lat  Lon ]
ps = ps([1 2]+[1 -1]*rev,:);

for ii = 1 : 2

     h = copyobj(apd.obj.Title,apd.axe.Data);

     set(h,'units','normalized');

     p = get(h,'position'); 

     p(1) = pp(ii);

     [v,s] = str2val(ps(ii,:),'gpos2');

     set(h,'position',p,'string',s,'horiz',ht{ii},'visible','on');

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
   n = size(str,1);
   if n == 1
      str = strrep(str,char(32),char(1));
   end
   str = cellstr(str);
   if n == 1
      str = strrep(str,char(1),char(32));
   end
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


