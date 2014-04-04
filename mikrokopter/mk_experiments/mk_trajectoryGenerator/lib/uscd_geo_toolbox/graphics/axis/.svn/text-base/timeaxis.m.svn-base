function [Msg,axt,lbt] = timeaxis(varargin)

% TIMEAXIS  Sets Axes-TickLabels in Time/Date-Format
%
% The Ticks and Labels will set depending on the actual Limit of the Axis.
% The AxisLimits will interpreted in units days !!!
%
%
% TIMEAXIS( AxeHandle , XY , Number , BaseDate )
%
%   AxeHandle   Handle of Axes, which TickLabels are to set
%                default: gca
%
%   XY         'x'  |  'y'  |  'z', specifies Axis
%                default: 'x'
%
%              An uppercase letter for XY will use a long Format for TickLabels
%               i.e. add the Date to hourly and the Year to daily TickLables.
%               No extra Labelext will displayed.
%
%  Number      Requested Number of Intervalls between Labels
%                 [ LabelNumber  TickNumber ]
%               default:  [8], the TickNumber will set automaticly
%               maximum LabelNumber: 50
%
%  BaseDate    BaseDate of Axis: [ YY MM DD ]
%                                [ YY MM DD hh mm ss ]
%               default: [ 0000 01 00  00 00 00 ]
%                          00-Jan-0000 00:00:00, i.e. this Date has day == 0 
% 
%
% [Msg,LabelAxe,LabelText] = TIMEAXIS( ... )
%
%  Msg         ErrorMessage from Zoom-Functionality, empty if everything Ok
%
%  LabelAxe    AxesHandle for a NEW created Axes, which contains the Labels.
%               Take care, that the function "gca"  can returns this Handle instead
%               of the former "AxeHandle".
%              Tag: TIMEAXIS
%
%  LabelText  TextHandles for NEW created AxesLabels, 
%              which contain the strings for Year or Date
%             Tag: TIMEAXIS_LABEL
%
%  Sometimes the LabelTicks and Minorticks looks unequal spaced, that occured on
%   Intervalls in units days, example: 1 5 10 15 20 25  1 5 10 ...
%                                      1   10    20     1   10 ...
%   In this case the Day 30 gets no tick !
%
%  TIMEAXIS( ... , 'zoom' ) activates a ZoomFunction, the TickLabels will
%     updated automaticly after each Zoom, requested m-file: AXE_ZOOM
%
%      Click and Drag with Left MouseKey to zoom in
%      Click with right MouseKey to zoom back
%
%      DoubleClick with left MouseKey zooms to origin
%
%  TIMEAXIS( AxeHandles , ... , 'xzoom' ) activates a ZoomFunction
%     for multiple axes along the X-Axes, the TickLabels will
%     updated automaticly after each Zoom, requested m-file: XZOOM
%
%      Click and Drag with Left MouseKey to zoom in
%      Click with right MouseKey to zoom back
%
%      DoubleClick with left MouseKey zooms to origin
%
%
%  All Inputs are optional, they can follow in any order.
%
%---------------------------------------------------------
%
% See also: DATENUM, AXE_ZOOM, XZOOM
%
%---------------------------------------------------------
%
% Examples:  
%
%% -----------------------------------------------
%  figure
%  axis([0 1000 0 1]);
%  axe = gca;
%  timeaxis(10,'zoom')
%
%  % timeaxis(axe,4,'zoom');
%
%  % timeaxis(axe,[6 32],[1999 01 01],'zoom');
%
%% -----------------------------------------------
%  figure
%  n   = 3;
%  axe = zeros(n,1);
%  for ii = 1 : n
%    axe(n-ii+1) = subplot(n,1,ii);
%    plot( 100 * ii * rand(100,1) );
%  end
%
%  timeaxis(axe,'xzoom')
%  % timeaxis(axe,'xzoom','y')
%
%  

nl = char(10);

Msg = '';

axt = [];
lbt = [];

app = mfilename;
fsp = ( app == filesep );
if any(fsp)
    app = app(max(find(fsp))+1:end);
end
app = upper(app);


%***********************************************
% Defaults

DelForm = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 );
DelForm = [ 'try, delete(' DelForm ');end' ];


ud0 = struct( 'Int'  , { 8   }     , ...
              'XY'   , { 'x' }     , ...
              'Base' , { [0 1 0 0 0 0] } , ...
              'Mode' , {  0 }      , ...
         'LabelAxe'  , { [] }      , ...
         'LabelText' , { [] }      , ...
         'Children'  , { [] }              );


[ud,axe,is_zoom,cl] = check_in(ud0,app,varargin);

%***********************************************
% Benutze folgende Intervalle

% Sekunden/MinutenIntervalle

ds = [ 1 2 5 10 15 20 30 ];

 
% StundenIntervalle

dh = [ 1  2  3  6  12 ];


% TageIntervalle

dd = [ 1  2  5  10  15 ];


% MonatsIntervalle

dm = [ 1  2  3  6  ];


% JahresIntervalle

dy = [  1  2  5  ];

dy = [ dy*1e0  dy*1e1  dy*1e2  dy*1e3  ...
       dy*1e4  dy*1e5  dy*1e6  dy*1e7  dy*1e8  dy*1e9 ];


mint = 50;  % Maximum LabelIntervall

%************************************************
% Concatinate Intervalls, units: days !!!

d0 = cat( 2 , ds , ds , dh , dd , dm , dy );

d  = cat( 2 , ds/24/3600 , ds/24/60 , dh/24 , dd , dm*30 , dy*360 );

% dy*360 !!! for LastCommonMultiple of MinorTicks, see below

%************************************************
% Find Optimal Intervall

int  = ud.Int;
base = ud.Base;
xy   = ud.XY;

base = datenum(base(1),base(2),base(3),base(4),base(5),base(6));


limp   = cat(2,xy,'lim');
tickp  = cat(2,xy,'tick');
ticklp = cat(2,xy,'ticklabel');
labp   = cat(2,xy,'label');

lim = get(axe,limp) + base;

dl = lim(2) - lim(1);

%------------------------------------------------
% Get LabelIntervall

int(1) = min( int(1) , mint );

ii = sum( d < dl/int(1) );
ii = ii + ( ii < size(d,2) );  

% [h,ii]=min( abs( dl/int(1) - d ) );

il  = d(ii);
il0 = d0(ii);

%------------------------------------------------
% Get TickIntervall 

if prod(size(int)) == 1
 
  % Automaticly

  jj = ii - 1 * ( ii > 1 );

else

  if int(2) <= int(1)
 
     jj = ii;

  else

     jj = sum( d < dl/int(2) ) + 1;

     jj = min(jj,ii);

  end

end

% Search that il is Least Common Multiple

acc = 1e-10;

while ( abs( lcmtpl(d(jj),il) - il ) > acc )  &  ( jj > 1 )
      
      jj = jj - 1;

end

% Special for 10-days LabelIntervall &  2-days TickIntervall !!!

jj = jj + ( ( il == 10 ) & ( d(jj) == 2 ) );

it  = d(jj);
it0 = d0(jj);


%------------------------------------------------
% Ticks & Labels

%------------------------------------------
if il < 1
% Stunden

  [t,tl,l,lp,lt] = hour_tick(ud.Mode,lim,it,il);

%------------------------------------------
elseif il < 360
% Tage des Monats

  [t,tl,l,lp,lt] = day_tick(ud.Mode,lim,it,il,it0,il0);

%------------------------------------------
else
% Year

  [t,tl,l] = year_tick(lim,it,il,it0,il0);

  lp = [];
  lt = [];

end

t  = t  - base;
tl = tl - base;
lp = lp - base;

%********************************************************
% Create Labels

fig      = get(axe(1),'parent');
curr_axe = get(fig,'currentaxes');

%--------------------------------------------------------
% Use DummyAxes for Labels, Check existing LabelAxe

axt = ud.LabelAxe;
 
ok = ~isempty(axt);
if ok
   ok = ishandle(axt);
   if ok
      ok = ( get(axt,'parent') == fig );
      if ~ok
         delete(axt)
      end
   end
end

if ~ok

   fcn0 = get(axe,'DeleteFcn');

   % Check old DeleteFcn

   if ~isempty(axt)
      fcn1 = sprintf(DelForm,axt);
      ii = findstr(fcn0,fcn1);
      if ~isempty(ii)
          ind = ( 1 : size(fcn1,2) );
          for jj = ii(end:1)' 
              fcn0( jj+ind-1 ) = [];
          end
      end
   end

   % Create new LabelAxe

   axt = axes( 'parent' , fig , ...
                'box'    , get(axe,'box') );

   % Expand DeleteFcn

   set( axe , 'DeleteFcn' , cat( 2 , fcn0 , ',' , sprintf(DelForm,axt) ) );  

end

%----------------------------------------------------------
% Set TickLabels now !!! (Matlab7 changes AxesPosition)
%

set( axe ,  tickp , t      , ...
           ticklp , {}     , ...
          'color' , 'none' , ... 
       'nextplot' , 'add'        );

%----------------------------------------------------------

set( axt ,   'units' , get(axe,'units')      , ...
          'position' , get(axe,'position')   , ...
          'fontname' , get(axe,'fontname')   , ...
         'fontunits' , get(axe,'fontunits')  , ...
          'fontsize' , get(axe,'fontsize')   , ...
        'fontweight' , get(axe,'fontweight') , ...
         'fontangle' , get(axe,'fontangle')  , ...
           'tickdir' , get(axe,'tickdir')    , ...
         'linewidth' , get(axe,'linewidth')  , ...
     'gridlinestyle' , get(axe,'gridlinestyle')   , ...
        'ticklength' , 2 * get(axe,'ticklength')  , ... 
        'CameraPosition'      , get(axe,'CameraPosition')      , ...
        'CameraPositionMode'  , get(axe,'CameraPositionMode')  , ...
        'CameraTarget'        , get(axe,'CameraTarget')        , ...
        'CameraTargetMode'    , get(axe,'CameraTargetMode')    , ...
        'CameraUpVector'      , get(axe,'CameraUpVector')      , ...
        'CameraUpVectorMode'  , get(axe,'CameraUpVectorMode')  , ...
        'CameraViewAngle'     , get(axe,'CameraViewAngle')     , ...
        'CameraViewAngleMode' , get(axe,'CameraViewAngleMode') , ...
              'View'          , get(axe,'View')                , ...
        'DataAspectRatio'     , get(axe,'DataAspectRatio')     , ...
        'DataAspectRatioMode' , get(axe,'DataAspectRatioMode') , ...
            'tag'    , 'TIMEAXIS' , ...
          'userdata' ,  axe       , ...
             'color' , 'none'     , ...
          'nextplot' , 'add'                 );

for c = 'xyz'

      set( axt ,  [ c 'color' ] , get(axe,[c 'color']) , ...
                  [ c 'dir'   ] , get(axe,[c 'dir'  ]) , ...
                  [ c 'grid'  ] , get(axe,[c 'grid' ]) , ...
                  [ c 'lim'   ] , get(axe,[c 'lim'  ]) , ...
                  [ c 'tick'  ] , []                   , ...
             [ c 'ticklabel'  ] , {}                         );

    if ~isequal(c,'z')

      set( axt , [ c 'axislocation' ]  , get(axe,[c 'axislocation'  ]) );

    end

end

set( axt , tickp  , tl , ...
           ticklp , l       );

%---------------------------------------------------
% Label below

% Delete old Labels

lbt = ud.LabelText;

if ~isempty(lbt)
    ii = ishandle(lbt);
    if ~any(ii)
        lbt = [];
    elseif ~all(ii)
        ii = find(ii);
        lbt = lbt(ii);
    end
end

if ~isempty(lbt)
    label_off = all( strcmp( get(lbt,'visible') , 'off' ) );
    for hh =  lbt(:)'
        delete(hh);
    end
else
    label_off = 0;
end


% Build new Labels

if isempty(lp)

    lbt = [];

else

  xl  = get( axt , labp );
  uni = get(xl,'uni');

  set( xl , 'fontname' , get(axe,'fontname')   , ...
           'fontunits' , get(axe,'fontunits')  , ...
            'fontsize' , get(axe,'fontsize')   , ...
          'fontweight' , get(axe,'fontweight') , ...
           'fontangle' , get(axe,'fontangle')         );

  nl = prod(size(lp));

  lbt = zeros(nl,1);

  np = find( xy == 'xyz' );

  lvis = { 'on'  'off' };
  lvis = lvis{ 1 + label_off };

  for ii = 1 : nl

     hh = copyobj(xl,axt);

     set( hh , 'units' , 'data' );

     pos = get(hh,'position');

     pos(np) = lp(ii);

     set( hh , 'position' , pos    , ...
               'string'   , lt{ii} , ...
               'visible'  , lvis   , ...
               'tag'      , 'TIMEAXIS_LABEL'    );

     set( hh , 'units' , uni );

     lbt(ii) = hh;

  end

end


%--------------------------------------------------------
% Set Axe

ud.LabelText = lbt;

ud.LabelAxe = axt; 

setappdata(axe,app,ud);

set( fig , 'currentaxes' , curr_axe , ...
           'toolbar'     , 'none'          );

fcn = get(fig,'resizefcn'); 
      set(fig,'resizefcn','');

drawnow

     set(fig,'resizefcn',fcn);

%-------------------------------------------
% Set FigureChildren instead of: axes(axe)

shh = get(0,'showhiddenhandles');
      set(0,'showhiddenhandles','on');

ch = get(fig,'children');

      set(0,'showhiddenhandles',shh);
      
i0 = find( ch == axe  );
i1 = find( ch == axt );

if i0 > i1

    ch([i0 i1]) = ch([i1 i0]);

    set(fig,'children',ch);

end

%-------------------------------------------
% Check for Call by AXE_ZOOM or 
%  empty Children and NO Zoom

if   strcmp(cl,'axe_zoom') | ...
   ( strcmp(is_zoom,'none') & isempty(ud.Children) )
   return
end

%*************************************************************
% Activate Zoom

 fcn = is_zoom;
 fpp = upper(fcn);

 new = 0;           % Default for new AXE_ZOOM

switch is_zoom

 %*************************************************************
 % AXE_ZOOM

 case 'axe_zoom'

  axn = cat(1,axe,ud.Children);
  nax = size(axn,1);
  new = ones(nax,1);

  %-------------------------------------------------------------
  for ii = 1 : nax
  %-------------------------------------------------------------

   %-------------------------------------------------------------
   % Check if Zoom allready exist

   hst = [];

   new(ii) = ~isappdata(axn(ii),fpp);

   if ~new(ii)

       Msg = feval(fcn,'on',axn(ii));

       new(ii) = ~isempty(Msg);  % Zoom On not successfull

       if new(ii)
          % Try to get History
          try
             hst = getappdata(axn(ii),fpp);
             hst = hst.History;
          end
       end

   end

   %-------------------------------------------------------------
   % Activate New Zoom

   if new(ii)
 
      Msg = feval(fcn,'New',axn(ii),mfilename);

      new(ii) = isempty(Msg);
 
      if new(ii)

         Msg = feval(fcn,'on',axn(ii));

         new(ii) = isempty(Msg);

         if new(ii)

            % Use old History
            if ~isempty(hst)

                zud = getappdata(axn(ii),fpp);
 
                zud.History = hst;

                setappdata(axn(ii),fpp,zud)

            end

         end

      end

   end

  %-------------------------------------------------------------
  end  % ii
  %-------------------------------------------------------------

 %*************************************************************
 % XZOOM

 case 'xzoom'

  %-------------------------------------------------------------
  % Check if Zoom allready exist

  hist = [];

  ok = isappdata(axe,fpp);

  if ok

     Msg = xzoom(axe,'Active','on');

     ok = isempty(Msg);

     if ~ok
     % Zoom On not successfull, get History
        try
           hist = getappdata(axe,fpp);
           hist = hist.XLim;
        end
     end
     
  end

  %-------------------------------------------------------------
  % Activate New Zoom

  if ~ok

     Msg = xzoom(cat(1,axe,ud.Children),'New',mfilename);

     if isempty(Msg)

        % Use old History

        if ~isempty(hist)

           zud = getappdata(axe,fpp);
 
           zud.XLim = hist;

           setappdata(axe,fpp,zud)

        end

     end

  end

end
% switch is_zoom

%*************************************************************
% Set Children

sets = { 'on'  'off' };
sets = sets{ 1 + strcmp(ud.XY,'x') };

for h = ud.Children(:)'

    [m,at,lt] = timeaxis( h , ud.Int , ud.XY , ud.Base );

    if isempty(m)

       set( lt , 'visible' , sets );

       axt = cat( 1 , axt , at );

    end

end

%-------------------------------------------------------------
% Activate new AXE_ZOOM for LabelAxe

if ~( strcmp(is_zoom,'axe_zoom') & any(new) )
    return
end

new = find(new);

for ii = new(:)'

    tud = getappdata(axn(ii),app);

    set( tud.LabelAxe , 'buttondownfcn' , get(axn(ii),'buttondownfcn') );

end

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function   [t,tl,l,lp,lt] = hour_tick(mode,lim,it,il);

lp = [];
lt = {};

mon = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';
       'Aug';'Sep';'Oct';'Nov';'Dec'];

% Ticks

n = round( il / it );

t = ( ceil(lim(1)/it) : floor(lim(2)/it) );

ok = find( mod(t,n) == 0 );

t  = it * t;


% LabelTicks

tl = t(ok);

nl = size(tl,2);

l  = cell(1,nl);

d  = datetime(tl);

%-------------------------------------------------
% Formats

%--------------------------------------------
if il < 1/24/60
% Sekunden  hh:mm:ss

      ind = [ 4  5  6 ];
     form = '%.0f:%2.2d:%2.2d';
 
%--------------------------------------------
elseif il < 1
% Stunden  hh:mm

      ind = [ 4  5 ];
     form = '%.0f:%2.2d';


end


for ii = 1 : nl

  l{ii} = sprintf(form,d(ii,ind)');

  if mode

     ld = datetime(datenum(d(ii,[1 2 3])));

     l{ii} = sprintf('%2.2d.%s.%4.4d %s',ld(3),mon(ld(2),:),ld(1),l{ii});

  end

end

%-------------------------------------------------
% DayLabel

if mode
   return
end

ld = ( floor(lim(1)) : ceil(lim(2)) )';

nd = size(ld,1);

lp = ld;

lp([1 nd]) = lim(:);

i01 = zeros(1,2);

i01(1) =  1 + ( ( lp(2)  - lp(1)    ) < il );
i01(2) = nd - ( ( lp(nd) - lp(nd-1) ) < il );

% Check if i0 <= i1 !!!

i01 = i01( [1 2] + [1 -1] * ( i01(1) > i01(2) ) );
i01 = i01 + [ -1  1 ] * ( i01(1) == i01(2) );

i01 = ( i01(1) : i01(2) );

ld = ld(i01);
lp = lp(i01);

nd = size(ld,1)-1;

lp = ( lp(1:nd) + lp(2:(nd+1)) ) / 2; 

lt = cell(nd,1);

ld = datetime(ld);

for ii = 1 : nd

  form = cat(2,'%2.2d-',mon(ld(ii,2),:),'-%4.4d');

  lt{ii} = sprintf(form,ld(ii,[3 1]));

end

if ( nl == 0 ) & ( nd == 1 ) 

       d      = datevec(lim(1));

       d(:,6) = floor(d(:,6));  % Floor seconds !!!

      ind = [ 4  5  6 ];
     form = '%.0f:%2.2d:%2.2d';

    lt{1} =  sprintf('%s %s',lt{1},sprintf(form,d(ind)));

end
 

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function   [t,tl,l,lp,lt] = day_tick(mode,lim,it,il,it0,il0);

% TageIntervalle  =  [ 1  2  5  10  15 ];
% MonatsIntervalle  =  [ 1  2  3  6  ];

lp = [];
lt = {};

mon = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';
       'Aug';'Sep';'Oct';'Nov';'Dec'];


%-------------------------------------------
% LabelTicks

t = ( ceil(lim(1)) : floor(lim(2)) )';

d = datetime(t);

mm = d(:,2);  % MonNumber
dd = d(:,3);  % DayNumber


if il < 30

   ok = find( ( mod(dd-1*(il<=2),il) == 0 ) | ( dd == 1 ) );

   if prod(size(ok)) > 1

    il_bad = min(il,3);

    bad = find( ( diff(t(ok)) < il_bad )  &  ~( dd(ok(1:(end-1))) == 1 ) );

    ok(bad) = [];

   end

else

  ok = find( ( mod(mm-1,il0) == 0 )  &  ( dd == 1 )  );

end


tl = t(ok);
dl = d(ok,[1 2 3]);

nl = size(tl,1);

 l = cell( nl , 1 );

for ii = 1 : nl

  form = cat( 2 , '%2.2d.' , mon(dl(ii,2),:) );
 
  l{ii} = sprintf(form,dl(ii,3));

  if mode

     l{ii} = sprintf('%s.%4.4d',l{ii},dl(ii,1));

  end

end
 

%-------------------------------------------
% Ticks

if it < 1

  t = it * ( ceil(lim(1)/it) : floor(lim(2)/it) )';

else

  if it < 30

    ok = find( ( mod(dd-1*(it<=2),it) == 0 ) | ( dd == 1 ) );

    if prod(size(ok)) > 1

      it_bad = min(it,3);

      bad = find( ( diff(t(ok)) < it_bad ) & ~( dd(ok(1:(end-1))) == 1 ) );

      ok(bad) = [];

    end

  else

    ok = find( ( mod(mm-1,it0) == 0 )  &  ( dd == 1 ) );

  end

  t = t(ok);
    
end


%-------------------------------------------------
% YearLabel

if mode
   return
end

ld = ( d(1,1) : ( d(end,1) + 1 ) )';

nd = size(ld,1);

lp = datenum(ld,1,1);

ls = lp(2:nd-1);
ns = nd - 2;

lp([1 nd]) = lim(:);

i0 =  1 + ( ( lp(2)  - lp(1)    ) < il );
i1 = nd - ( ( lp(nd) - lp(nd-1) ) < il );

ld = ld(i0:i1);
lp = lp(i0:i1);

nd = size(ld,1)-1;

%%% ls = lp(2:nd);   % see above

lp = ( lp(1:nd) + lp(2:(nd+1)) ) / 2; 

%%% nt = 2 * nd - 1;

nt = nd + ns;         % see above

lt = cell(nt,1);

for ii = 1 : nd

  lt{ii} = sprintf('%4.4d',ld(ii));

end

lt(nd+1:nt) = {'|'};

lp = cat( 1 , lp , ls );

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function   [t,tl,l] = year_tick(lim,it,il,it0,il0);


%-------------------------------------------
% LabelTicks

dl = datetime(lim);

tl = il0 * ( ceil(dl(1,1)/il0) : floor(dl(2,1)/il0) )';

nl = size(tl,1);

 l = cell(nl,1);

form = cat(2,char(32*ones(1,10)),'%4.4d');

for ii = 1 : nl

  l{ii} = sprintf(form,tl(ii));

end

tl = datenum(tl,1,1);


%-------------------------------------------
% Ticks

if it < 1

   t = it * ( ceil(lim(1)/it) : floor(lim(2)/it) )';

elseif it < 360

   t = ( ceil(lim(1)) : floor(lim(2)) )';

   d = datetime(t);

  mm = d(:,2);  % MonNumber
  dd = d(:,3);  % DayNumber

  if it < 30

    ok = find( ( mod(dd-1*(it<=2),it) == 0 ) | ( dd == 1 ) );

    if prod(size(ok)) > 1

      it_bad = min(it,3);

       bad = find( ( diff(t(ok)) < it_bad ) & ~( dd(ok(1:(end-1))) == 1 ) );

       ok(bad) = [];

    end

  else

    ok = find( ( mod(mm-1,it0) == 0 )  &  ( dd == 1 ) );

  end

     t = t(ok);
   
else

  t = it0 * ( ceil(dl(1,1)/it0) : floor(dl(2,1)/it0) )';

  t = datenum(t,1,1);
      
end


%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ud,axe,is_zoom,cl] = check_in(ud0,app,VarIn);

% CHECK_IN  Checks Inputs-Arguments

axe     = [];

is_zoom = 'none';

ud1 = struct( 'Int'  , { [] } , ...
              'XY'   , { [] } , ...
              'Base' , { [] } , ... 
              'Mode' , { [] } , ...
          'Children' , { [] }       );

for vv = VarIn(:)'

    val = vv{1};

    nv = prod(size(val));

    is_d = strcmp(class(val),'double');
    is_c = strcmp(class(val),'char');
   
    %----------------------------------------
    % Check for AxeHandle

    ok = ( is_d & ( nv >= 1 ) );
    if ok
       ok = all(ishandle(val));
       if ok
          ok = all( strcmp( get(val,'type') , 'axes' ) );
          if ok
             axe = val;
             if prod(size(axe)) > 1
                axe = axe(:);
                ud1.Children = axe(2:end);
                axe = axe(1);
             end
          end
       end
    end

    %----------------------------------------
    % Check for Intervall

    ok = ( ~ok & is_d &  ( nv > 0 )  &  ( nv <= 2 ) );
    if ok
       ok = all( ( mod(val,1) == 0 )  &  ( val > 0 ) );
       if ok
          ud1.Int = val(:)';
       end
    end

    %----------------------------------------
    % Check for xy

    ok = ( ~ok & is_c & ( nv == 1 ) );
    if ok
       ud1.Mode = strcmp(val,upper(val));
       val = lower(val);
       ok  = any( strcmp( val , { 'x' 'y' 'z' } ) );
       if ok
          ud1.XY = val;
       end 
    end

    %----------------------------------------
    % Check for Base 

    ok = ( ~ok & is_d  &  any( nv == [3 6] ) );
    if ok

       ok = all( mod(val,1) == 0 );

       if ok

         val = val(:)';
         val = cat( 2 , val , zeros(1,6-nv) );

         ud1.Base = val;

       end

    end

    %----------------------------------------
    % Check for ZOOM

    ok = ( ~ok & is_c & ( nv > 0 ) & ( nv == size(val,2) ) );
    if ok
       val = lower(val);
       switch val
         case  'zoom'
             is_zoom = 'axe_zoom';
         case 'xzoom'
             is_zoom = 'xzoom';
       end
    end

end

%-----------------------------------------------
% Check if TIMEAXIS allready active in axe

if isempty(axe)
   axe = gca;
end
   
%-----------------------------------------------
% Check, if requested active AxeHandle is the LabelAxes
if isequal(axe(1),gca)  &  strcmp( get(axe(1),'tag') , 'TIMEAXIS' )

   axep = get(axe(1),'userdata');

   if isnumeric(axep) & ( prod(size(axep)) == 1 )
      if ishandle(axep)
         if isappdata(axep,app)
            udp = getappdata(axep,app);
            if isequal( udp.LabelAxe , axe(1) )
               axe(1) = axep;
            end
         end
      end
   end

end

ok = isappdata(axe(1),app);
if ok
   ud0 = getappdata(axe(1),app);
end

%-----------------------------------------------
% Check for call from XZOOM or AXE_ZOOM

cl = caller;

n = min(2,size(cl,1));

cl = lower(cl{n,1});

if prod(size(ud1.Int)) == 2

   if isequal( get(axe(1),'xlim') , ud1.Int )

      if strcmp( cl , 'xzoom' )

         ud1.Int = [];

      end

   end

end    

%-----------------------------------------------
% Set Inputs in ud1 to ud0

field1 = fieldnames(ud1);

for ff = field1(:)'
    
    val = getfield( ud1 , ff{1} );
 
    if ~isempty(val);
 
        ud0 = setfield( ud0 , ff{1} , val );
   
    end

end
   
%-----------------------------------------------
% Append UserData

ud = ud0;



%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function d = datetime(t);

% DATETIME  returns correct DateTime
%
% Takes care on accuraccy-problems of DATEVEC, 
%    which returns seconds == 59.99999999999272
%

d      = datevec(t);

d(:,6) = round(d(:,6));  % Round seconds

dd = d(:,3);  % Original DayNumber

quot = [ 60 60 24 ]; %  [ ss-->mm  mm-->hh  hh-->dd ]
ind  = [ 6  5  4  ];

for ii = 1 : 3

    p = fix( d(:,ind(ii)) / quot(ii) );

 d(:,ind(ii)-0) = d(:,ind(ii)-0) - p * quot(ii);
 d(:,ind(ii)-1) = d(:,ind(ii)-1) + p;
  
end

% Check if DayNumber has changed

ii = find( d(:,3) > dd );

if isempty(ii)
   d(:,1:3) = d(:,1:3) + all( d(:,1:3) == 0  , 2 ) * [ -1 12 31 ];
   return
end

% New Date

[d(ii,1),d(ii,2),d(ii,3)] = datevec( datenum(d(ii,1),d(ii,2),d(ii,3)) );


d(:,1:3) = d(:,1:3) + all( d(:,1:3) == 0  , 2 ) * [ -1 12 31 ];

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [m,g,c,d,e,f] = lcmtpl(a,b,acc)

%LCMTPL    Least common multiple.
%
%   LCMTPL(A,B) is the least common multiple of corresponding elements of
%   A and B.  The arrays A and B must be the same size (or either can be scalar).
%
%   LCMTPL(A,B,Accuracy) gives Value of Accuracy to stop loop using GCDIV,
%                         default: EPS
%
%   [M,G,C,D] = LCMTPL(A,B) returns C and D so that G = A.*C + B.*D
%   These are useful for solving Diophantine equations and computing
%   Hermite transformations.
%
%   [M,G,C,D,E,F] = LCMTPL(A,B) returns E and F so that A = E.*G and B = F.*G
%
%   See also GCDIV.
%

%   Copyright (c) 1984-98 by The MathWorks, Inc.  &  cbegler 2001
%   $Revision: 5.6 $  $Date: 1997/11/21 23:45:41 $


if nargin < 3
   acc = eps;
end

[g,c,d,e,f] = gcdiv(a,b,acc);

m = abs(a.*f);
  

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [g,c,d,e,f] = gcdiv(a,b,acc)

%GCDIV    Greatest common divisor.
%
%   GCDIV(A,B) is the greatest common divisor of corresponding
%   elements of A and B.  The arrays A and B must be the same size 
%    (or either can be scalar).
%
%   GCDIV(0,0) is 0 by convention; all other GCDs are positive integers.
%
%   GCDIV(A,B,Accuracy) gives Value of Accuracy to stop loop, default: EPS
%
%   [G,C,D] = GCDIV(A,B) returns C and D so that G = A.*C + B.*D
%   These are useful for solving Diophantine equations and computing
%   Hermite transformations.
%
%   [G,C,D,E,F] = GCDIV(A,B) returns E and F so that A = E.*G and B = F.*G
%
%   See also LCMTPL.
%

%   Algorithm: See Knuth Volume 2, Section 4.5.2, Algorithm X.
%   Author:    John Gilbert, Xerox PARC
%   Copyright (c) 1984-98 by The MathWorks, Inc.  &  cbegler 2001
%   $Revision: 5.10 $  $Date: 1998/02/17 18:40:23 $

 
% Do scalar expansion if necessary
if length(a) == 1
   a = a(ones(size(b)));
elseif length(b) == 1
   b = b(ones(size(a)));
end

if ~isequal(size(a),size(b))
    error('Inputs must be the same size.')
end;

if nargin < 3
   acc = eps;
end

n = prod(size(a));

e = 0*a;
f = e;
c = e;
d = e;
g = e;
 
for k = 1 : n

   u = [1 0 abs(a(k))];
   v = [0 1 abs(b(k))];

   dg = acc + 1;

   while ( dg > acc )  &  v(3)

       q = round( ( u(3) - mod(u(3),v(3)) ) / v(3) );

       t = u - v*q;

       u = v;
       v = t;
 
       dg = abs( abs( a(k) / ( v(2) + (v(2)==0) ) ) - ...
                 abs( b(k) / ( v(1) + (v(1)==0) ) )       );

   end

   e(k) = abs(v(2)) * sign(a(k));
   f(k) = abs(v(1)) * sign(b(k));

   c(k) = u(1) * sign(a(k));
   d(k) = u(2) * sign(b(k));

   g(k) = ~( v(2) == 0 ) * abs( a(k) / ( v(2) + (v(2)==0) ) ) + ...
          ~( v(1) == 0 ) * abs( b(k) / ( v(1) + (v(1)==0) ) );

   g(k) = g(k) / sum( ~( v([1 2]) == 0 ) );


end

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function c = caller;

% CALLER   returns the History of the calling functions
%
%  C = CALLER;   
%
%  C = { FunctionName  LineNumer };  2-Column-CellArray
%
%  CALLER, executed from the Workspace or a function which is
%   executed directly from the Workspace, returns an empty CellArray. 
%
%  see  also:  DBSTACK
%

c = cell(0,2);   % { Name LineNr }

s = dbstack;

s = s( 3 : end );  % without caller.m and function

n = prod(size(s));

c = cell( n , 2 );

for ii = 1 : n

    [p,c{ii,1},n] = fileparts(s(ii).name);
       c{ii,2}    = s(ii).line;

end
