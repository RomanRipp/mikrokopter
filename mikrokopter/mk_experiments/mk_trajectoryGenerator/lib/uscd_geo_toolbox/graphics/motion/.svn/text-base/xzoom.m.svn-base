function Msg = xzoom(axe,action,in1,in2);

% XZOOM   ZoomFunction for unique X-Limits for multiple Axes
%
% Msg = XZOOM( Handle , Action , varargin )
%
%----------------------------------------
% 1. Initialize XZOOM
%
%  XZOOM(AxeHandles,'NEW',ZoomFcn,MiddleFcn)
%
% ZoomFcn will evaluated after successfull Zoom-Action, using FEVAL:
%
%   Msg = FEVAL( ZoomFcn{:} , AxeHandle , XLim );  % where AxeHandle == AxeHandles(1)
%
% MiddleFcn will evaluated after press the Middle MouseButton, using FEVAL:
%
%   Msg = FEVAL( ZoomFcn{:} , AxeHandle , XLim );  % where AxeHandle == AxeHandles(1)
%
%  The current XAxesLimits can''t be exeded !!!
%
%----------------------------------------
% 2. Activate XZOOM and show the XZOOM-Menu
%    ( Automaticly done if called XZOOM( ... 'New' ... ) )
%
%  XZOOM(AxeHandle,'ON'),  where AxeHandle == AxeHandles(1)
%
%----------------------------------------
% 3. DisActivate XZOOM and hide the XZOOM-Menu
%
%  XZOOM(AxeHandle,'OFF'), where AxeHandle == AxeHandles(1)
%
%----------------------------------------
% 4. Delete XZOOM
%
%  XZOOM(AxeHandle,'DELETE'), where AxeHandle == AxeHandles(1)
%
%----------------------------------------
% 5. Activate/Disactivate  XZOOM-Activity, XZOOM-Display, XZOOM-Help, XZOOM-Menu
%    Redefine the ZoomFcn or MiddleFcn
%
%  XZOOM(AxeHandle,Action,Sets), where AxeHandle == AxeHandles(1)
%
%  Action = 'Active' | 'Display' | 'Help' | 'Menu'
%  Sets   = 'on' | 'off' 
%
%  Action = 'ZoomFcn'  |  'MiddleFcn'
%  Sets   = { Fcn  varargin }
%  
%----------------------------------------------------------------   
%
% required m-files:   MHELP
%
%----------------------------------------------------------------   
%
% example:
%
%  n   = 3;
%  axe = zeros(n,1);
%  for ii = 1 : n
%    axe(ii) = subplot(n,1,ii);
%    plot( ii * rand(100,1) );
%  end
%
%  xzoom(axe,'New');
%
%  % xzoom(axe(1),'off');  % DisActivates XZOOM
%  % xzoom(axe(1),'on');   %    Activates XZOOM
%

 
Msg = '';
nl  = char(10);

app = mfilename;
fsp = ( app == filesep );
if any(fsp)
    app = app(max(find(fsp))+1:end);
end
app = upper(app);


Msg0 = sprintf('%s: ',app);

Nin = nargin;

%-----------------------------------------------------------
if Nin < 1
   axe = gca;
end

if Nin < 2
   action = 'New';
end

%-----------------------------------------------------------
% Check AxesHandles

ok = all(ishandle(axe));
if ok
   ok = all(strcmp(get(axe,'type'),'axes'));
end

if ~ok
   Msg = sprintf('%sFirst Input must be AxesHandles.',Msg0); 
  return
end

%-----------------------------------------------------------
% Check Action

if ~chkstr(action,1)
    Msg = sprintf('%sAction must be a String.',Msg0);
    return
end

%------------------------------------------------------

Msg0 = sprintf('%s(%s): ',app,action);

action = upper(action);

%------------------------------------------------------
% Check for Menu-Sub or Fcn's

switch action

   case  { 'ACTIVE'  'HELP'  'DISPLAY' }

      if Nin >= 3
         in2 = in1;
         Nin = Nin + 1;
      end

      in1 = cat( 2 , upper(action(1)) , lower(action(2:end)) );

      action = 'Menu';

   case { 'ZOOMFCN'  'MIDDLEFCN' }

      if Nin >= 3
         in2 = in1;
         Nin = Nin + 1; 
      end

      in1 = cat( 2 , upper(action(1)) , lower(action(2:(end-3))) , 'Fcn'  );

      action = 'Fcn';

end

%------------------------------------------------------
% Check Handle

if ~strcmp(action,'NEW')

    if ~isappdata(axe(1),app)
        Msg = sprintf('%sInput H must be an %s-Handle.',Msg0,app);
        return
    end

    ud = getappdata(axe(1),app);

end

%-----------------------------------------------------------

fig = get(axe(1),'parent');

switch upper(action)

%*********************************************
case 'NEW'

  if Nin < 3
    ZoomFcn = {};
  else
    ZoomFcn = in1;
  end

  if Nin < 4
    MiddleFcn = {};
  else
    MiddleFcn = in2;
  end

  %------------------------------------------------------------

axe = axe(:);


Menu = struct( 'Menu'    , { [] } , ...
               'Active'  , { [] } , ...
               'Display' , { [] } , ...
               'Help'    , { [] }         );

ud = struct('Axe'           , { axe }   , ...
            'Line'          , { 0*axe } , ...
            'Menu'          , { Menu  } , ...
            'Display'       , { []    } , ...
            'Help'          , { []    } , ...
            'XLim'          , { get(axe(1),'xlim') } , ...
            'Point1'        , { zeros(0,0)         } , ...
            'Form'          , { '%g  '             } , ...
            'DownFcn'       , { get(fig,'WindowButtonDownFcn'  ) } , ...
            'MotionFcn'     , { get(fig,'WindowButtonMotionFcn') } , ...
            'ZoomFcn'       , { {} } , ...
          'MiddleFcn'       , { {} }       ); 

     
ud.Form = [ getform(ud.XLim(end,:)) '  ' ];


%-----------------------------------------------------------------
% Lines for ZoomIntervall

for ii = 1 : size(axe,1)

  ud.Line(ii) = line('parent'    , axe(ii)       , ...
                     'xdata'     , NaN*ones(1,5) , ...
                     'ydata'     , NaN*ones(1,5) , ... 
                     'linestyle' , '-'    , ...
                     'linewidth' , 0.5    , ...
                     'marker'    , 'none' , ...
                     'erasemode' , 'xor'  , ...
                     'visible'   , 'off'  , ...
                     'tag'       , sprintf('%s_LINE',app)        );                         
end


%-----------------------------------------------------------------
% HelpMouse

   nl = char(10);
   txt = struct('Left'  ,{['1. Click to Mark first  Limit of new X-Intervall' nl  ... 
                           '2. Click to Mark second Limit of new X-Intervall' nl  ... 
                           '      or Click right to cancel first Mark '    nl nl  ...
                           'DoubleClick to Reset.'                                  ]}, ...
                'Right' ,{['Click to Zoom Out' nl ...
                           'Click after first Mark (with LeftButton) to cancel' ]} );
 
 
  % Create HelpMouse, default: 2 Buttons
   [Msg,ud.Help] = mhelp(fig,'New','ButtonNumber'    , 2         , ...
                                   'String'          , txt       , ...
                                   'backgroundcolor' , [ 230 230 200 ]/255 );


%-----------------------------------------------------------------
% Menu

  % StringFormat for Handle
   clear eps
   Hform = sprintf('%%.%.0ff',floor(abs(log(eps)/log(10))) );

     Axe = sprintf(Hform,axe(1));

   CBMenu = [ 'xzoom('  sprintf(Hform,axe(1))   ',''Menu'');' ];

  ud.Menu.Menu = uimenu( 'parent'   , fig      , ...
                         'label'    , '&XZoom'  , ...
                         'checked'  , 'off'    , ...
                         'tag'      , sprintf('%s_MENU',app)  );

  fields = fieldnames(ud.Menu);
  fields = fields(:)';

  for ff = fields(2:end)

     CB = [ 'xzoom('  Axe   ',''Menu'','''   ff{1}  ''');' ];
  
     hm = uimenu( 'parent'   , ud.Menu.Menu   , ...
                  'label'    , [ '&' ff{1} ]  , ...
                  'checked'  , 'on'   , ...
                  'callback' , CB     , ...
                  'tag'      , sprintf( '%s_MENU_%s', app , upper(ff{1}) )     );

    ud.Menu = setfield(ud.Menu,ff{1},hm);

  end


  ud.Display = uimenu( 'parent'   , fig      , ...
                       'label'    , ''       , ...
                       'checked'  , 'off'    , ...
                       'visible'  , get(ud.Menu.Display,'checked') , ...
                       'tag'      , sprintf('%s_DISPLAY',app)  );


  setappdata(axe(1),app,ud)

%-----------------------------------------------------------------
% Fcn's

  Msg1 = '';
  Msg2 = '';

  if ~isempty(ZoomFcn)
     Msg1 = xzoom(axe(1),'ZoomFcn',ZoomFcn);
  end

  if ~isempty(MiddleFcn)
     Msg2 = xzoom(axe(1),'MiddleFcn',MiddleFcn);
  end

  if ~( isempty(Msg1) & isempty(Msg2) )
     Msg = cat( 2 , Msg0 , Msg1 , nl(1:(end*(~isempty(Msg1)))) , Msg2 );
  end

%-----------------------------------------------------------------
% Activate

  xzoom(axe(1),'on');



%*********************************************
case 'MENU'

  if ~chkstr(in1,1)
     Msg = [ Msg0 'Menu-Input must be a String.' ];
     return
  end

    sets = { 'on'  'off' };

  %-------------------------------------------------
  % Check for Menu

  if any( strcmp( lower(in1) , sets ) )
     set( ud.Menu.Menu , 'visible' , in1 );
     return
  end

  %-------------------------------------------------
  % Check for Fields

  if ~any( strcmp( in1 , fieldnames(ud.Menu) ) )
     Msg = [ Msg0 'Invalid Menu-Input.' ];
     return
  end
  
  hm = getfield(ud.Menu,in1);

  %-------------------------------------------------
  % Check for Menu-Callback

  if Nin  < 4

    jj = find( strcmp( sets , get(hm,'checked') ) );

    sets = sets{3-jj};   % Switch

  %-------------------------------------------------
  else

    ok = chkstr(in2,1);
    if ok
       ok = any( strcmp( in2 , sets ) );
    end

    if ~ok 
       Msg = [ Msg0  'Last Input must be ''on'' or ''off''.' ];
       return
    end

    sets = in2;

  end

  set( hm , 'checked' , sets );

  %--------------------------------------
  switch upper(in1)

    %--------------------------------------
    case 'ACTIVE'

      xzoom(axe,sets,0);

    %--------------------------------------
    case 'DISPLAY'

      set( ud.Display , 'visible' , sets );

    %--------------------------------------
    case 'HELP'

      mhelp(ud.Help,'Set','visible',sets);

  end

%*********************************************
case 'FCN'

  if ~chkstr(in1,1)
     Msg = [ Msg0 'Fcn-Input must be a nonempty String.' ];
     return
  end

  field = in1;

  if ~any( strcmp( field , { 'ZoomFcn'  'MiddleFcn'} ) )
      Msg = [ Msg0  'Following Input must be ''ZoomFcn'' or ''MiddleFcn''.' ];
      return
  end

  %-----------------------------------------------------------

  if Nin < 4
     return
  end

  %-----------------------------------------------------------
  % Empty

  if isempty(in2)

     fcn = {};

  else

     fcn = in2;

     if ischar(fcn)
        fcn = cellstr(fcn);
     end

     ok = iscell(fcn);
     if ok
        ok = chkstr(fcn{1},1);
     end
 
     if ~ok
         Msg = [  Msg0  'Value for ' field ' must be a CharacterArray' ...
                        ' or CellArray with a String in the 1. Element.' ];
         return
     end

  end

  %-----------------------------------------------------------
  % Set

     ud = setfield( ud , field , fcn );

     setappdata(axe,app,ud);

  %-----------------------------------------------------------
  % UpDate HelpMouse

  if ~strcmp(field,'MiddleFcn')
      return
  end

  nr  = 2 + ~isempty(fcn);
  str = '';

  if ~isempty(fcn)
      str = cat( 2 , 'Call ' , upper(fcn{1}) );
  end

  str = struct( 'Middle' , {str} );

  mhelp( ud.Help , 'Set' , 'ButtonNumber' , nr , 'String' , str );

%*********************************************
case 'ON'
 
  axe = axe(1);

  % StringFormat for Handle
   clear eps
   Hform = sprintf('%%.%.0ff',floor(abs(log(eps)/log(10))) );

  CBMotion = [ 'xzoom('  sprintf(Hform,axe)   ',''Move'' );' ud.MotionFcn ];
  CBDown   = [ 'xzoom('  sprintf(Hform,axe)   ',''Down'' );' ];

  CBDown0  = get(fig,'WindowButtonDownFcn');

  if ~isequal(CBDown,CBDown0)

     ud.DownFcn = get(fig,'WindowButtonDownFcn');

  end

set(fig,'WindowButtonDownFcn'    , CBDown         , ...
        'WindowButtonMotionFcn'  , CBMotion               );


xl = get(axe,'xlim');

if ~isequal(ud.XLim(end,:),xl)

  ud.XLim = cat(1,ud.XLim,xl);

end

  set(ud.Axe,'xlim',ud.XLim(end,:));

  ud.Form = [ getform(ud.XLim(end,:)) '  ' ];

  setappdata(axe,app,ud);

  setline(ud.Axe,ud.Line,mean(ud.XLim(end,:)),NaN);


    set( ud.Line         , 'visible' , 'on' );
    set( ud.Menu.Menu    , 'visible' , 'on' );

    set( get(ud.Menu.Menu,'children') , 'visible' , 'on' , ...
                                        'enable'  , 'on'        );

    set( ud.Menu.Active , 'checked' , 'on' );

    set( ud.Display      , 'visible' , get(ud.Menu.Display,'checked') );
  mhelp( ud.Help , 'Set' , 'visible' , get(ud.Menu.Help   ,'checked') );


%*********************************************
case 'OFF'

 if nargin < 3
   in1 = 1;     % Hide the XZOOM-Menu
 end

  set(fig,'WindowButtonDownFcn'   , ud.DownFcn    , ...
          'WindowButtonMotionFcn' , ud.MotionFcn          );

    set( ud.Line    ,      'visible' , 'off' );
    set( ud.Display ,      'visible' , 'off' );
  mhelp( ud.Help , 'Set' , 'visible' , 'off' );

    set( [ ud.Menu.Display ud.Menu.Help ] , 'enable' , 'off' );

    set( ud.Menu.Active , 'checked' , 'off' );

 if in1 
   set(ud.Menu.Menu,'visible','off');
 end


%*********************************************
case 'DELETE'

  set(fig,'WindowButtonDownFcn'   , ud.DownFcn    , ...
          'WindowButtonMotionFcn' , ud.MotionFcn          );

  delete(ud.Line)
 
  delete(ud.Menu.Menu)

  delete(ud.Help)

  rmappdata(axe(1),app);
   
%*********************************************
case 'MOVE'

  cp = get(axe,'currentpoint');

  cp = cp(1,[1 2]);

  if  ( cp(1) < ud.XLim(end,1) )  |  ( cp(1)  > ud.XLim(end,2) )
    return
  end

  setline(ud.Axe,ud.Line,cp(1));

  set( ud.Display,'label',sprintf(ud.Form,[ ud.Point1  cp(1) ]));


%*********************************************
case 'DOWN'

 h = get(fig,'currentobject');

 if ~strcmp(get(h,'type'),'axes')
    h = get(h,'parent');
    if ~strcmp(get(h,'type'),'axes')
        return
    end  
 end

 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 % Check for LabelAxes of TIMEAXIS
 if ~any( h == ud.Axe )
    if strcmp( get(h,'tag') , 'TIMEAXIS' )
       h = get(h,'userdata');
    end
 end
 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if ~any( h == ud.Axe )
    return
 end

 select = get(fig,'selectiontype');

 switch upper(select) 

   %-------------------------------------------
   case 'NORMAL'

     cp = get(h,'currentpoint');
     cp = cp(1,[1 2]);
         
     if isempty(ud.Point1)
     % Store Point as FirstPoint
     % Fix Line

         ud.Point1 = cp(1);

         setappdata(axe,app,ud);

         setline(ud.Axe,ud.Line,cp(1),cp(1));

         return

     end
 
     xl = cat(2,ud.Point1,cp(1));
     xl = xl( ( [ 1  2 ] + [ 1 -1 ] * ( xl(1) > xl(2) ) ) );

     ud.Point1 = [];

     if xl(1) < xl(2)

       ud.XLim = cat( 1 , ud.XLim , xl );

     end

   %-------------------------------------------
   case 'ALT'

     if ~isempty(ud.Point1)
     % Cancel 1. Point
     % Delete fixed Line

       setline(ud.Axe,ud.Line,NaN,NaN);

       ud.Point1 = [];

       setappdata(axe,app,ud);

       return

     end


     if size(ud.XLim,1) <= 1
        return
     end

     % Zoom 1 Time up
     ud.XLim = ud.XLim(1:end-1,:);


   %-------------------------------------------
   case 'OPEN'
   % Zooms to first
    
     if size(ud.XLim,1) <= 1
         return
     end

     ud.XLim = ud.XLim(1,:);

     ud.Point1 = [];

   %-------------------------------------------
   case 'EXTEND'

      if isempty(ud.MiddleFcn)
         return
      end

      MsgF = '';

      try
        MsgF = feval( ud.MiddleFcn{:} , axe , ud.XLim(end,:) );
      catch
        MsgF = lasterr;
      end
 
      if ~isempty(MsgF)

         MsgF = [ Msg0 ' Error in UserZoomFcn. ' nl ...
                         MsgF nl ];

         fprintf([ MsgF  nl ]);

      end

  end
  %-------------------------------------------
 
     % Save YLimMode auto
     set(ud.Line,'visible','off');
   
     setline(ud.Axe,ud.Line,NaN,NaN);

     set(ud.Axe,'xlim',ud.XLim(end,:));

%     drawnow

     ud.Form = [ getform(ud.XLim(end,:)) '  ' ];

     setappdata(axe,app,ud);

     xzoom(ud.Axe(1),'Move');

     set(ud.Line,'visible','on');

  %------------------------------------------------

  if isempty(ud.ZoomFcn)
     return
  end

  MsgF = '';

  try
    MsgF = feval( ud.ZoomFcn{:} , axe , ud.XLim(end,:) );
  catch
    MsgF = lasterr;
  end
 
  if ~isempty(MsgF)

     MsgF = [ Msg0 ' Error in UserZoomFcn. ' nl ...
                     MsgF nl ];

     fprintf([ MsgF  nl ]);

  end

end



%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function setline(axe,line,x1,x0);

Nin = nargin;

emodes  = { 'xor'  'background' };

 for ii = 1 : size(axe,1)

    yl = cat(2,get(axe(ii),'ylim'),NaN,NaN);

    xx = get(line(ii),'xdata');
    if Nin == 4
      xx([1 2]) = x0;
    end

%    if ~isnan(x1)
      xx([4 5]) = x1;
%    end

    yy = get(line(ii),'ydata');
    yy([1 2 4 5]) = yl([1 2 1 2]+2*isnan(xx([1 2 4 5])));

    set(line(ii),'xdata',xx, ...
                 'ydata',yy, ...
                 'erasemode' , emodes{1+(xx(1)==xx(4))}      );
%  disp([ xx ; yy ])

 end

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function form = getform(lim);


nk = floor( log( (lim(2)-lim(1))/1e3 ) / log(10) + 2*eps );

vk = floor( log(max(abs(lim))) / log(10) + 2*eps ) + 1;

form = sprintf('%%.%.0fg',vk-nk);


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

