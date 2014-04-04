function [Msg,hm,hf,ht] = msg_logo(fig,action,txt,varargin);

% MSG_LOGO  Creates a Logo and MessageListBox in top of figure,
%             using  MSG_LIST  and  TEXTLOGO
%
%-----------------------------------------------------------------
% 1.  Create a Logo and MessageListBox
%
% [ Msg, MessageListHandle, LogoFrameHandle, LogoTextHandle] = ...
%    MSG_LOGO( FigureHandle , 'New' , LogoText , ...
%               LogoForeGroundColor , LogoBackGroundColor )
%
% FigureHandle  specifies the Figure, in which the Objects will
%                created
%
% LogoText      CharacterArray or CellStringArray for Logo, containing N Strings
%                 like: { 'O'  'OceanData'  'ToolBox' }
%
%
% [ Msg, MessageListHandle, LogoFrameHandle, LogoTextHandle] = ...
%    MSG_LOGO( FigureHandle , 'New' , LogoText , ...
%               Property1 , Value1 , Property2 , Value2 , ... )
%
% Use specified Properties for MessagList and Logo:
%
%  MessageList (see MSG_LIST):
%
%          'PixelOffset'         [ Left Right Top ]
%          'RowNumber'           
%          'RecentNumber'        
%          'AppendString'        
%          'NumberForm'          
%          'FontSize'            
%          'FontName'            
%          'FontWeight'          
%          'ForeGroundColor'  |  'MsgFGC'
%          'BackGroundColor'  |  'MsgBGC'
%
%  Logo (see TEXTLOGO):
%
%          'Offset'             
%          'FAspect'            
%          'HAspect'            
%          'FName'              
%          'FWeight'            
%          'FAngle'             
%          'FontSize'           
%          'FontName'           
%          'FontWeight'         
%          'ForeGroundColor'  |  'LogoFGC' 
%          'BackGroundColor'  |  'LogoBGC'    
%
% More Inputs can be UIControl-Properties-PropertyValue-Pairs.
%
%------------------------------------------------------------------------
% Output
%
% Msg                contains ErrorMessages, empty if all ok
%
% MessageListHandle  Handle, refering to MessageListBox, created by
%                     MSG_LIST
% 
% LogoFrameHandle    Handle of LogoFrameUIControl, created by TEXTLOGO
% LogoTextHandle     Handle of  LogoTextUIControl, created by TEXTLOGO
%
%
% Note:  
%
%  The Hight of the Logo is the same like the Height of the MessageList,
%    the Width of the Logo will fitted to the Width of the
%      maxium Length of the Text
%  DefaultSettings for the MessageList and the Logo are defined  
%   at the BEGIN of this File, that make's sure, that the ToolBox's, 
%   which use MSG_LOGO, have an unique style.
%  Edit this File to change the Settings.
%
%  The     Units of the UIControl's are set to 'pixels'.
%  The FontUnits of the UIControl's are set to 'points'.
%
%  To use the MessageList, use MSG_LIST with the Input MessageListHandle.
%
%  If the figure is resized, use following command to Resize the Logo
%    and MessageList.
%
%
%-----------------------------------------------------------------
% 2. Resize the Logo and MessageList after resizing the figure
%
%  Msg = MSG_LOGO( LogoFrameHandle , 'Resize' , ResizeMode );
%
%  ResizeMode             'normalized'  or  'absolut' 
%                 short:  'n'               'a'
%               default:  'normalized'
%
% See also MSG_LIST.
% This command could be added to the ResizeFunction of the Figure.
%
%
%-----------------------------------------------------------------
% example:
%
%  % Create a figure
%   fig = figure('position',[100 200 400 300]);
%
%  txt = { 'O'  'OceanData'        'ToolBox' };  % 3-Line Text
%
%  % txt = { 'E'  'ExperimentSetup'  'ToolBox' };  % Try this too
%
%  [Msg,HList,HLogo,HText] = msg_logo(fig,'New',txt);
%
%  % Check the ErrorMessage
%  if ~isempty(Msg)
%     fprintf([ 'Error using MSG_LOGO: '  Msg ]);
%     return
%  end
%
%  % New Message
%  Msg = msg_list( HList, 'Message', 'Read Data' );
% 
%  % Append a Message
%  Msg = msg_list( HList, 'Message', 'ok', 'append' );
%
%  % Resize the figure and the Logo, MessageList
%  %  Executes following commands step by step !!!
%
%    set(fig,'position',[100 200 500 500]);
%  % STOP Look what's happen
%
%    Msg = msg_logo(HLogo,'Resize','absolut');
%  % STOP Look what's happen
%
%    Msg = msg_logo(HLogo,'Resize','normalized');
%  % STOP Look what's happen
%   
%  % or set the ResizeFunction of the Figure
%
%    % StringFormat for Handle
%     clear eps
%     Hform = sprintf( '%%.%.0fg' , ceil(abs(log(eps)/log(10)))+1 );
%   
%   ResizeFcn = ['msg_logo(' sprintf(Hform,HLogo) ',''Resize'',''absolut'');' ];
%
%   set(fig,'ResizeFcn',ResizeFcn); 
%
%   % Test it by Resize the Figure using the Mouse !!!
%
%



Msg = '';
hm = [];
hf = [];
ht = [];

nl = char(10);

Msg0 = 'MSG_LOGO: ';

nm0 = size(Msg0,2);

nl0 = char([ 10 32*ones(1,nm0+0) ]);

mode = 'absolut';   % default ResizeMode


%*********************************************************
% Check Inputs

Nin = nargin;

if Nin < 2
  Msg = [ Msg0  'Not enough InputArguments.' ];
  return
end

ok = ( isnumeric(fig) & ( prod(size(fig)) == 1 ) );
if ok
  ok = ishandle(fig);
  if ok
     ok = any( strcmp( get(fig,'type') , { 'figure' 'uicontrol' } ) );
  end
end

if ~ok
 Msg = [  Msg  nl0(1:(end*(~isempty(Msg)))) ...
         'First Input must be a Handle of a Figure or UIControl.' ];
end

ok = ( ischar(action)  & ( prod(size(action)) == size(action,2) ) );
if ok
  action = upper(action);
  ok = any( strcmp( action , { 'NEW' 'RESIZE' 'VISIBLE' } ) );
end

if ~ok
 Msg = [  Msg  nl0(1:(end*(~isempty(Msg)))) ...
         'Second Input must be a String: ''New'', ''Resize'' or ''Visible''.' ];
end

if ~isempty(Msg)
  Msg = [ Msg0   Msg ];
  return
end 

%--------------------------------------------------------------
% Check FrameHandle

if any( strcmp( action , {'RESIZE' 'VISIBLE'} ) )

  field0 = { 'Handle'
             'Offset'
             'FramePosition'
             'TextPosition'
             'FontSize'      };
  hf = fig;

  ok = strcmp(get(hf,'type'),'uicontrol');
  if ok
    ok = strcmp(get(hf,'style'),'frame');
    if ok
      ud = get(hf,'userdata');
      ok = isstruct(ud);
      if ok
        ok = isequal( fieldnames(ud) , field0 );
      end
    end
  end

  if ~ok
    Msg = [ 'First Input must be a valid Handle of LogoFrame, created by MSG_LOGO.' ];
  end
    
end

%--------------------------------------------------------------

switch action

%******************************************************************
case 'NEW'
  
 if ~strcmp(get(fig,'type'),'figure')
   Msg = [  Msg  nl0(1:(end*(~isempty(Msg)))) ...
           'First Input must be a FigureHandle.' ];
 end

 if Nin < 3

   Msg = [ Msg0  'Input LogoText is undefined.' ];

   return

 end

 %----------------------------------------------------------- 
 % LogoText
 
 notxt = isempty(txt);

 ok = notxt;
 if ~ok
    [ok,txt] = chkcstr(txt,0);
 else
     txt = {''};
 end

 if ~ok
     Msg = [  Msg0  'LogoText must be a nonempty CharArray or CellStringArray.' ];
     return
 end

 txt = txt(:);

 nt = size(txt,1);

 [m,MsgIn,LogoIn,prop] = checkin(nt,varargin);

 if ~isempty(m)
     frm = sprintf('%s%s','%s',nl0);
     Msg = [ Msg  nl0(1:(end*(~isempty(Msg)))) ...
             sprintf(frm,m{:}) ];
 end

 if ~isempty(Msg)
     Msg = [ Msg0   Msg ];
     return
 end 

%******************************************************************

     MsgOffset = MsgIn{1,2};                    % [ Left  Right Top ]
    LogoOffset = cat( 2 , MsgOffset(1) , 0 );   % [ Left  Right ]

   %******************************************************************
  % MessageListBox in Top

    MsgIn = permute(MsgIn,[2 1]);

    [Msg,hm,pos] = msg_list( fig , 'New' , MsgIn{:} , prop{:} );

    if ~isempty(Msg)
      Msg = [ Msg0 'Error using MSG_LIST( New ).' nl Msg ];
      return
    end

  %******************************************************************
  % Logo in upper left Corner, same Height like MessageList

    LogoIn = permute(LogoIn,[2 1]);

    % Width == NaN   ==>  automaticly fitted  
    pos = [ LogoOffset(1)+1   pos(2)   NaN   pos(4) ];

    [Msg,hf,ht] = textlogo( fig , txt , pos , LogoIn{:} , prop{:} );
 
    if ~isempty(Msg)
      Msg = [ Msg0 'Error using TEXTLOGO.' nl Msg ];
      return
    end

    LogoPos = get(hf,'position');

  %------------------------------------------------------------
  % First LogoText in BOLD

    hh = struct( 'List'  , { hm } , ...
                 'Frame' , { hf } , ...
                 'Text'  , { ht }       );

    offs = struct( 'Msg' , {  MsgOffset } , ...
                  'Logo' , { LogoOffset }        );

    % Store original  Positions and FontSize in UserData,
    %  need them for Resize,
    %   set correct FieldNames for Check under RESIZE !!!
   
    ud = struct( 'Handle'        , { hh }       , ...
                 'Offset'        , { offs }     , ...
                 'FramePosition' , { LogoPos  } , ...
                  'TextPosition' , { get(ht,'position')  }  , ...
                  'FontSize'     , { get(ht,'fontsize')  }        );

    set(hf,'userdata',ud);

    if notxt
       set([hf ht],'visible','off');
    end

  %******************************************************************
  % Shift the MessageList in Left to Make space for Logo
  % Modify the UserData in PixelOffset(1) 

    mud = get(hm,'userdata');

    if notxt
       mud.PixelOffset(1) = MsgOffset(1);
    else
       mud.PixelOffset(1) = LogoOffset(1) + LogoPos(3) + LogoOffset(2) + ...
                            MsgOffset(1) ;
    end

    set(hm,'userdata',mud);

  %------------------------------------------------------------------
  % Resize the MessageList with modified Userdata
 
    [ Msg , pos ] = msg_list(hm,'Resize','absolut');

  %------------------------------------------------------------------
  % UpDate NormalizedPosition of MessageList

    mud              = get(hm,'userdata');
    mud.NormPosition = mud.PixelPosition ./ mud.FigPos0([3 4 3 4]);

    set(hm,'userdata',mud);
  


%******************************************************************
case 'RESIZE'

 if Nin == 3
   mode = txt;
   ok = ( ischar(mode)  &  ~isempty(mode)  &  ...
          ( prod(size(mode)) == size(mode,2) ) );
   if ok
    mode = lower(mode(1));
    ok = any( strcmp( mode , { 'n' 'a' } ) );
   end

   if ~ok
     Msg = [  Msg  nl0(1:(end*(~isempty(Msg)))) ...
              'ResizeMode must be a String, begin with ''n'' or ''a''.' ];
   end
 end

 if ~isempty(Msg)
   Msg = [ Msg0  Msg ];
   return
 end
  


  MsgOffset = ud.Offset.Msg;
 LogoOffset = ud.Offset.Logo;
 LogoPos    = ud.FramePosition;


  fig = get(hf,'parent');

  %------------------------------------------------------------------
  % FigurePosition in Pixels

    figuni = get(fig,'units');
             set(fig,'units','pixels');
    figpos = get(fig,'position');
             set(fig,'units',figuni);

  %------------------------------------------------------------------
  % Resize the MessageList
 
    [ Msg , pos ] = msg_list(ud.Handle.List,'Resize',mode);

    if ~isempty(Msg)
      Msg = [ Msg0 'Error using MSG_LIST( New ).' nl Msg ];
      return
    end

  %------------------------------------------------------------------
  % Make space for Logo

    mud = get(ud.Handle.List,'userdata');
 
    LogoWidth = LogoPos(3) * pos(4) / LogoPos(4);
   
    if strcmp(mode,'n')
       scw = figpos(3) / mud.FigPos0(3);  % Scaling for Width
       MsgOffset([1 2]) =  MsgOffset([1 2]) * scw;
      LogoOffset([1 2]) = LogoOffset([1 2]) * scw;
    end
 
    % Left Distance of MessageListBox
    llm  = MsgOffset(1) + LogoWidth + LogoOffset(1) + LogoOffset(2) + 1;

    % Correct Position of MessageList
    pos(3) = pos(3) + ( pos(1) - llm );
    pos(1) = llm;

    pos([3 4]) = pos([3 4]) + ( 1+2*mud.BorderWidth - pos([3 4]) ) .* ...
                              ( 1+2*mud.BorderWidth > pos([3 4]) );

    set( mud.ListHandle , 'units'    , 'pixels' , ...
                          'Position' , ( pos + [ 1  1   -2  -2 ] * mud.BorderWidth ) );

    set( mud.AxeHandle , 'units'      , 'pixels' , ...
                         'position'   ,  pos     , ...
                         'visible'    , 'off'          );

    set( mud.FrameHandle , 'units'      , 'pixels' , ...
                           'position'   ,  pos           );
 
  %------------------------------------------------------------------
  % New LogoPosition

    pos = [ LogoOffset(1)+1  pos(2)  LogoWidth  pos(4) ];

    % Scaling of Position and FontSize
    sc = pos(4) ./ LogoPos(4);

    pos([3 4]) = pos([3 4]) + ( [ 1  1 ] - pos([3 4]) ) .* ...
                              ( pos([3 4]) <  1 );

    set( ud.Handle.Frame , 'units'    , 'pixels' , ...
                           'position' ,  pos            ); 
 
    for ii = 1 : size(ud.Handle.Text)

       pos1        = ud.TextPosition{ii};
       pos1([1 2]) = pos1([1 2]) - ud.FramePosition([1 2]);
       pos1        = pos1 * sc;
       pos1([1 2]) = pos1([1 2]) + pos([1 2]);

       pos1([3 4]) = pos1([3 4]) + ( [ 1  1 ] - pos1([3 4]) ) .* ...
                                   ( pos1([3 4]) <  1 );
  
       set( ud.Handle.Text(ii) , ...
          'units'     , 'pixels' , ...
          'position'  ,  pos1     , ...
          'fontunits' , 'points' , ...
          'fontsize'  , floor( ud.FontSize{ii} * sc ) );

    end

%*********************************************************
case 'VISIBLE'

  ok = ( Nin >= 3 );

  if ok

    sets = txt;

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
   

   ud = get(hf,'userdata');

   set([ ud.Handle.Frame ; ud.Handle.Text ] , 'visible' , sets  );

   Msg = msg_list(ud.Handle.List,'visible',sets);


end
% switch

%***********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,MsgIn,LogoIn,vin] = checkin(n,vin)

msg = cell(0,1);


PixelOffset = [ 3  3  3 ];   % PixelOffset:  [ Left Right Top ]

if n == 1
   cf = 18; 
else
   cf = polyfit(linspace(1,n,3),[18 12 10],2);
end

LogoFont  = round(polyval(cf,(1:n)));  % FontSize for N Strings

LogoWeight    = { 'normal' };
LogoWeight    = LogoWeight(ones(1,n));
LogoWeight(1) = {'bold'};

LogoHeight    = ones(1,n);
LogoHeight(1) = 1.1;

% Determine RowNumber depending on ScreenHeight
% >= [  640 x  480 ]  ... 3 Rows
% >= [  800 x  600 ]  ... 4 Rows
% >= [ 1240 x 1024 ]  ... 5 Rows

uni = get(0,'units');
      set(0,'units','pixels');
ssi = get(0,'screensize');
      set(0,'units',uni');

RowNumber = 3 + ( ssi(4) >= 600 );


% Inputs for MessageListBox  and  Logo

 MsgIn = { 'PixelOffset'         PixelOffset     
           'RowNumber'           RowNumber             
           'RecentNumber'        100              
           'AppendString'        ' ... '         
           'NumberForm'          '%4.3d: '
           'FontSize'            8-strcmp(computer,'PCWIN')             
           'FontName'            'helvetica'       
           'FontWeight'          'normal'        
           'ForeGroundColor'     [ 0  0  0   ]     
           'BackGroundColor'     [ 1  1  0.8 ]   };

LogoIn = { 'Offset'              [ 0  0 ]            
           'FAspect'             ( LogoFont / min(LogoFont) ) 
           'HAspect'             LogoHeight   
           'FName'               'helvetica'         
           'FWeight'             LogoWeight
           'FAngle'              'normal'    
           'FontSize'            LogoFont(end)       
           'FontName'            'helvetica'       
           'FontWeight'          'normal'        
           'ForeGroundColor'     [ 0     0.3    1.0 ]     
           'BackGroundColor'     [ 0.95  0.99   1.0 ]   };


if isempty(vin)
   return
end

nin = prod(size(vin));

%------------------------------------------------------------
% Old Syntax: msg_logo(fig,'New',txt,LogoFGC,LogoBGC)

if isnumeric(vin{1}) 
   LogoIn{end-1,2} = vin{1};
   if nin == 1
      vin = {};
      return
   end
   vin = vin(2:nin);
   nin = nin - 1;
   if isnumeric(vin{1})
       LogoIn{end,2} = vin{1};
       if nin == 1
          vin = {};
          return
       end
       vin = vin(2:nin);
       nin = nin - 1;
   end
end

%%% vin

%------------------------------------------------------------
        
ok = ( mod(nin,2) == 0 );
if ok
   nin = nin/2;
   vin = reshape(vin,2,nin);
   ok = chkcstr(vin(1,:));
end

if ~ok
    msg = cat(1,msg,{'Additional Inputs must be Property-Value-Pairs.'}, ...
                    {'Properties must be Strings.'} );
    return
end

%------------------------------------------------------------

ok = zeros(1,nin);

mopt = lower( MsgIn(:,1));
lopt = lower(LogoIn(:,1));

fini = { 0  'fname'    'fontname' 
         0  'fweight'  'fontweight'
         0  'fangle'   'fontangle'  };

for ii = 1 : nin

    p = lower(vin{1,ii});
    v = vin{2,ii};

    im = strcmp( mopt , p );
    il = strcmp( lopt , p );

    ok(ii) = 1*any(im) + 2*any(il);

    for jj = 1 : size(fini,1)
        fini{jj,1} = ( fini{jj,1} | strcmp(p,fini{jj,2}) );
        if strcmp(p,fini{jj,3}) & ~fini{jj,1}
                  kk    = find( strcmp(lopt,fini{jj,2}) );
           LogoIn{kk,2} = v;
        end
    end

    if any( ok(ii) == [ 1  3 ] )
             im    = find(im);
       MsgIn{im,2} = v;
    end

    if any( ok(ii) == [ 2  3 ] )
              il   = find(il);
       LogoIn{il,2} = v;
    end

    if ok(ii) == 0

        switch p

          case { 'fore'  'back' }

               isf = strcmp(p,'fore');
               MsgIn{end-isf,2} = v;
              LogoIn{end-isf,2} = v;

              ok(ii) = 3;
 
          case { 'msgfgc' 'msgbgc' }

               isf = strcmp(p,'msgfgc');
               MsgIn{end-isf,2} = v;

               ok(ii) = 1;

          case { 'logofgc' 'logobgc' }

               isf = strcmp(p,'logofgc');
               LogoIn{end-isf,2} = v;

               ok(ii) = 2;

        end
        % p
    end
    % ok(ii)
end
% ii

ok = ( ok == 0 );

if ~any(ok)
    vin = {};
else
    ok  = find(ok);
    vin = vin(:,ok);
end

%------------------------------------------------------------
% Check PixelOffset

v = MsgIn{1,2};

nv = prod(size(v));

ok = isnumeric(v) & any( nv == [ 1  3 ] );
if ~ok
    msg = cat(1,msg,{'PixelOffset must be a single or 3-Element numeric.'});
elseif nv == 1
    MsgIn{1,2} = v(ones(1,3));
end


%***********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

