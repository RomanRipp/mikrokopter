function  [Msg,out1,out2] = mhelp(h,action,varargin);

% MHELP  shows a Mouse in Figure, containing HelpStrings in ToolTips
%
%
% [ Msg, ... ] = MHELP( Handle , Action , ... );
%
%  Note:  The first Output contains ErrorMessages for invalid Inputs.
%         If all is ok, Msg is empty.
%
%-------------------------------------------------------------------------
% 1. Create a new HelpMouse in a Figure, Action 'NEW'
%
%   
% [ Msg, MouseHandle, PropertyList ] = MHELP( ...
%        FigureHandle, 'New', Property1, PropertyValue1, ... );
%
%
% Creates a new HelpMouse in the Figure, defined by 
%   FigureHandle, build from PushButtons and a Frame arround,
%   and returns the Handle (Frame) of the HelpMouse and 
%   their PropertyList (UserData of MouseHandle). 
%   The Units for the Position of the HelpMouse 
%   in the Figure are set to 'pixels'.
% 
% Additional Inputs are Properties and their Values:
%
%    'ButtonNumber' , [ Number of MouseButtons: 2 or 3 ]
%    'Title'        , HelpTitleStructure  for the MouseButtons
%    'String'       , HelpStringStructure for the MouseButons
%    'CData'        , [ TrueColorMatrice ]  or  'default'  or  FileName
%    'CallBack'     , CallBackStructure for PushButtons
%    'Position'     , [ Left  Bottom ]  , PixelPosition in Figure
%    'PixelOffset'  , [ Width Height ]  , PixelOffset arround CData
%
% You can define a 2-Button or 3-ButtonMouse,
%
% The HelpTitle- and HelpStringStructure could contains the fields:
%
%   'Left'   : HelpForLeft
%   'Middle' : HelpForMiddle
%   'Right'  : HelpForRight
%
% Multiple Rows in the Strings must be separated by NewLineCharacters (char(10))
% When the Pointer is over the Buttons, the ToolTipString will shown,
%  build like:
%
%   <HelpTitle>
%   -----------
%   <HelpString> 
%
% The CallBackStructure could contains the fields:
%
%   'Center' : CallbackForCenterButton
%   'Left'   : CallbackForLeftButton
%   'Middle' : CallbackForMiddleButton
%   'Right'  : CallbackForRightButton
%
% When this Buttons are pressed, the Callbacks will evaluated.
%
% The CData will set to the CenterButton by the CDataProperty for a 
%  PushButton. The Value for CData could be:
%    - a 3-dimensional TrueColorImage, with 0 <= CData <= 1  
%       NaN's which will filled with the BackGroundColor,
%    - the string 'default', then the DefaultCData will used
%    - a FileName for an XPM-ImageFile, the Data will read from the 
%       File using READXPM !!!
%
%
% More Inputs are UIControlProperties and their Values, like
%  'BackGroundColor', 'Visible' , ...
% Following UIControlProperties are not allowed:
%   'Style', 'Tag', 'UserData'
%   'ToolTipString', 'ButtonDownFcn', 'DeleteFcn'
% (and the Properties 'String', 'CallBack', 'Position' are occupied
%  by the MouseHandleProperties)
%
% Note: The Position of the HelpMouse in Width and Height will defined
%       automaticly by the Size of the Image
%
%-------------------------------------------------------------------------
% 2. Change the HelpMouseProperties, Action 'SET'
%
% [ Msg, PropertyList ] = MHELP( ...
%        MouseHandle, 'Set', Property1, PropertyValue1, ... );
%
% All HelpMouseProperties, you could define when you create the HelpMouse,
%  can be changed later using Action 'SET'.
% May be to change the ButtonNumber, Position, CData, HelpString,
%  or the UIControlProperties like BackGroundColor, Visible
% 
% Note: Please use only MHELP to change the UIControlProperties,
%        that makes sure, that all works correct.
%        (MHELP changes the Properties of 5 UIControls)
%      Never change the Property 'Tag' of the UIControls.
%
% 
%-------------------------------------------------------------------------
% 3. Tips
%
% If you Press with the RightMouseButton on the Center, 
%  you can move the HelpMouse in the Figure.
%
% If you have resized the Figure ( New FigureWidth and FigureHeight )
%   make sure, that the HelpMouse is indide the Figure, call:
%   >> mhelp(MouseHandle,'Resize')
%
% To hide the HelpMouse, use 
%   >> mhelp(MouseHandle,'Set','visible','off')  
%
% To delete the HelpMouse, just use:
%   >> delete(MouseHandle)
%
%-------------------------------------------------------------------------
% 5. Example, execute the commands Step by Step
%
% % HelpText
%  nl = char(10);
%  txt = struct('Left'  ,{['Press left Button to:' nl ' ... ']}, ...
%               'Right' ,{['Press right Button and ' nl ' ... ' nl 'etc' ]} );
%
%  fig = figure;
%
% % Create HelpMouse, default: 2 Buttons
%  [Msg,HMouse,ud] = mhelp(fig,'New','String',txt,...
%                              'backgroundcolor',[240 240 220]/255); 
% 
% if ~isempty(Msg)
%   % An Error occured
%   fprintf([ Msg  char(10) ]);
%   return
% end
%
% % Switch to 3-ButtonMouse
%  [Msg,ud]  = mhelp(HMouse,'Set','ButtonNumber',3);
%
% % We've forgotten to set the HelpString for the MiddleButton
%  [Msg,ud]  = mhelp(HMouse,'Set','String', ...
%                     struct('Middle',{['Nothing to do']})   );
%
% % Change the Title and ImageData of the HelpMouse
%  titel = struct('Left'  ,{['Button1']}, ...
%                 'Right' ,{['']} );
%
% % 30 by 30 CData, 10 Pixel oversized (5 each Side)
%  [Msg,ud]  = mhelp(HMouse,'Set','Title', titel , ...
%                                 'CData', rand(30,30,3), ...
%                                 'PixelOffset',[10 10]       );
%
% % Back to default
%  [Msg,ud]  = mhelp(HMouse,'Set','CData','default');
%
% % The best image in town ...
%  [Msg,ud]  = mhelp(HMouse,'Set','CData','mouse.xpm');
%
%    
%   
%

Msg  = '';
out1 = [];
out2 = [];

nl = char(10);


Msg0 = ' MHELP: ';

Nin = nargin;

%-----------------------------------------------------
% Check Basic Inputs

if Nin < 2
  Msg = [ Msg0  'Inputs H and Action are undefined.' ];
  return
end

ok = ( isnumeric(h)  &  ( prod(size(h)) == 1 ) );
if ok
 ok = ishandle(h);
end

if ~ok
 Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
         ' First Input must be a single Handle.' ];
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

%-----------------------------------------------------
% Check Handle

if strcmp(action,'NEW')

 if ~strcmp(get(h,'type'),'figure')
   Msg = [ Msg0  'Handle before NEW-Action must be a Figure.' ];
   return
 end
 fig = h;
   h = [];

%-----------------------------------------------------
else

 ok = strcmp(get(h,'type'),'uicontrol');
 if ok 
   if any(strcmp(action,{'MOVE' 'UP'}))
     ok = ( strcmp(get(h,'style'),'pushbutton')  &  ...
            strcmp(get(h,'tag'),'MHELP_CENTER')           );
   else
     ok = ( strcmp(get(h,'style'),'frame')      &  ...
            strcmp(get(h,'tag'),'MHELP_FRAME')            );
   end
 end
 if ~ok
   Msg = [ Msg0  'Handle before ' action '-Action must be a ' ...
           'valid HelpButton.' ];
   return
 end
 fig = get(h,'parent');

end
%-----------------------------------------------------


  VarArg = varargin;
  VarArg = VarArg(:);
 


switch action

%*********************************************************
case 'NEW'


try

 [ MsgC,DefaultCData ] = readxpm('mouse.xpm');

catch

   MsgC = lasterr;

end

%----------------------------------------------------

if ~isempty(MsgC)

  DefaultCData = [ ...

'..............'
'..............'
'....?????.....'
'...?...????...'
'..??....????..'
'.????...?????.'
'.????...?????.'
'.????...?????.'
'..??....????..'
'........????..'
'.......????...'
'.......???....'
'......???.....'
'......??......'
'......?.......'
'......?.......'
'......?.......'
'..............'
'..............'
'.....???......'
'....?????.....'
'....?????.....'
'....?????.....'
'.....???......'
'..............'
'..............' ];

  DefaultCData                                        = abs(DefaultCData);
  DefaultCData( find( DefaultCData == double('.') ) ) = NaN;
  DefaultCData( find( DefaultCData == double('?') ) ) =  0 ;

  DefaultCData = DefaultCData(:,:,[1 1 1]);

end
%----------------------------------------------------


  % Title of ToolTip
  Title  = struct( 'Left'    , { 'Left Button'   } , ...
                   'Middle'  , { 'Middle Button' } , ...
                   'Right'   , { 'Right Button'  }       );

  % String of ToolTip
  String = struct( 'Left'    , { '' } , ...
                   'Middle'  , { '' } , ...
                   'Right'   , { '' }       );


  % CallBack for Buttons
  CallBack = struct( 'Center'    , { '' } , ...
                     'Left'    , { '' } , ...
                     'Middle'  , { '' } , ...
                     'Right'   , { '' }       );


  % Handle of Buttons
  Handle  = struct( 'Frame'    , { [] } , ...
                    'Center'     , { [] } , ...
                    'Left'     , { [] } , ...
                    'Middle'   , { [] } , ...
                    'Right'    , { [] }        );

  % Not allowed Properties in VarArg, for New and Set
  NotAllwd = { 'style'
               'tooltipstring' 
               'userdata' 
               'tag' 
               'buttondownfcn'
               'deletefcn'         };

  ud = struct( 'ButtonNumber'   , {  2  }          , ...
               'Title'          , { Title  }       , ...
               'String'         , { String }       , ...
               'CData'          , { DefaultCData } , ...
               'CallBack'       , { CallBack     } , ...
               'Position'       , { [ 10 10 ] }    , ...  % [ Left  Bottom ]
               'PixelOffset'    , { [  0  0 ] }    , ...  % [ Width Height ] arround CData
               'BorderWidth'    , {  3  }          , ...  % PixelWidth of PushButtonBorder
               'FrameWidth'     , {  1  }          , ...  % PixelWidth of FrameBorder
               'ButtonMinWidth' , {  4  }          , ...  % MinWidth of Button (3 Buttons)
               'DefaultCData'   , { DefaultCData } , ...
               'Handle'         , { Handle }       , ...  % HandleStructure
               'HH'             , { []     }       , ...  % Vector of Handles
               'SetFields'      , { [ 1 2 3 4 5 6 7 ] } , ...
               'NotAllwd'       , { NotAllwd }                 );
           

  ud.ButtonMinWidth = ud.ButtonMinWidth + mod(ud.ButtonMinWidth,2); 

  fields = fieldnames(ud.Handle);

  nh = size(fields,1);

  ud.HH    = zeros(nh,1);
  ud.HH(1) = uicontrol( 'style'    , 'frame'       , ...
                        'parent'   , fig           , ...
                        'tag'      , 'MHELP_FRAME' , ...
                      'busyaction' , 'cancel'      , ...
                   'interruptible' , 'off'                   );


  ud.Handle.Frame = ud.HH(1);

  set(ud.Handle.Frame,'userdata',ud);

  [Msg,IsSet,VarIn] = mhelp(ud.Handle.Frame,'Check',VarArg{:});

  if ~isempty(Msg) 
    delete(ud.Handle.Frame)
    return
  end

  ud0 = ud;
  ud  = get(ud0.Handle.Frame,'userdata');

  fig = get(ud.Handle.Frame,'parent');


  for ii = 2 : nh
      
    ud.HH(ii) = uicontrol( ...
                   'style'   ,'pushbutton' , ...
                   'parent'  , fig         , ...
                   'tag'     , [ 'MHELP_' upper(fields{ii}) ] , ...
                'busyaction' , 'cancel'    , ...
             'interruptible' , 'off'             );

    ud.Handle = setfield(ud.Handle,fields{ii},ud.HH(ii));

  end

  set( ud.Handle.Center,'tooltipstring' , 'Click with Right Button and Drag to Move.');
  
  % Update UserData

  set(ud.Handle.Frame,'userdata',ud);


  % Set all Values of IsSet to 1
   fields = fieldnames(IsSet);

   is_set = cell(2,size(fields,1));
   is_set(1,:) = fields(:)';
   is_set(2,:) = { 1 };
 
   IsSet = struct(is_set{:});

  Msg = mhelp(ud.Handle.Frame,'SetProp',IsSet,VarIn{:});


  if ~isempty(Msg)
     delete(ud.HH)
     return
  end


  % StringFormat for Handle
   clear eps
   Hform = sprintf( '%%.%.0fg' , ceil(abs(log(eps)/log(10)))+1 );

   HFrame = sprintf(Hform,ud.Handle.Frame);

   DelFcn  = ['mhelp(' HFrame ',''Delete'');'];
   DownFcn = ['mhelp(' HFrame ',''Down'');'];

   set( ud.Handle.Frame , 'DeleteFcn'     , DelFcn  );
   set( ud.Handle.Center, 'ButtonDownFcn' , DownFcn );


  out1 = ud.Handle.Frame;
  out2 = ud;

%**********************************************************************
case 'SET'

  if isempty(VarArg)
    return
  end

  [Msg,IsSet,VarIn] = mhelp(h,'Check',VarArg{:});

  if ~isempty(Msg) 
    return
  end

  Msg = mhelp(h,'SetProp',IsSet,VarIn{:});

 
  if isempty(Msg)
    out1 = get(h,'UserData');
  end

%**********************************************************************
case 'CHECK'

  if ( mod(size(VarArg,1),2) ~= 0 )
    Msg = [ 'Additional Inputs must contain Property-Value-Pairs.' ];
  else 
    VarArg = reshape(VarArg,2,size(VarArg,1)/2)';
    if ~iscellstr(VarArg(:,1))
     Msg = [ 'Additional Inputs must contain Property-Value-Pairs.' nl ...
             'Properties must be Strings.' ];
    end
  end

  if ~isempty(Msg)
     Msg = [ Msg0  Msg ];
     return
  end

  ud = get(h,'userdata');

  fields = fieldnames(ud);

  %-----------------------------------------------------
  % Search for Inputs, refering UserData of ud.SetFields
  %
  %    ButtonNumber: 2
  %           Title: [1x1 struct]
  %          String: [1x1 struct]
  %           CData: [24x14x3 double]
  %        CallBack: [1x1 struct]
  %        Position: [10 10]
  %     PixelOffset: [0 0]

  fields = fields(ud.SetFields);

  nf = size(fields,1);

  is_set      = cell(2,nf);
  is_set(1,:) = fields(:)';
  is_set(2,:) = { 0 };       % Will Set to 1 if Property found in VarArg

 
  nv = size(VarArg,1);

  is_ud  = zeros(nv,1);
 
  for ii = 1 : nv
     jj = find(strcmp(lower(VarArg{ii,1}),lower(fields)));
     if ~isempty(jj)

       val = VarArg{ii,2};
       msg = '';

       %------------------------------------------------------------------------
       if strcmp(fields{jj},'ButtonNumber')

         if ~( isequal(val,2)  |  isequal(val,3) )
           msg = [ 'Value for '  fields{jj} ' must be the Number 2 or 3.' ];
         end

       %------------------------------------------------------------------------
       elseif   any(strcmp(fields{jj},{'Title' 'String' 'CallBack'}))

         if ~isstruct(val)

           msg = [ 'Value for '  fields{jj} ' must be the Structure.' ];

         else
           if prod(size(val)) ~= 1
             msg = [ 'Size of Structure for ' fields{jj} ' must be [ 1 by 1 ].' ]; 
           end
           fn0 = lower(fieldnames(getfield(ud,fields{jj})));
           fn  = fieldnames(val);
           nf  = size(fn,1);
           ok  = zeros(nf,1);
           for kk = 1 : nf
             ok(kk) = any(strcmp(lower(fn{kk}),fn0));
             if ok(kk)
               ok(kk) = find(strcmp(lower(fn{kk}),fn0));
             end
           end
           if any(~ok)
              kk = find(~ok);
              msg = [  msg  nl(1:(end*(~isempty(msg)))) ...
                      'Invalid FieldNames for Structure ' fields{jj} ': ' nl ...
                        strhcat(fn(kk),', ',length(kk)+1)  ];
           else
             for kk = 1 : nf
               val1 = getfield(val,fn{kk});
               ok(kk) = ok(kk) * ( ischar(val1) &  ...
                                   ( prod(size(val1)) == size(val1,2) )  );
             end
             if any(~ok)
                kk = find(~ok);
                msg = [  msg  nl(1:(end*(~isempty(msg)))) ...
                        'Invalid Values for Fields of Structure ' fields{jj} ': ' nl ...
                          strhcat(fn(kk),', ',length(kk)+1)  nl ...
                        'Values must be a String.'   ];
             else
               val0 = getfield(ud,fields{jj});
               fn0  = fieldnames(val0);
               for kk = 1 : nf
                 val0 = setfield(val0,fn0{ok(kk)},getfield(val,fn{kk}));
               end
               val = val0;
             end             
           end 
         end           

       %------------------------------------------------------------------------
       elseif   any(strcmp(fields{jj},{'Position'  'PixelOffset'}))

         if ~isnumeric(val)  |  ~isequal(size(val),[1 2])
            msg = [ 'Value for '  fields{jj} ' must be a 2-dimensional, ' ...
                    'Row Vector.' ];
         else
           if any( ( mod(val(:),1) ~= 0 )  |  ~isfinite(val(:))  )
             msg = [ 'Values for ' fields{jj} ' must be Integers.' ]; 
           end
         end 
          
       %------------------------------------------------------------------------
       elseif strcmp(fields{jj},'CData')
         if ischar(val)
           ok = ( prod(size(val)) == size(val,2) );
           if ok
              ok = strcmp(lower(val),'default');
              if ok
                val = ud.DefaultCData;
              else
              % Check for XPM-Image
                [msg,val] = readxpm(val);
                ok = isempty(msg);
              end
           end                
           if ~ok
              msg = [ 'String Value for '  fields{jj} ' must be "default" ' ...
                       'or a valid FileName of an XPM-Image.' ];
           end
         elseif isnumeric(val)  
            ok = ( ( ndims(val) == 3 )  &  ~isempty(val) );
            if ok
               ok = all( ( ( val(:) >= 0 )  & ...
                           ( val(:) <= 1 )  & ...
                            isfinite(val(:))      ) | isnan(val(:))   ); 
            end
            if ~ok
              msg = [ 'Numeric Value for '  fields{jj} ' must be a 3-dimensional, ' ...
                      'nonemtpy TrueColorMatrice with Values between 0 and 1 or NaN.' ];
            end
         else
             msg = [ 'Values for ' fields{jj} ' must be a String "default" or ' ...
                     ' a valid XPM-FileName or a 3-dimensional TrueColorMatrice.' ]; 
         end 
          
       end
       %------------------------------------------------------------------------

       if isempty(msg)
         if strcmp(fields{jj},'CData')
           v1 = val;
           v1(find(isnan(v1))) = 0;
           v2 = getfield(ud,fields{jj});
           v2(find(isnan(v2))) = 0;
           is_set(2,jj) = { ~isequal(v1,v2) };
         else
           is_set(2,jj) = { ~isequal(val,getfield(ud,fields{jj})) };
         end
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
  
  VarArg(find(is_ud),:) = [];


  % Check if BackGroundColor changes

  bc0 = get(h,'backgroundcolor');  

  %------------------------------------------------------------------------
  % Check UIControlProperties
  if ~isempty(VarArg) 

    nn = size(ud.NotAllwd(:),1);
    ok = zeros(nn,1);
    for ii = 1 : nn
      ok(ii) = ~any(strcmp(ud.NotAllwd{ii},lower(VarArg(:,1))));
    end

    if any(~ok)
      kk = find(~ok);
      Msg = [  Msg0 'Not allowed Properties: ' nl ...
                strhcat(ud.NotAllwd(kk),', ',length(kk)+1) ];
      return
    end             

    VarArg = VarArg';

    try 
      set(h,VarArg{:});
    catch
      Msg = lasterr;
    end

    if ~isempty(Msg)
      Msg = [ Msg0  Msg ];
      return
    end

  end

  bc  = get(h,'backgroundcolor');


  % Store updated UserData

  set(h,'userdata',ud);

  out1 = struct(is_set{:});
  out2 = VarArg;

  out1.CData = ( out1.CData | ~isequal(bc0,bc) );

  
%**********************************************************************
case 'SETPROP'

  nv       = size(VarArg,1);

  IsSet    = VarArg{1};

  VarArg = VarArg(2:nv);

  % IsSet =           0  |  1 
  % 
  %    ButtonNumber: 
  %           Title:
  %          String:
  %           CData:
  %        CallBack:
  %        Position:
  %     PixelOffset:

  ud = get(h,'userdata');

 is3 = ( ud.ButtonNumber == 3 );

 %---------------------------------------------------------------------
 % New Positioning 
 if IsSet.ButtonNumber | IsSet.Position | IsSet.PixelOffset  | ...
     ~isequal(size(ud.CData),size(get(ud.Handle.Center,'cdata')))

    csi = [ size(ud.CData,2)  size(ud.CData,1) ] ;  % [ Width Height ] of Image    
    csi = csi + ud.PixelOffset;

    % Minimum Width of Center, without Border
    Width = 3*(ud.ButtonMinWidth+2*ud.BorderWidth) - 2*ud.BorderWidth ;

    Width  = max([Width csi(1)]) + 2*ud.BorderWidth;

    Width  = ud.ButtonNumber * ceil( Width / ud.ButtonNumber ); 
    
    Height = max([Width-2*ud.BorderWidth csi(2)]) + 2*ud.BorderWidth;

    % Buttons
    BWidth  = Width / ud.ButtonNumber;
    BHeight = ceil( Width / 2 );

    % Frame
    FWidth  = Width + 2 * ud.FrameWidth;
    FHeight = Height - ud.BorderWidth + BHeight + 2*ud.FrameWidth;

    pos      = zeros(5,4);
    pos(1,:) = [ ud.Position                FWidth  FHeight ];  % Frame
    pos(2,:) = [ ud.Position+ud.FrameWidth   Width   Height ];  % Center
    pos(3,:) = [ pos(2,1)  pos(2,2)+Height-ud.BorderWidth  ...
                                     BWidth  BHeight ];  % Left 
    pos(4,:) = [ pos(3,1)+BWidth*is3  pos(3,[2 3 4]) ];  % Middle
    pos(5,:) = [ pos(4,1)+BWidth      pos(3,[2 3 4]) ];  % Right    


   % Check if Center inside Figure
    fig = get(h,'parent');
    uni = get(fig,'units');
    set(fig,'units','pixels');
    figpos = get(fig,'position');
    set(fig,'units',uni);
 

     [dx,dy] = chk_pos( pos(2,:)-[ 0  0  0 ud.BorderWidth ] , ...
                        figpos , 3*ud.BorderWidth );


    pos(:,1) = pos(:,1) + dx;
    pos(:,2) = pos(:,2) + dy;

    % Update Position
      ud.Position = pos(1,[1 2]);

    % Update UserData
      set(h,'userdata',ud);

    for ii = 1 : 5
      set(ud.HH(ii),'units'    ,'pixels', ...
                    'position' ,pos(ii,:)      );
    end

 
 end
 % New Positioning & CData
 %---------------------------------------------------------------------

 %---------------------------------------------------------------------
 % New CData
 if IsSet.CData    

    % Set undefined ImageColors to BackGroundColor
    %  BackGroundColor of Frame allready updated during  mhelp( Check )

      bc = get(ud.Handle.Frame,'backgroundcolor');

      ii = find(isnan(sum(ud.CData,3)));
 
    si12 = size(ud.CData,1)*size(ud.CData,2);

      ud.CData(ii+0*si12) = bc(1);
      ud.CData(ii+1*si12) = bc(2);
      ud.CData(ii+2*si12) = bc(3);

    set(ud.Handle.Center,'cdata',ud.CData);

 end
 %---------------------------------------------------------------------

 %---------------------------------------------------------------------
 % New ToolTip
 if IsSet.Title  |  IsSet.String  
    fields = fieldnames(ud.Title);
    nf     = size(fields,1);
    for ii = 1 : nf
       tt = getfield(ud.Title ,fields{ii});
       ss = getfield(ud.String,fields{ii});
       hh = getfield(ud.Handle,fields{ii});
       
        nd = 0;   % DelimiterNumber
      if ~isempty(tt) & ~isempty(ss)
         i10 = [ 0  find(abs(tt)==10) size(tt,2)+1 ];
         nd = max([ nd  diff(i10) ]);
      end

      nld = char( 10 * ones(1,( nd ~= 0 ))); % NewLineDelimiter

      set(hh,'tooltipstring' , [ tt nld char(double('-')*ones(1,nd)) nld ss ]); 

    end
    % ii
 end
 % New TooTip
 %---------------------------------------------------------------------

 
 %---------------------------------------------------------------------
 % New CallBack
 if IsSet.CallBack
    fields = fieldnames(ud.CallBack);
    nf     = size(fields,1);
    for ii = 1 : nf
       set( getfield(ud.Handle,fields{ii}) , 'callback' , ...
            getfield(ud.CallBack,fields{ii}) );
    end
 end
 %---------------------------------------------------------------------

 %---------------------------------------------------------------------
 % Set other Properties
 if ~isempty(VarArg)
   set( ud.HH , VarArg{:} );
 end

 %---------------------------------------------------------------------
 % Set Visibility of MiddleButton

   Sets = { 'off'  'on' };

   ii = ( is3  &  strcmp(get(h,'visible'),'on') );

   set(ud.Handle.Middle,'visible',Sets{1+ii});


    
%**********************************************************************
case 'DOWN'

   ud = get(h,'userdata');
  fig = get(h,'parent');

  uni = get([ud.HH;fig],'units');


  set([ud.HH;fig],'units','pixel');

  % StringFormat for Handle
   clear eps
   Hform = sprintf( '%%.%.0fg' , ceil(abs(log(eps)/log(10)))+1 );

   HCenter = sprintf(Hform,ud.Handle.Center);

   MoveFcn = ['mhelp(' HCenter ',''Move'');'];
     UpFcn = ['mhelp(' HCenter ',''Up'');'];

   pos = get(ud.HH,'position');
   pos = cat(1,pos{:});

  
   ud1 = struct( 'Figure'    , { fig   }   , ...
                 'HH'        , { ud.HH }   , ...
                 'Units'     , { uni   }   , ...
                 'FigPos'    , { get(fig,'position')     } , ...
                 'StartPos'  , { pos                     } , ...
                 'StartPoint', { get(fig,'CurrentPoint') } , ...
                 'BorderWidth', { ud.BorderWidth         } , ...
                 'MotionFcn' , { get(fig,'windowbuttonmotionfcn') } , ...
                 'UpFcn'     , { get(fig,'windowbuttonupfcn')     }       );

   set(ud.Handle.Center,'userdata',ud1);

   set( fig , 'windowbuttonmotionfcn' , MoveFcn , ...
              'windowbuttonupfcn'     ,   UpFcn        );


%**********************************************************************
case 'MOVE'

   ud = get(h,'userdata');

   cp = get(ud.Figure,'currentpoint');

   pos      = ud.StartPos;

   pos(:,1) = ud.StartPos(:,1) + cp(1)-ud.StartPoint(1) ;
   pos(:,2) = ud.StartPos(:,2) + cp(2)-ud.StartPoint(2) ;

   [dx,dy] = chk_pos( pos(2,:)-[ 0  0  0 ud.BorderWidth ] , ...
                      ud.FigPos , 3*ud.BorderWidth );


    pos(:,1) = pos(:,1) + dx;
    pos(:,2) = pos(:,2) + dy;

   for ii = 1 : size(pos,1);
     set(ud.HH(ii),'position',pos(ii,:));
   end

   eval(ud.MotionFcn,'0;');

%**********************************************************************
case 'UP'

   ud = get(h,'userdata');

   set( fig , 'windowbuttonmotionfcn' , ud.MotionFcn , ...
              'windowbuttonupfcn'     , ud.UpFcn        );
 
   set(h,'userdata',[])

   % Update Position in FrameUserdata

     ud1 = get(ud.HH(1),'userdata');
     pos = get(ud.HH(1),'position');

     ud1.Position  = pos([1 2]);


   % Update UserData

     set(ud.HH(1),'userdata',ud1);


   % Set Units to default when MouseDown

   for ii = 1 : size(ud.HH,1)
     set(ud.HH(ii),'units',ud.Units{ii});
   end

   set(ud.Figure,'units',ud.Units{end});


%**********************************************************************
case 'RESIZE'

    ud = get(h,'userdata');
 
   % Update FigurePosition

    fig    = get(h,'parent');
    figuni = get(fig,'units');

     set(fig,'units','pixels');
   
     FigPos = get(fig,'Position');

     set(fig,'units',figuni);

   % Get HandlePositions

     uni = get([ud.HH;fig],'units');

     set([ud.HH;fig],'units','pixel');

     pos = get(ud.HH,'position');
     pos = cat(1,pos{:});

   % Correct Positions

     [dx,dy] = chk_pos( pos(2,:)  , ...
                        FigPos , 3*ud.BorderWidth );


     pos(:,1) = pos(:,1) + dx;
     pos(:,2) = pos(:,2) + dy;


     for ii = 1 : size(pos,1);
       set(ud.HH(ii),'position',pos(ii,:));
       set(ud.HH(ii),'units',uni{ii});
     end

   set(h,'userdata',ud);


%**********************************************************************
case 'DELETE'
 
 ud = get(h,'userdata');

 for ii = 2 : size(ud.HH(:),1) 
   try
    delete(ud.HH(ii));
   end
 end


%**********************************************************************
otherwise


 Msg = [ Msg0  'Invalid action.' ];


end

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function  str = strhcat(str,del,n)

% STRHCAT  Concatenates Strings into ONE
%
% STRHCAT( StringArray , Delimiter )

nl = char(10);

str = str(:);

str( : , 2 ) = { del };

str(n:n:size(str,1),2) = {nl};

str(    size(str,1),2) = {''};


str = str';

str = cat(2,str{:});


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function [ dx , dy ] = chk_pos(tpos,fpos,lim);

% CHK_POS   Checks, if a Target Position is inside a FramePosition
%
% [ DX , DY ] = CHK_POS( TargetPosition , FramePosition , Limit )
%
%  Positions are: [ Left Bottom Width Height ]
%
% TargetPosition is relative to FramePosition
%
%  Limit is the PixelLimit for the Target to be inside the Frame
%


  % FrameCorners             

    ini = [ 1  1 ; ...  % Upper right  First !!!
            0  1 ; ...  % Upper left
            1  0 ; ...  % Lower right
            0  0   ...  % Lower left  
          ];

      ni = size(ini,1);  % Number of Points
      ii = ones(1,ni);

      op = ( 1 - 2*ini );  % Operator  [ 0  1 ] --> [ 1  -1 ]

   % Compare this Points with opposite Corner of Target (1-ini);
   %  TargetCorner

    tpos = tpos(ii,[1 2]) + (1-ini) .* (tpos(ii,[3 4])-1);

   %  FigureCorner, lim Pixel less

    fpos = ones(ni,2) + ini  .* (fpos(ii,[3 4])-1) + op * lim;


   % Compare Positions

      ok = double( op.*tpos >= op.*fpos );

   % Index of First Point Outside

      ok = sum(cumprod(ok,1),1)+1;

   % ZERO, if First OutsidePoint not in Data

       f = ( ok <= ni );          
      ok = ni + f .* ( ok - ni );
      
      dx = f(1) * ( fpos(ok(1),1) - tpos(ok(1),1) );
      dy = f(2) * ( fpos(ok(2),2) - tpos(ok(2),2) );
