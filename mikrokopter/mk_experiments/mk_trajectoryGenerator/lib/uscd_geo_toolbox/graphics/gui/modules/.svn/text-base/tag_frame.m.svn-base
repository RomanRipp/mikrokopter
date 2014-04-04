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
