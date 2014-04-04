function [msg,out]=tab_list(H,action,varargin)

% TAB_LIST  Creates a Tabular
%
% [Msg , Out ] = TAB_LIST( Handle , Property1 , Property1Value , ...
%                                   Property2 , Property2Value , ... )
%
% Valid Handles are:  
%
%  FigHandle    to create a new Tabular in Figure of FigHandle
%               Out = TabHandle
%
%  TabHandle    to set properties of the Tabular of TabHandle
%
%  TabHandle,'delete' to delete the Tabular of TabHandle
%
% Valid Properties and there Values are:
%
%  'position' , [ Left Bottom ] , ...    PixelPosition of lower, left Corner 
%               [-Right -Top  ]             in Figure
%  'visible'  , 'on' | 'off'    , ...    Set's Visibility of Tabular
%  'enable'   , 'on' | 'off'    , ...    Set's Enability of Tabular
%
%  'vspace'     , V                      VerticalPixelSpace between UIControls
%                                          (default: 1 )
%  'sliderwidth', N , ...                PixelWidth of the ListboxSlider 
%                                          (default: 15)
%
% Note: TAB_LIST didn't use the ListBoxSlider, they are hidden  
%       by the following ListBox or an SliderUIControl.
%       If it's not perfect hidden, set the SLIDERWIDTH !
%
%  'ntab' , N , ...                      Number of TabularColumns 
%  'nrow' , R , ...                      Number of CharakerRows
%  'ncol' , [ C1 .. CN ] , ...           Number of CharakterColumns 
%                                          per TabularColumn   [ N by 1 ] 
%  'string' , { String1 .. StringN }     CellArray of N Strings (each String  
%             { Vrow Vcol String# ... }   with Vector for Rows and TabularColumns
%             { Vrow String# ... }
%  'form'   , { Format1 .. FormatN }     Format of the Strings, using by VAL2STR
%             { Vcol  Format#  }          with Vector for TabularColumns
%  'title'  , { String1 ,, StringN }     Title for TabularColumns
%             { Vcol String# }            with Vector for TabularColumns
%
%  'callback'  CallBack      Set User defined CallBack, after successfull
%                             'SLIDER', 'LIST', 'SET' , 'INS', or 'DEL' - Action
%               using:  feval( CBFcn{:} , Action , 1 , varargin )
%                  Actions are:  'TabSlider'  'TabList'  'TabSet'  'TabIns' 
%                                'TabDel'     'TabReset' 
%
%  'setquestfcn', SetQuestionFunction    Function to evaluate before 
%                                         'SET' and 'INS' - Action
%                                        by ControlButton, using by FEVAL, 
%                                        following action is: 'TabSetQuest' , 1 , 
%                                        function must answer with 0 or 1 or STR
%
%  'delquestfcn', DeleteQuestionFunction  Function to evaluate before 
%                                          'DEL' and 'RESET' - Action
%                                         by ControlButton, using by FEVAL, 
%                                        following action is: 'TabDelQuest' , 1
%                                        function must answer with 0 or 1
%
%  'msgfcn', MessageFunction  Function to evaluate for Messages, using by FEVAL, 
%                                       
%  'button' ,   Adds a Button to ControlBar, following Inputs are
%                   UIControl-Property-Value-Pairs
%
%  'visible' , 'on' | 'off'   sets Visibility of TabList
%  'enable'  , 'on' | 'off'   sets  Enability of TabList
%
%  'userdata' { UserData }               CellArray, Length == R
%             { Vrow {UserData} }         with Vector for Rows
%
%  'value'   , V                         Set ListValues to V
%  'setedit' , RowNumber                 Set ListRow into EditFields
%  'delete'  , RowNumbers                Deletes Rows
%  'delete'  , 0                         Deletes TabList
%
%
%  'val2str' , TabCol , Values           Gives formatedStrings and corresponding Values,
%  'val2str' , TabCol , {Strings}           in the Format of the TabularColumns
%                                        Out = {  [ 0 | 1 ]   Value   String  Msg1  Msg2 }
%                                          ( see VAL2STR )
%
%  'vedit'   , { 'on' | 'off' }          Set the Visibility of the EditButtons
%              { Vcol   'on' | 'off' }    with Vector of TabularColumns
%                                        [ 1 by 1 ]  CellStringArray
%                                        [ 1 by N ]  CellStringArray
%                                        [ {VC}  {1 by length(VC)} ] CellArray
%
%  To set the Visibility of the ControlButtons:
%  'cedit'   ,  'on' | 'off'
%  'cset'    ,  'on' | 'off'
%  'cins'    ,  'on' | 'off'
%  'cundo'   ,  'on' | 'off'     NOTE:  UNDO still under construction
%  'credo'   ,  'on' | 'off'     NOTE:  REDO still under construction
%  'cdel'    ,  'on' | 'off'
%  'creset'  ,  'on' | 'off'
% 
%
% Any Other UIControlProperties like 'BackGroundColor' or 'Font...'  are allowed.
%  For instance  >> tab_list(TabHandle,'min',0,'max',2)
%   allows multiple ListBox selection.
%
% An "ZERO". Column with the ListCounter is set automaticly. The User can't
%  change this Column !!!
%
% 
% The TabHandle is the surrounding FrameUIControl, it's UserData contains
%  all Informations about the other Handles and PropertyValues in a StructArray.
%  So  >> ud = get(TabHandle,'userdata');  % gives this StructArray
%      >> ud.string  % are the current Strings of the Tabular
%  Discover this StructArray for you own.
%
% Following Properties you must get directly from the Tabular(Frame)Handle:
%
%  get(TabHandle,'max')     Number of ListEntries   
%  get(TabHandle,'value')   Value  of ListBoxes (Marked Entries)
% 
% How to work ineractively with the Tabular?
%  If you write in the EDIT-Fields below and press the SET-Button,
%   the strings from the EDIT-Fields are Set into the List. The Position
%   in the list is defined by the Number in the first ("ZERO") Field.
%   If this Number is empty, the strings will be append to the old List.
%  By pressing the EDIT-Button, the Strings of the first marked ListEntries
%   will set to the EDIT-Fields.
%  By pressing the DEL-Button, the first marked ListEntry will be deleted
%   after request.
%
% required function: RMBLANK
%
%
%  Valid Formats are:
%                                                   '#'   String
%  '%##.##'   for using by SPRINTF
%
%  'geo#'     for Geographic Coordinates  Longitude: 1     ##°##.##' <E/W>
%                                         Latitude : 2     ##°##.##' <N/S>
%
%  'time#'    for TimeConversation, first Value Day: 1    DD HH:MM
%                                              Hour: 2       HH:MM
%                                              Hour: 3       HH:MM:SS
%                                            Minute: 4          MM:SS
%
%
% requires functions VAL2STR, STR2DAY, DAY2DHMS, RMBLANK
%

% Internal Actions:
%
%  'set_prop'  Set ListboxProperties
%  'list'     ListBoxCallback, set ListBoxValue and ListBoxTop 
%  'slider'   SliderCallback, vertical Slider for ListBox's
%  'edit'     EditBox
%  'editbutton'  EditControl
%  'setbutton'       SetControl
%  'insbutton'       InsControl
%  'delbutton'      DeleteControl 
%  'resbutton'      DeleteControl 
%



Nin  = nargin;
Nout = nargout;


msg = '';
out = [];


msg0 = ' TAB_LIST: ';   % Add "action"   directly before  switch action


nl = char(10);



if Nin == 0
 H = gcf;
 action = 'new';
end


%---------------------------------------------------------------------
% Check Handle
ok = ( isnumeric(H)  &  ( prod(size(H)) == 1 ) );

if ok
   ok = ishandle(H);
   if ok
  
    typ = get(H,'type');
    tag = get(H,'tag');
  
    ok = ( strcmp(typ,'figure')  |  ...
           ( strcmp(typ,'uicontrol') & strcmp(tag,'TAB_FRAME') )  );
   end
end

if ~ok
  msg = [ msg0 'Handle must be a Figure or valid TabularHandle.' ];
  return
end

if Nin < 2
   if strcmp(typ,'figure')
      action = 'new';
   end
end



%---------------------------------------------------------------------

if Nin < 3
  VarIn = cell(0,1);
else
  VarIn = varargin;
  VarIn = VarIn(:);
end


%---------------------------------------------------------------------
% Check Action
%           action   min(Nin)   set_prop

actions = { 'list'        1          0
            'slider'      0          0
            'edit'        1          0
            'editbutton'  0          0
            'setbutton'   0          0
            'insbutton'   0          0
            'undobutton'  0          0
            'redobutton'  0          0
            'delbutton'   0          0
            'resetbutton' 0          0
            'undo'        0          0
            'redo'        0          0
            'delete'      1          0
            'reset'       1          0
            'setedit'     1          0
            'value'       1          0
            'val2str'     2          0
            'set_prop'    2          0
            'visible'     1          0
            'enable'      1          0
            'string'      1          1
            'userdata'    1          1
            'setquestfcn' 1          1
            'delquestfcn' 1          1
            'msgfcn'      1          1
            'button'      0          0
            'vedit'       1          1
            'cedit'       1          1
            'cset'        1          1
            'cins'        1          1
            'cundo'       1          1
            'credo'       1          1
            'cdel'        1          1
            'creset'      1          1      };
            


ok = ( ischar(action)  &  ~isempty(action)  &  ...
        ( prod(size(action)) == size(action,2) )      );

if ok

   action = lower(action);

   if ( strcmp(typ,'figure')  &  ~strcmp(action,'new') )
      VarIn = cat(1,{action},VarIn);
      action = 'new';
   end

else

   msg = [ msg0 'Second Input must be an Action or PropertyName (String).' ];
   return

end



%--------------------------------------------------
if ~strcmp(action,'new');

  ud = get(H,'userdata');

  is_act = [];

  % Search External Action
  is_act = find(strcmp(action,actions(:,1)));

  if isempty(is_act)
     VarIn  = cat(1,{action},VarIn);
     action = 'set_prop';
  else
   if actions{is_act,3}
   % a SetProp-Action !!!!!!!!
     VarIn  = cat(1,{action},VarIn);
     action = 'set_prop';
   end
  end


  jj = find(strcmp(action,actions(:,1))); 
  % SET_PROP  Value
  if size(VarIn,1) < actions{jj,2}
    msg = [ msg0 'Not enough InputArguments.' ];
    return
  end

end


Nin = size(VarIn,1);

%**************************************************************************
if any( strcmp( action , { 'new'  'set_prop' } ) )

% Search for SET - Properties

  if ~isempty(VarIn)

    ok = ( mod(Nin,2) == 0 );
    if ok
      Nin   = Nin/2;
      VarIn = reshape( VarIn , 2 , Nin );
      ok    = iscellstr(VarIn(1,:));
    end

    if ~ok
      msg = [ msg0 ' Inputs must be Property-Value-Pairs, Properties must be Strings.' ];
      return
    end

  end

  %!!!!!!!!!!!!!!!!!
  % 5 CONTROL-Button on bottom: SET  EDIT  UNDO  REDO  DEL
  % 1 Column == PorpertyName
  % 2 Column == Index in VarIn
  % 3 Column == DefaultValues

  slw = 15 + ( 3 + 3 * ( get(0,'screenpixelsperinch') > 96 ) ) * ...
             strcmp( upper(computer) , 'PCWIN' );

  set_prop = { 'ntab'          0    2        ; ...
               'nrow'          0    5        ; ...
               'ncol'          0   10        ; ...
               'position'      0   get(0,'defaultuicontrolposition')  ; ...
               'value'         0    1        ; ...
               'string'        0   ''        ; ...
               'callback'      0   {}        ; ...
               'border'        0   5         ; ...
               'vedit'         0   'on'      ; ...
               'cedit'         0   'on'      ; ...
               'cset'          0   'on'      ; ...
               'cins'          0   'on'      ; ...
               'cundo'         0   'off'     ; ...
               'credo'         0   'off'     ; ...
               'cdel'          0   'on'      ; ...
               'creset'        0   'on'      ; ...
               'form'          0   'char'    ; ...
               'title'         0   ''        ; ...
               'sliderwidth'   0   slw       ; ...
               'userdata'      0   cell(1,1) ; ...
               'setquestfcn'   0   {}        ; ...
               'delquestfcn'   0   {}        ; ...
               'msgfcn'        0   {}        ; ...
               'vspace'        0    1              };


  Nsprop = size(set_prop,1);


  if ~strcmp(action,'new')
  % Default Values from Existing Handle

    set_prop(:,3) = { ud.ntab ; ...
                      ud.nrow ; ...
                      ud.ncol ; ...
                      ud.position      ; ...
                      get(H,'value')   ; ...
                      ud.string        ; ...
                      ud.callback      ; ...
                      ud.border        ; ...
                      ud.edit{2}(2:length(ud.edit{2})) ; ...
                      ud.cntr{2}{1}    ; ...
                      ud.cntr{2}{2}    ; ...
                      ud.cntr{2}{3}    ; ...
                      ud.cntr{2}{4}    ; ...
                      ud.cntr{2}{5}    ; ...
                      ud.cntr{2}{6}    ; ...
                      ud.cntr{2}{7}    ; ...
                      ud.form          ; ...
                      ud.title         ; ...
                      ud.sliderwidth   ; ...
                      ud.userdata      ; ...
                      ud.setquestfcn   ; ...
                      ud.delquestfcn   ; ...
                      ud.msgfcn        ; ...
                      ud.vspace               };
  end



  %-----------------------------------
  % Select SET - Properties from VarIn
  %...................................


  if isempty(VarIn)

     is_prop = 0;

     VarIn = { '' ; [] };

  else

    is_prop = ones(Nin,1);      % VarIn is UIControl-Property

    VarIn(1,:) = lower(VarIn(1,:));

  end

  for ii = 1 : Nin

     kk = find( strcmp( VarIn{1,ii} , set_prop(:,1) ) );
  
     if ~isempty(kk)

        is_prop(ii) = 0;

        kk = kk(end);

        prop = set_prop{kk,1};

        val = VarIn{2,ii};

        is_char = ischar(val); 

        is_cell = iscell(val);


        prop_clb  = { 'callback'  'setquestfcn'  'delquestfcn'  'msgfcn' };
        prop_sets = { 'cedit'  'cset'  'cins'  'cundo'  'credo'  'cdel'  'creset' };

       chk_clb  = any(strcmp(prop,prop_clb));
       chk_sets = any(strcmp(prop,prop_sets));

       chk_usd  = strcmp(prop,'userdata');


        if is_char & chk_clb

           val = { val };

           is_cell = 1;

        end

        if is_cell

          switch prop

            %------------------------------------------------------------------
            case 'string'
  
              val = val(:)';
              is_char1 = 1;
              for jj = 1:length(val)
                 is_char1 = ( is_char1 &  ...
                     (  ( isnumeric(val{jj}) & any(jj==[1 2]) ) | ...
                            ischar(val{jj})  |  iscellstr(val{jj})      ) );
              end

              is_char = ( is_char |  is_char1 | iscellstr(val) );

            %------------------------------------------------------------------
            case 'vedit'

              sets = {'off' ; 'on'};
              val = val(:)';

              try
                is_char = ( isnumeric(val{1}) &  ...
                            iscellstr( val(2:length(val)))  &  ...
                            ( length(val{1}(:)) == length(val)-1  ) );
              catch
                is_char = 0;
              end
            
              is_char = ( is_char | iscellstr(val) );

            %------------------------------------------------------------------
            case { 'title'  'form' }

              try
                is_char = ( isnumeric(val{1}) & iscellstr(val(2:length(val))) );
                if isnumeric(val{1}) & ~is_char & strcmp(prop,'title')
                   is_char = 1;
                   for jj = 2 : prod(size(val))
                     is_char = ( is_char & ( ischar(val{jj}) | iscellstr(val{jj}) ) );
                   end
                end
              catch
                is_char = 0;
              end

              is_char = ( is_char | iscellstr(val) );

              if  ~is_char & strcmp(prop,'title')
                 is_char = 1;
                 for jj = 1 : prod(size(val))
                   is_char = ( is_char & ( ischar(val{jj}) | iscellstr(val{jj}) ) );
                 end
              end

            %------------------------------------------------------------------
            case 'userdata'

              if size(val,2) == 2
                try
                   is_cell = ( isnumeric(val{1}) & iscell(val{2}) &  ...
                             ( length(val{1}(:)) == length(val{2}(:)) ) );
                catch
                   is_cell = 0;
               end
             else
               val = val(:); 
             end
 
             is_cell = ( is_cell  &   ~isempty(val) );

            %------------------------------------------------------------------
            case prop_clb
            
              % 1. CellElement must be a String, using FEVAL

              is_cell = isempty(val);
              if ~is_cell
                 is_cell = ( ischar(val{1}) &  ~isempty(val{1}) & ...
                             ( prod(size(val{1})) == size(val{1},2) )  );
              end

          end
          % switch

        end
        % is_cell


        %------------------------------------------------------------------

        sets = {'off' ;  'on'};

        is_sets = 1;

        if chk_sets

             try
               is_sets = any(strcmp( val , sets ));
             catch
               is_sets = 0;
             end

             if  is_sets 
               val = sets{find(strcmp( val , sets ))};
             end

        end

        %------------------------------------------------------------------
        if any( strcmp( prop , { 'ntab' 'ncol' 'nrow' ...
                                 'position'  'value'  ...
                                 'border'  'sliderwidth'  'vspace' } ) )

          if isnumeric(val)

             is_pos = strcmp(prop,'position');
             is_bor = strcmp(prop,'border');
 
             val_ok =  ( all(isfinite(val)) & ...
                        ( ( all(val>=1) & ~(is_pos | is_bor) )    | ...
                          ( all(val>=0) &   is_bor           )    | ...
                          ( ( size(val,2) >= 2 )  &  is_pos  )           )  );
          else
             val_ok = 0;
          end

        %------------------------------------------------------------------
        else 

          val_ok = ( (    chk_sets                       & is_sets )  | ...
                     (  ( chk_usd  | chk_clb )           & is_cell )  | ...
                     ( ~( chk_sets | chk_usd | chk_clb ) & is_char )         );
 
        end

        if val_ok 

          if strcmp(prop,'vedit')
             val0 = set_prop{kk,3};
             if ischar(val)
               val = [ {val}  val0(2:length(val0)) ]; 
             elseif iscellstr(val)
               val = [  val   val0(length(val)+1:length(val0)) ];
             else
                vv       = round(val{1}(:));
                   jj    = find( vv < 1  |  length(val0) < vv );
                vv(jj)   = [];
               val(jj+1) = [];
  
               val0(vv)  = val(2:length(val)); 
               val       = val0;
             end
          end

          % Check Value with old Value 
          if ~isequal( set_prop{kk,3} , val )  |  strcmp(prop,'position') 

             set_prop{kk,2} = ii;
             set_prop{kk,3} = val;

          end

        else

          msg = [ msg nl(1:(end*(~isempty(msg)))) ...
                  'Invalid Value for Property '''  set_prop{kk,1} ''''  ];

        end

     end
     % ~isempty(kk)

  end
  % ii = 1 : Nin


  if ~isempty(msg)
    msg = [ msg0 msg ];
    return
  end


  %-------------------------------
  % Select  FONT - Properties
  %...............................

  is_font = strmatch( 'font' , VarIn(1,:) );

  is_prop(is_font) = 0;



  %-------------------------------
  % Other Properties
  %...............................


  is_prop = find(is_prop);


   % StringFormat for Handle

   clear eps

   Hform = sprintf( '%%.%.0fg' , ceil(abs(log(eps)/log(10)))+1 );


end
% 'set_prop' or 'new'






% CALLBACK in ud.callback

callback_ok = 0;
callback_in = [];

% 1 if SET or DEL(Yes) or ListBoxAction



msg0 = [ ' TAB_LIST( ' upper(action) ' ): ' ];



%******************************************************************
switch action

case 'new'

   is_win = strcmp(upper(computer),'PCWIN');


   % Colors
   cc = [ 1     1      1
          0.9   1      1
          1.0   0.95   0.95  ];

   % ControlButtons:
   cntr = { 'Edit' cc(2,:)  'Edit first selected ListEntry'
            'Set'  cc(2,:)  'Set edited Entry to List'
            'Ins'  cc(2,:)  'Insert edited Entry in List'
            'Undo' cc(2,:)  ''
            'Redo' cc(2,:)  ''
            'Del'  cc(3,:)  'Delete selected ListEntry(s)'
          'Reset'  cc(3,:)  'Delete ALL ListEntry(s)'   }; 


   nc = size(cntr,1);
   
   ud = struct( ...
     'ntab'         , { 0 } , ...
     'nrow'         , { 0 } , ...   
     'ncol'         , { 0 } , ...   
     'border'       , { 0 } , ...   
     'position'     , { 0  } , ...  
     'sliderwidth'  , { 15+3*is_win  } , ... 
     'borderwidth'  , {  3-0*is_win  } , ... 
     'lineoffset'   , {  2-[0 2]*is_win  } , ... 
     'vspace'       , {  1  } , ...  
     'Hframe'        , { NaN  } , ...    % FrameHandle
     'Haxe'          , { NaN  } , ...    % DummyAxeHandle
     'Htext'         , { NaN  } , ...    % DummyText
     'Htitle'        , { NaN  } , ...    % TitleHandle of each TabularColumn
     'Hlist'         , { NaN  } , ...    % ListBoxHandles
     'Hslider'       , { NaN  } , ...    % RightSliderHandle
     'edit'         , { {  { NaN }   { 'on' } } } , ...   % EditBoxes below ListBox
     'cntr'         , { { NaN*ones(nc,1) cellstr(char(ones(nc,1)*'on'))  } } , ...
     'form'         , { { 'char' } } , ...
     'title'        , { {''}  }      , ... 
     'string'       , { {''}  }      , ... 
     'callback'     , {  ''   }      , ... 
     'userdata'     , {  {}   }      , ... 
     'setquestfcn'  , {  ''   }      , ... 
     'delquestfcn'  , {  ''   }      , ... 
     'msgfcn'       , {  ''   }            );

     
   % default: Undo and Redo 'off'

    ud.cntr{2}([4 5]) = {'off'};



   
    ud.Hframe =  uicontrol('style','frame', ...
                          'units'   , 'pixels' , ... 
                          'string',{''} , ...
                          'min'   ,0 , ...
                          'max'   ,0 , ...
                          'value' ,1 , ...
                          'tag'  ,'TAB_FRAME' , ...
                          'parent' ,H, ...
                          'visible','off');

    HFrame = sprintf( Hform , ud.Hframe );

    DelCB = [ 'tab_list(' HFrame ',''Delete'',-1);' ];

    set( ud.Hframe , 'DeleteFcn' , DelCB);


    % DummyAxe for DummyText

    ud.Haxe =  axes('parent'  , H , ...
                      'units'   , 'pixels' , ... 
                      'position', get(ud.Hframe,'position') , ...
                      'xlim'    , [ 0  1 ] , ...
                      'ylim'    , [ 0  1 ] , ...
                      'xtick'   , [] , ...
                      'ytick'   , [] , ...
                      'visible' , 'off'  , ...
                     'nextplot' , 'add' , ...
                      'tag'     , 'TAB_AXES'  , ...
             'handlevisibility' , 'callback'          );

    ud.Htext = text('parent' , ud.Haxe , ...
          'units' , 'pixels' , ...
          'position',[ 1 1 0 ] , ...
          'string'  , '' , ...
          'visible' , 'off' , ...
          'tag'     , 'TAB_TEXT'         );



    ud.Htitle = uicontrol('style','text', ...
                           'units'   , 'pixels' , ... 
                           'tag'  , 'TAB_TITLE' , ...
                           'parent' ,H, ...
                           'visible','off');



    CB = [ 'tab_list(' HFrame ',''list'',1);' ];

    % NumberTabularColumn
    ud.Hlist = uicontrol('style','listbox', ...
                        'units'   , 'pixels' , ... 
                        'min',0,'max',2, ...
                        'string','', ...
                        'callback', CB , ...
                        'tag' , 'TAB_LIST' , ...
                        'parent' ,H, ...
                        'visible','off');

    % Edit below each TabularColumn, First is NumberCounter
    %  UserData:  { Form  OldString }
    CB = [ 'tab_list(' HFrame ',''edit'',1);' ];

    ud.edit{1} = uicontrol('style','edit', ...
                           'units'   , 'pixels' , ... 
                           'tag'  , 'TAB_EDIT' , ...
                           'userdata',{'char',''}, ...
                           'callback', CB , ...
                           'userdata',{ '%3.0f'  ''  [] } , ...
                           'parent' ,H, ...
                           'visible','off');

    % Slider on the Right Side
    CB = [ 'tab_list(' HFrame ',''slider'');' ];
    ud.Hslider = uicontrol('style','slider', ...
                          'units'   , 'pixels' , ... 
                          'min'   ,1 , ...
                          'max'   ,2 , ...
                          'sliderstep',[0.5 1], ...
                          'value'     ,1 , ...
                          'userdata'  , 0 , ...
                          'callback', CB  , ...
                          'tag'  ,'TAB_SLIDER' , ...
                          'parent' ,H, ...
                          'visible','off');

    for jj = 1 : size(cntr,1)

     CB = [ 'tab_list(' HFrame ','''  cntr{jj,1}   'Button'',1);' ];
 
     ud.cntr{1}(jj) = uicontrol('style','pushbutton', ...
                           'units'   , 'pixels' , ... 
                           'tag'  , ['TAB_' upper(cntr{jj,1}) ] , ...
                           'callback', CB , ...
                           'string',cntr{jj,1}, ...
                           'backgroundcolor', cntr{jj,2} , ...
                           'tooltipstring'  , cntr{jj,3} , ...
                           'parent' ,H, ...
                           'visible','off');
    end

    %----------------------------------------------------
    % Disable Insert-Button first 

    set( ud.cntr{1}(3) , 'enable' , 'off' );

    %----------------------------------------------------

   set(ud.Hframe,'userdata',ud)


   set_prop = set_prop(:,[1 3])';

   msg = tab_list( ud.Hframe , set_prop{:}      , ...
                               'visible' , 'on' , ...
                               VarIn{:,is_font} , ...
                               VarIn{:,is_prop}       );


  out = ud.Hframe;



%******************************************************************************
case 'set_prop'


  fig = get(ud.Hframe,'parent');


  % Store old 'ud'
  ud0 = ud;


  % Set Fields of 'ud' to values of set_prop{:,3}

   is_set = cat(1,set_prop{:,2});
   is_set = find( is_set );

   for ii = is_set(:)'

       prop = set_prop{ii,1};

       switch prop

        case 'vedit' 

         n = max(length(ud.edit{2}),length(set_prop{ii,3})+1);
         ud.edit{2}(2:n) = set_prop{ii,3};


        case { 'cedit' 'cset' 'cins' 'cundo' 'credo' 'cdel'  'creset' }
        % Cntr

          ud.cntr{2}{ii-9} = set_prop{ii,3};

        otherwise
 
          if ~any( strcmp( prop , {'string' 'form' 'title' 'userdata'} ) )

            ud = setfield(ud,set_prop{ii,1},set_prop{ii,3}); 

          end

       end
       % switch

   end
   % ii 

   %--------------------------------------------
   % No  Set | Insert  ==>  No Edit

   if ~any( strcmp(  ud.cntr{2}([2 3]) , 'on' ) )
        ud.cntr{2}{1} = 'off';
   end


    ud.ntab        = round( ud.ntab(1) );
    ud.nrow        = round( ud.nrow(1) );
    ud.ncol        = round( ud.ncol(:) );

    ud.ncol        = [ ud.ncol ; ud.ncol(size(ud.ncol,1))+zeros(ud.ntab,1) ];
    ud.ncol        = ud.ncol( 1 : ud.ntab );

    ud.position    = round( ud.position(1:2) );
    ud.border      = round( ud.border(1)     );

    ud.sliderwidth = round( ud.sliderwidth(1) );

            sloffs = ud.sliderwidth;
   

 %---------------------------------------------------------------
 % Selected Properties

  sel = set_prop(:,[1 2])';

  sel = struct(sel{:});


  if isempty(is_set)

    prop = struct( 'Nothing' , {[]} );

  else

    prop = set_prop(is_set,[1 3])';

    for ii = 1 : size(prop,2)
        prop{2,ii} = { prop{2,ii} };
    end
 
    prop = struct(prop{:});

  end


  
 %---------------------------------------------------------------
 % Set TabularColumns
 %..............................................................

 ntab = ud.ntab + 1;  % Absolut Number of TabularColumns

 if sel.ntab 


  % Delete Columns

  ind = ( ntab+1 : ud0.ntab+1 );

  delete( ud.Htitle(ind)  )
  delete( ud.Hlist(ind)   )
  delete( ud.edit{1}(ind) )


  % Create new Columns


   Nmax = get(ud.Hframe,'max');

     n0 = max(ud.nrow+1,Nmax);

    str = char( 32 * ones(n0,1) );

    val = get(ud.Hlist(1),'value');
    top = get(ud.Hlist(1),'listboxtop');

   lmin = get(ud.Hlist(1),'min');
   lmax = get(ud.Hlist(1),'max');
                             
   tagl = get(ud.Hlist(1),'tag');
   tagt = get(ud.Htitle(1),'tag');
   tage = get(ud.edit{1}(1),'tag');

  HFrame = sprintf( Hform , ud.Hframe );


  for ii = ud0.ntab+2 : ntab

    n = n0 + ( Nmax - n0 ) * ( ii == ntab );

    ud.Htitle(ii) = uicontrol('parent'  , fig    , ...
                              'style'   , 'text' , ...
                              'tag'     , tagt   , ...
                              'visible' , 'off'        );

    CB = [ 'tab_list(' HFrame ',''list'','  sprintf('%.0f',ii) ');' ];

    ud.Hlist(ii) = uicontrol('style'      , 'listbox'  , ...
                             'min'        , lmin       , ...
                             'max'        , lmax       , ...
                             'string'     , str(1:n,:) , ...
                             'value'      , val        , ...
                             'listboxtop' , top        , ...
                             'callback'   , CB         , ...
                             'tag'        , tagl       , ...
                             'parent'     , fig        , ...
                             'visible'    , 'off'            );

    CB = [ 'tab_list(' HFrame ',''edit'','  sprintf('%.0f',ii) ');' ];

    ud.edit{1}(ii) = uicontrol('parent'   , fig        , ...
                               'style'    , 'edit'     , ...
                               'tag'      , tage       , ...
                               'userdata' , { 'char' '' [] }, ...
                               'callback' , CB         , ...
                               'visible'  , 'off'      );

    if length(ud.edit{2}) < ii
      ud.edit{2}(ii) = {'on'};
    end

      
    % Copy Slider to Set it into Front
    
    hsl = ud.Hslider;
    
    ud.Hslider = copyobj( hsl , fig );
    
    delete( hsl );
    
    
 end

  ud.Htitle  = ud.Htitle(1:ntab);
  ud.Hlist   = ud.Hlist(1:ntab);
  ud.edit{1} = ud.edit{1}(1:ntab);
  ud.edit{2} = ud.edit{2}(1:ntab);
    


  % Set String to TabularHandle ud.Hframe

  str = get(ud.Hlist(2:ntab),'string');

  for ii = 1:ud.ntab
    str{ii} = str{ii}(1:Nmax,:);
  end
  
  set(ud.Hframe,'string',str);



 end
 % Set TabularColumns
  

 %-----------------------------------------------------------
 % Set Formats 
 %...........................................................

 if sel.form

  forms = { 'geo1'  ; 'geo2'  ; 
            'time1' ; 'time2' ; 'time3' ; 'time4' ;
            'char'  ; 'none'                        };

  str = prop.form;
  nr2 = [];

  if isnumeric(str) 
      nr2 = str(:);
      str = '';
  elseif ischar(str)
      str = { str };
  elseif iscell(str)
     str = [ str(:)'  {''} ];
   if isnumeric(str{1});
         nr2    = str{1}(:);
   end
   
   if ~isempty(nr2)
    str(1) = [];
   end

   str = str(:)';

   if length(str) == 1
    if ischar(str)
     str = cellstr(str)';
    elseif iscellstr(str)
     str = str(:)';
    else
     str = '';
    end
   else
    for ii = 1 : size(str,2)
     if      iscellstr(str{ii})
      str{ii} = str{ii,1};
     elseif ~ischar(str{ii})
      str{ii} = ''; 
     end
    end
   end

  end
 

 Nstr2 = size(str,2);


  % Fill Missing TabColumns 
  for ii = Nstr2+1 : ud.ntab 
     str{ii} = '';
  end   
        
 % NumberCounter 'nr2'
 if ~isempty(nr2)
   hilf    = ones(Nstr2-size(nr2,1)+1,1);
   hilf(1) = max(nr2)+1;
   nr2(size(nr2,1)+1:Nstr2) = cumsum( hilf(1:Nstr2-size(nr2,1)) );
   nr2 = nr2(1:Nstr2);
 else   
  nr2 = ( 1 : Nstr2 );
 end    
   nr2( find( ( nr2 < 1 )  |  ( ud.ntab < nr2 ) )  ) = [];
        

   form = str;

  for ff = 1:length(nr2)
   form{ff} = deblank(form{ff});
   ok = any(strcmp(form{ff},forms));
   if ~ok
    ok = 1;
    tests = [ -1 0 1 ];
    for jj = 1:length(tests)
     [v,s,m1,m2] = val2str(tests(jj),form{ff});
     ok = ( ok | ( isempty(m1) & isempty(m2) ) );
    end
    if ~ok
     msg = [ msg nl  msg0 ' Invalid Format: '  form{ff} ','  ];
     form{ff} = 'char';
    end
   end
   % ~ok
   ud_ed = get(ud.edit{1}(nr2(ff)+1),'userdata'); % { form  old_string  old_value }
   ud_ed{1} = form{ff};
 
   set(ud.edit{1}(nr2(ff)+1),'userdata',ud_ed);

  end
  % ff
  
 end
 % Set Formats


 %-----------------------------------------------------------
 % Set Strings
 %...........................................................

 if sel.string 

  % The First Column of Cell STR could be the RowNumber 
  % The Second Column could be the ColumnNumber

  str = prop.string;
  nr1 = [];
  nr2 = [];

 if isempty(str)
  str = {''};
 end

  if isnumeric(str) 
      nr1 = str(:);
      str = '';
  elseif ischar(str)
     str = { str };
  elseif iscell(str)
     str = [ str(:)' {'' ''} ];
   if isnumeric(str{1});
         nr1    = str{1}(:);
   end
   if isnumeric(str{2});
         nr2    = str{2}(:);
   end
   
   if ~isempty(nr2)
    str(2) = [];
   end
   if ~isempty(nr1)
    str(1) = [];
   end
 
   % form Character Arrays
   for ii = 1 : size(str,2)
    if      iscellstr(str{ii})
     str{ii} = char(str{ii});
    elseif ~ischar(str{ii})
     str{ii} = ''; 
    end
   end

  end

  Nstr1 = 0;
  % Maximum RowNumber of given Strings
  for ii = 1 : size(str,2) 
    Nstr1 = Nstr1 + ( size(str{ii},1) - Nstr1 ) * ( size(str{ii},1) > Nstr1 );
  end
 

  Nstr2 = size(str,2);


  % NumberCounter 'nr1'
  if ~isempty(nr1)
   hilf    = ones(Nstr1-size(nr1,1)+1,1);
   hilf(1) = get(ud.Hframe,'max')+1;
   nr1(size(nr1,1)+1:Nstr1) = cumsum( hilf(1:Nstr1-size(nr1,1)) );
   nr1 = nr1(1:Nstr1);
   nr1( find( nr1 < 1 ) ) = [];
  else
   nr1 = ( 1 : Nstr1 ) + get(ud.Hframe,'max');
  end

  % NumberCounter 'nr2'
  if ~isempty(nr2)
   hilf    = ones(Nstr2-size(nr2,1)+1,1);
   hilf(1) = max(nr2)+1;
   nr2(size(nr2,1)+1:Nstr2) = cumsum( hilf(1:Nstr2-size(nr2,1)) );
   nr2 = nr2(1:Nstr2);
  else
   nr2 = ( 1 : Nstr2 );
  end

   nr2( find( ( nr2 < 1 )  |  ( ud.ntab < nr2 ) ) ) = [];


  Nstr1 = length(nr1);  % +isempty(nr1);
  Nstr2 = length(nr2);  % +isempty(nr2);


  % Convert with Format
  for ii = 1 : Nstr2 
   ud_ed = get(ud.edit{1}(nr2(ii)+1),'userdata');
   form = ud_ed{1};
   if ~strcmp(form,'char')
    str1 = cell(size(str{ii},1),1);
    for jj = 1:size(str{ii},1)
     [val,str1{jj},msg1,msg2] = val2str(str{ii}(jj,:),form);
     if ~isempty(msg1)  |  ~isempty(msg2)
      str1{jj} = deblank(str{ii}(jj,:));
     end 
    end
    % jj
    str{ii} = char(str1);
   end
   % ~'char'
  end
  % ii

  % Fill Empty Rows  
  for ii = 1:Nstr2
    str{ii}(size(str{ii},1)+1:Nstr1,:) = char( 32 + ...
                         zeros( Nstr1-size(str{ii},1) , size(str{ii},2) ) );
  end
  
   if isempty(nr1)
      mnr1 = 0;
   else
      mnr1 = max(nr1);
   end
  
   Nmax0 = get(ud.Hframe,'max');          % Old Number of Entries
   Nmax1 = max( mnr1 , Nmax0 );       % New Number of Entries

   Noffs = max( ud.nrow+1 , Nmax1 ) - Nmax1;  

   sl_max = Nmax1 - ud.nrow + 2;
   sl_max = sl_max + (2-sl_max)*( sl_max < 2 );

   sl_stp = [ 1/sl_max  max( 2/sl_max , round(0.1*sl_max)/sl_max ) ];
   sl_stp(2) = 1 + (sl_stp(2)-1) * (sl_stp(2)>sl_stp(1)); 

   top  = get(ud.Hlist(1),'listboxtop');

   val1 = 1 + sl_max - top;
   val2 = val1 + (  1     - val1 ) * ( val1 <  1     );
   val2 = val2 + ( sl_max - val2 ) * ( val2 > sl_max );

   top  = 1 + sl_max - val2;
 
   set( ud.Hslider , 'max'        , sl_max , ...
                     'sliderstep' , sl_stp , ...
                     'value'      , val2         )

  %--------------------------------------------------------
  % Set given Strings

  for ii = 1 : length(nr2) 

    str0 = get(ud.Hlist(nr2(ii)+1),'string');
    str0 = str0(1:min(Nmax0,size(str0,1)),:);
    str1 = str{ii};

    % Fill Columns
    str0 = [ str0 char( 32*ones(size(str0,1),size(str1,2)-size(str0,2)) ) ]; 
    str1 = [ str1 char( 32*ones(size(str1,1),size(str0,2)-size(str1,2)) ) ];

    % Fill Rows;

    nf   =  Nmax1 - Nmax0 + Noffs * ( nr2(ii) < ud.ntab );
 
    str0 = cat( 1 , str0 , ...
                    char( 32 * ones( nf , size(str0,2) ) ) );

    if ~isempty(str0)
      if isempty(str1) & ~isempty(nr1)
        str0(nr1,:) = char( 32*ones(1,size(str0,2)) );
      else
        str0(nr1,:) = str1;
      end
    end

    set( ud.Hlist(nr2(ii)+1) , 'string' , str0 )

  end
  
  %--------------------------------------------------------
  % Set NumberString

   str00 = zeros(0,3);

   if Nmax1 > 0
     form0 = '''%3.0f'';';
     str00 = sprintf(form0,(1:Nmax1));
     str00 = cat( 2 , '[' , str00(1:(end-1)) , ']' );
     str00 = eval(str00);
   end
   
   str00 = cat( 1 , str00 , char( 32*ones(Noffs,size(str00,2)) ) );


   set( ud.Hlist(1) , 'string' , str00 )

  %--------------------------------------------------------

   val = getval( ud.Hlist(1) , Nmax1 );

   set( ud.Hlist , 'value'      , val , ...
                   'listboxtop' , top       );

   set( ud.Hframe , 'max'   , Nmax1 , ...
                    'value' , val          );

   % Fill UserData  

   ud_userdata          = cell(Nmax1,1);
   if ~isempty(ud.userdata)
     ud_userdata(1:Nmax0) = ud.userdata(1:Nmax0);
   end
   ud.userdata          = ud_userdata;

 end
 % String

 


 %----------------------------------------------------------------
 % Set Title
 %...............................................................

 if sel.title

  str = prop.title;
  nr2 = [];

  if isnumeric(str) 
      nr2 = str(:);
      str = '';
  elseif ischar(str)
      str = { str };
  elseif iscell(str)
     str = [ str(:)'  {''} ];
   if isnumeric(str{1});
         nr2    = str{1}(:);
   end
   
   if ~isempty(nr2)
    str(1) = [];
   end

   str = str(:)';

   if length(str) == 1
    if ischar(str)
     str = cellstr(str)';
    elseif iscellstr(str)
     str = str(:)';
    elseif ~iscellstr(str{1})
     str = '';
    end
   else
    for ii = 1 : size(str,2)
     if ~( ischar(str{ii}) | iscellstr(str{ii}) )
      str{ii} = ''; 
     end
    end
   end

  end
 
 Nstr2 = size(str,2);

  % Fill Missing TabColumns 
  for ii = Nstr2+1 : ud.ntab 
   str{ii} = '';
  end
     
 % NumberCounter 'nr2'
 if ~isempty(nr2)
   hilf    = ones(Nstr2-size(nr2,1)+1,1);
   hilf(1) = max(nr2)+1;
   nr2(size(nr2,1)+1:Nstr2) = cumsum( hilf(1:Nstr2-size(nr2,1)) );
   nr2 = nr2(1:Nstr2);
 else
  nr2 = ( 1 : Nstr2 );
 end

   nr2( find( ( nr2 < 1 )  |  ( ud.ntab < nr2 ) ) ) = [];
  
  
   for ii = 1:length(nr2)
    set(ud.Htitle(nr2(ii)+1),'string',str{ii});
   end


   ud.title = get(ud.Htitle,'string');

 end
 % Set Title


 %----------------------------------------------------------------
 % Set UserData
 %...............................................................

 if sel.userdata

  val = prop.userdata;
  nr1 = [];
  
  Nmax = get(ud.Hframe,'max');

  if  size(val,1) == 1   &  ...
      size(val,2) == 2   &  ...
      isnumeric(val{1})

     nr1 = val{1}(:);

     val = val{2}(:);

   else
     val = val(:);
   end
 
   if isempty(nr1)
     nr1 = ( 1 : length(val) )';
   end

   nr1 = round(nr1);

   ii = find( ( nr1 < 0 )  |  ( Nmax < nr1 ) );

   nr1(ii) = [];
   val(ii) = [];

   ud.userdata(nr1) = val;

 end
 % Set UserData       


 %-----------------------------------------------------------------------
 % Check, if we need the SLIDER
 %.......................................................................

     is_win = strcmp( upper(computer) , 'PCWIN' );
 
     is_slider = ( ( get(ud.Hframe,'max') > ud.nrow ) | is_win );
     slider_on = get(ud.Hslider,'userdata');

     set(ud.Hslider,'userdata',is_slider);

  set_slider = +1 * ( is_slider & ~slider_on) + ...
               -1 * (~is_slider &  slider_on) ;

  % +1 .. Set Slider ON   to hide the ListBoxSlider
  % -1 .. Set Slider OFF

  if ~( set_slider == 0 )

    uni = get(ud.Hlist(ntab),'units');

    set(ud.Hlist(ntab),'units','pixel')

    pos    = get(ud.Hlist(ntab),'position');

    pos(3) = pos(3) + set_slider*sloffs;

    set(ud.Hlist(ntab),'position',pos);

    set(ud.Hlist(ntab),'units',uni);

  end 
 


 %---------------------------------------------------------------
 % Set Position if ntab, nrow, ncol, position, cntr  or  Font-Property
 %...............................................................


 if  ( any( [sel.ntab   sel.nrow     sel.ncol           ...
             sel.border sel.position sel.sliderwidth    ...
             sel.vspace sel.vedit    sel.cedit          ...
             sel.cset   sel.cundo    sel.credo sel.cdel sel.title ] )  | ...
       ~isempty(is_font) );

 
   % Set the FigureChildren into correct order for Visibility
   %  and store UNITS of UIControls  
   %  and set   UNITS to PIXELS      !!!!!!!!!!!!!!!!!

    HH  = flipud([ ud.Hframe  ud.Hlist ud.Htitle ud.edit{1} ...
                   ud.Hslider ud.cntr{1}(:)'                      ]');
 
    uni = get(HH,'units');

    set(HH,'units','pixels');


  if ~isempty(is_font)
   try  
     set(ud.Hframe,VarIn{:,is_font});
   catch
     is_font = [];
   end
  end

   ppi = get(0,'screenpixelsperinch');


   fontunits = get(ud.Hframe,'fontunits'); % !!!!!!!!!!!!!!!!!
   fontsize  = get(ud.Hframe,'fontsize');

   set(ud.Hframe,'fontunits',fontunits);    % !!!!!!!!!!!!!!!!!!
 

  % Use DummyTextHandle to Determine TextWidth 

   set(ud.Htext,  ...
          'units'      , 'pixels'  , ...
          'position'   , [ 1  1  0 ] , ...
          'string'     , '' , ...
          'visible'    , 'off' , ...
         'FontName'    , get( ud.Hframe , 'FontName'  ) , ...
         'FontUnits'   , get( ud.Hframe , 'FontUnits' ) , ...
         'FontSize'    , get( ud.Hframe , 'FontSize'  ) , ...
         'FontAngle'   , get( ud.Hframe , 'FontAngle' ) , ...
         'FontWeight'  , get( ud.Hframe , 'FontWeight' )       );


   ext = zeros(2,4);

   for ii = [ 1  2 ]
       set( ud.Htext , 'string' , char('H'*ones(ii,ii)) );
       ext(ii,:) = get(ud.Htext,'extent');
   end

   mm = ext(2,[3 4]) - ext(1,[3 4]);
   nn = ext(1,[3 4]) - mm;

   hh  = mm(2);  % Single UIControl
   hht = mm(2);

   for ii = 1 : size(ud.title,1)

     set(ud.Htext,'string',ud.title{ii});

     ext = get(ud.Htext,'extent');

     hht = max(hht,ext(4));

   end
   
   set(ud.Htext,'string','');


   % Determine Positions

   ncol = [ 3+is_win ; ud.ncol ];    % First TabColumn == NumberTab, 3+1 Characters

   ww = ceil( ncol * mm(1) + nn(1) );   
   ww = ww + ud.lineoffset(1) + 2*ud.borderwidth + sloffs;   % Width's
   ll = cumsum([ 1 ; ww-sloffs ]);         % Left's
  
   pp = ud.position(:)';
   
   is_cntr = any(strcmp('on',ud.cntr{2}));
   is_edit = any(strcmp('on',ud.edit{2}(2:length(ud.edit{2})))) & ...
                 strcmp('on',ud.cntr{2}(2));  % & SET 'on'

   hhc = ceil(is_cntr * (hh+2*ud.borderwidth+0) );   % Cntr
   hhe = ceil(is_edit * (hh+2*ud.borderwidth+2) );   % Edit
   hhl = ceil(ud.nrow * (hh+ud.lineoffset(2))+2*ud.borderwidth);   % ListBox


   % FigurePosition in Pixels

    figuni = get(fig,'units');
             set(fig,'units','pixels');
    figpos = get(fig,'position');
             set(fig,'units',figuni);
   

   fpos = zeros(1,4);
   fpos(3) = max(ll) + sloffs + 2*ud.border;
   fpos(4) = (hhc+hhe+hhl+hht+(1+is_cntr+is_edit)*ud.vspace) + 2*ud.border;
   fpos([1 2]) = ( pp([1 2]) + 1 ) .* ( pp([1 2]) >= 0 ) + ...
                 ( figpos([3 4]) - fpos([3 4]) + 1 + pp([1 2]) ) .* ...
                 ( pp([1 2]) < 0 );

   set( ud.Hframe , 'units' , 'pixels' , 'position' , fpos );

   set(ud.Haxe,'units'    , 'pixels' , ...
               'position' , fpos     , ...
               'visible'  , 'off'                              )

   
   pp = fpos([1 2])+ud.border;

   if is_cntr

      jj = find(strcmp('on',ud.cntr{2}));
     Njj = length(jj);
 
       ppc = linspace(1,max(ll)+sloffs,Njj+1);   % Width's
       wwc = diff(ppc);
 
     for kk = 1 : Njj

      set(ud.cntr{1}(jj(kk)),'units','pixels','position', ...
         [  pp+[ ceil(ppc(kk)) 0 ]  floor(wwc(kk))   hhc ])                 

     end

   end
  

  for ii =  1 : ntab;

    set(ud.edit{1}(ii),'units','pixels','position', ...
           [ pp+[ ll(ii)  hhc+is_cntr*ud.vspace   ]  ww(ii)-sloffs   hhe+(hhe==0) ]) 

    set(ud.Hlist(ii),'units','pixels','position', ...
           [ pp+[ ll(ii)  hhc+hhe+(is_edit+is_cntr)*ud.vspace ] ...  
             ww(ii)-sloffs*(ii==ntab)*(~(is_slider | is_win))   hhl ])
             
    set(ud.Htitle(ii),'units','pixels','position', ...
           [ pp+[ ll(ii)  hhc+hhe+hhl+(1+is_edit+is_cntr)*ud.vspace ] ...
             ww(ii)-sloffs  hht ])

  end

  pos_sl = zeros(1,4);
  pos_sl(1) = pp(1) + max(ll) - ud.borderwidth*is_win;
  pos_sl(2) = pp(2) + hhc+hhe+(is_edit+is_cntr)*ud.vspace + ud.borderwidth*is_win;
  pos_sl([3 4]) = [ sloffs  hhl-2*ud.borderwidth*is_win ];
  
  set(ud.Hslider,'units','pixels','position',pos_sl);

   
 end
 % Set Position


 %-------------------------------------------------------------------
 % Font and Other Properties
 %...................................................................

  sets = { 'off' ; 'on' };  

 %--------------------------------------------------------------------
 % Not allowed Property's 

  notallwd = {'value'
              'callback' 
              'enable'     };

  for ii = 1 : size(notallwd,1)
              jj   = find( strcmp( notallwd{ii} , VarIn(1,is_prop) ) );
     is_prop( jj ) = [];
  end


  % for FRAME, EDIT, SLIDER

  notallwd = {'min'
              'max'  
              'userdata'};

  is_prop0 = is_prop;

  for ii = 1 : size(notallwd,1)
               jj   = find( strcmp( notallwd{ii} , VarIn(1,is_prop0) ) );
     is_prop0( jj ) = [];
  end


  % for CONTROL

  notallwd = {'backgroundcolor' };

  is_prop1 = is_prop;

  for ii = 1 : size(notallwd,1)
               jj   = find( strcmp( notallwd{ii} , VarIn(1,is_prop1) ) );
     is_prop1( jj ) = [];
  end


  
 %--------------------------------------------------------------------
 % Check other Property's

 if ~isempty(is_prop) 
   try
     set( ud.Hlist(1) , VarIn{:,is_prop} );
   catch
     msg = [ msg nl  msg0  ' Error in UIControlProperty''s:'  nl ...
              lasterr   nl ];
     is_prop    = [];
     is_prop0   = [];
   end

 end


  set( ud.Hframe , VarIn{:,is_font}  , ...
                   VarIn{:,is_prop0} , ...
                   'userdata' , ud        );

  vis_off = strcmp( get(ud.Hframe,'visible') , 'off' );


  setv = sets{ 2 - vis_off };


  set( ud.Hlist  , VarIn{:,is_font} , ...
                   VarIn{:,is_prop} , ...
                   'visible' , setv        );
 
  set( ud.Htitle , VarIn{:,is_font} , ...
                   VarIn{:,is_prop} , ...
                   'visible' , setv       );



 set_edit = ud.edit{2};

 if all(strcmp('off',ud.edit{2}(2:length(ud.edit{2}))))

   set_edit(1) = {'off'};

 end

 if strcmp('off',ud.cntr{2}(2))  % SET off

   set_edit(1:length(set_edit)) = {'off'};

 end

 for jj = 1 : ntab

   set( ud.edit{1}(jj) , VarIn{:,is_font}         , ...
                         VarIn{:,is_prop0}        , ...
                         'visible' , set_edit{jj}      );

 end

  set( ud.Hslider ,  VarIn{:,is_prop0}                         , ...
                    'visible' , sets{ 1 + ( is_slider & (~vis_off) ) }  );


 for jj = 1 : size(ud.cntr{1},1)

  % For Controls only FontProperty's !!!

  set( ud.cntr{1}(jj) , VarIn{:,is_font}    , ...
                        VarIn{:,is_prop1}   , ...
                       'visible' , ud.cntr{2}{jj}  );
 
 end

 %------------------------------------------------
 % Global Visibility

 if vis_off

    tab_list(ud.Hframe,'visible','off');

 end

 %-------------------------------------------------------------------
 % Fill Strings Right  if STRING, NTAB, NROW
 %...................................................................

 if sel.string | sel.ntab | sel.nrow

    Nmax0 = get(ud.Hframe,'max');

   ud_str = cell(1,ntab);

   for ii = 1 : ntab

     str   = get(ud.Hlist(ii),'string');

     Nmax1 = size(str,1);

     Nmax2 = max(Nmax0,Nmax1);

     Nmax3 = ( ud.nrow + 1 ) * ( ii <= ud.ntab );

        nf = max(Nmax2,Nmax3) - Nmax1;

       str = cat( 1 , str , ...
                      char( 32 * ones( nf , size(str,2) ) ) );

       set( ud.Hlist(ii) , 'string' , str );

       ud_str{ii} = str(1:Nmax0,:);

   end

   ud.string = ud_str(2:ntab);

   top = get(ud.Hlist(1),'listboxtop');

   top = top + ( Nmax0-ud.nrow+1 - top ) * ( top > Nmax0-ud.nrow+1 );
   top = top + ( 1 - top ) * ( top < 1 );
  
   val = getval( ud.Hlist(1) , Nmax0 );

   set(ud.Hlist,'value'     ,val , ...
                'listboxtop',top       );

   set(ud.Hframe,'value'    , val , ...
                 'userdata' , ud          );

 end

 
 setenable( ud.Hframe );


%***********************************************************************
case 'list'

  Nmax = get(ud.Hframe,'max');

  if Nmax == 0
    return
  end

  h = ud.Hlist(VarIn{1});

  top = get(h,'listboxtop');

  top_max = Nmax - ud.nrow + 2;
  top_max = top_max + ( 1 - top_max ) * ( top_max <= 0 );

  top = top + ( top_max - top ) * ( top > top_max );

  val = get(h,'value');
  if ~isempty(val)
     val = val + ( Nmax - val ) .* ( val > Nmax );
  end

  val = getval( h , Nmax , val );

  set(ud.Hframe,'value',val)
  set(ud.Hlist ,'value',val,'listboxtop',top)
 
  set(ud.Hslider,'value',1+get(ud.Hslider,'max')-top);

  callback_ok = 1;
  callback_in = { 'TabList'  1  val };

  %---------------------------------------------------------
  % OPEN ==> SetEdit  

  sel = get( get(ud.Hframe,'parent') , 'selectiontype' );

  if strcmp( sel , 'open' )

     tab_list(ud.Hframe,'setedit',val(1));

  else

     setenable( ud.Hframe );

  end

%***********************************************************************
case 'slider'
 
  val = round(get(ud.Hslider,'value'));
  top = 1+get(ud.Hslider,'max')-val;

  if top > get(ud.Hframe,'max')
     return
  end
 
  set(ud.Hlist  ,'value'     , top , ...
                 'listboxtop', top   )
 
  set(ud.Hslider,'value',val)

  set(ud.Hframe,'value',top)

  callback_ok = 1;
  callback_in = { 'TabSlider'  1  top };

%***********************************************************************
case 'edit'

  h = ud.edit{1}(VarIn{1});

  ud_ed = get(h,'userdata');
  str   = get(h,'string');


  form = ud_ed{1};

  if strcmp(form,'char')

    val = str;

  else

    [val,str,msg1,msg2] = val2str(str,form);
  
    if ~isempty(msg1)  |  ~isempty(msg2)
      errordlg('Wrong Input','Warning')
      str = ud_ed{2};
      val = ud_ed{3};
    end 
  
  end
  
  ud_ed{2} = rmblank(str,2);
  ud_ed{3} = val;

  set(h,'string',ud_ed{2},'userdata',ud_ed)

  if ~( VarIn{1} == 1 )
     return
  end

  %--------------------------------------
  % Enable Insert-Button

  setenable( ud.Hframe );

 
%***********************************************************************
case 'editbutton'

  val = get(ud.Hframe,'value');

  if ~isempty(val)

    tab_list(ud.Hframe,'setedit',val(1));

  end


%***********************************************************************
case 'setedit'

  val = VarIn{1};

  if isempty(val)
    return
  end

  if ~( isnumeric(val)  &  ( prod(size(val)) == 1 )  ) 
   msg = [ msg0  'Value for ListNumbers to Edit must be Numeric.'];
   return
  end

  val = val(:);

  Nmax = get(ud.Hframe,'max');

  val( find( val < 1  | val > Nmax ) )  = [];

  if isempty(val)

     return

  end


   val = val(1);

   for ii = 1 : ud.ntab+1

     str = get(ud.Hlist(ii),'string');
     str = rmblank(str(val,:),2);

     if strcmp( get(ud.edit{1}(ii),'style'),'popupmenu')

        txt = cellstr(get(ud.edit{1}(ii),'string'));
        cmp = strcmp( txt , str );

        if ~any(cmp)
            cmp = strwcmp(txt,['*' str]);
            if ~any(cmp)
                cmp = strwcmp(txt,[str '*']);
                if ~any(cmp)
                    cmp = strwcmp(txt,['*' str '*']);
                end
            end
         end

         if ~any(cmp)
             txt = cat(1,txt,{str});
             set( ud.edit{1}(ii),'string' , txt , ...
                                 'value'  , size(txt,1) );
         else
             nr = find(cmp);
             set( ud.edit{1}(ii),'value',nr(1));
         end

     else

        ud_ed    = get( ud.edit{1}(ii) , 'userdata' );
        ud_ed{2} = str;

        if ii == 1
           ud_ed{3} = val;
        else
           ud_ed{3} = str;
        end
        
        set( ud.edit{1}(ii) , 'string'   , str , ...
                              'userdata' , ud_ed );
     end
   end

  %--------------------------------------------------
  % Enable Insert-Button

    setenable( ud.Hframe );


%***********************************************************************
case { 'setbutton'  'insbutton' }

  ntab = ud.ntab + 1;

  Nmax = get(ud.Hframe,'max');

  ud_ed = get(ud.edit{1}(1),'userdata');

  [nr,str] = val2str(get(ud.edit{1}(1),'string'),ud_ed{1});

  ok = ~isempty(nr);
  if ok
     nr = nr(1);
     ok = ( isfinite(nr) & ( nr > 0 ) );
  end

  if ok
     nr = ceil(nr);
  else 
     nr = Nmax+1;
  end

  str = cell(1,ntab);
  str{1} = nr;

  for ii = 2 : ntab
      str{ii} = get(ud.edit{1}(ii),'string');
      if strcmp(get(ud.edit{1}(ii),'style'),'popupmenu')
          val = get(ud.edit{1}(ii),'value');
          str{ii} = cellstr(str{ii}); 
          str{ii} = str{ii}{val};
      end
  end

  fcn_ok = isempty(ud.setquestfcn);

  if ~fcn_ok
        
       msg = ''; 

       try
         fcn_ok = feval(ud.setquestfcn{:},'TabSetQuest',1,str);
       catch
         msg = lasterr;
       end

       if ~isempty(msg);

          msg = [ msg0 ' Error in SetQuestFcn. ' nl ...
                  msg nl ];

          fprintf([ msg  nl ])   

       end

  end

  %----------------------------------------------------------
  % Check fcn_ok:  0 | 1
  %                str:  { Nr String1 String2 ... }  

  ok = ( isnumeric(fcn_ok)  &  ( prod(size(fcn_ok)) == 1 ) );

  if ok
     %----------------------
     %  0 | 1

     ok = isequal( fcn_ok , 1 );

  else
     %----------------------
     %  str

     ok = iscell(fcn_ok);

     if ok
     
        fcn_ok = fcn_ok(:)';
     
        ok = ( ( size(fcn_ok,1) == 1           ) & ...
               ( size(fcn_ok,2) == size(str,2) ) & ...
               isnumeric(fcn_ok{1})              & ...
               ( prod(size(fcn_ok{1})) <= 1 )    & ...
               iscellstr(fcn_ok(2:size(fcn_ok,2)))    );

     end
   
     if ok
       
         str = fcn_ok;
       
         if isempty(str{1})
            str{1} = Nmax+1;
         end

     else

        msg = [ msg0 ' Invalid Data in SetQuestFcnReturn.' ];
        fprintf([msg nl ])   

     end

  end
       
  if ~ok
     return
  end

  %----------------------------------------------------------
  % Check for Insert

  str0 = str;

  if strcmp( action , 'insbutton' )  & ( str{1} > get(ud.Hframe,'max') )

    % Append, NOT Insert !!!
    
        action = 'setbutton';

  end

  %----------------------------------------------------------
  % Insert

  if strcmp( action , 'insbutton' )

      % Insert String !!!

      nr      = str0{1};
 
      str     = cell(1,ntab);

      str{1}  = ( 1 : Nmax+1 )';
 
      ind = ( 1 : Nmax+1 );

      ind((nr+1):(Nmax+1)) = ind((nr+1):(Nmax+1)) - 1;     
     
      for ii = 2 : ntab
        str1    = get(ud.Hlist(ii),'string'); 
        str{ii} = cellstr(str1(ind,:));
        str{ii}(nr) = str0(ii);    
      end

      % Insert UserData !!!

      ud_userdata = ud.userdata(ind(:));
      ud_userdata(nr) = {[]};

  end

  tab_list(ud.Hframe,'set_prop','string',str);

  tab_list(ud.Hframe,'value',str0{1});

  %--------------------------------------
  % Reset Inserted UserData

  if strcmp( action , 'insbutton' )

      ud = get(ud.Hframe,'userdata');
      ud.userdata = ud_userdata;
           set( ud.Hframe , 'userdata' , ud  )

  end

  %--------------------------------------
  % Check Enability of Insert-Button

    setenable( ud.Hframe );

  %--------------------------------------

  callback_ok = 1;

  CallAct = cat( 2 , 'Tab' , upper(action(1)) , lower(action([2 3])) );

  callback_in = { CallAct  1  str0 };



%***********************************************************************
case { 'delbutton'  'resetbutton' }


  Nmax = get(ud.Hframe,'max');

  if Nmax == 0
     return
  end

  if strcmp( action , 'delbutton' )

     val = get(ud.Hframe,'value');

     if isempty(val)
       
       return

     end
 
     val = sort(val(:)');

  else

     val = ( 1 : Nmax );

  end


    if ~isempty(ud.delquestfcn)
        
       msg = ''; 
       ok  = 0;

       try
         ok = feval(ud.delquestfcn{:},'TabDelQuest',1,val);
       catch
         msg = lasterr;
       end

       if ~isempty(msg)

           msg = [ msg0 ' Error in DelQuestFcn. ' nl ...
                   msg nl ];
           
           fprintf([ msg  nl ])   

       end

    else

      if isequal( val , ( 1 : Nmax ) );
 
         txt = 'Delete ALL ListEntry''s?';

      else
   
         txt = ['Delete ListEntry''s?' nl ...
                ' ( ' sprintf('%.0f ',val) ')' ];

      end

      ok = questdlg( txt , 'Stop' , ...
                     'Delete','Cancel','Cancel');

      ok = strcmp(ok,'Delete');

    end          
 
    if ~ok 

       return

    end


    callback_ok = 1;

    CallAct = cat( 2 , 'Tab' , upper(action(1)) , ...
                   strrep(lower(action(2:end)),'button','')  );

    callback_in = { CallAct  1  val };


    tab_list( ud.Hframe , 'delete' , val );



%***********************************************************************
case 'delete'
  
  val = VarIn{1};

  if isempty(val)
     return
  end

  if isequal(val,0) | isequal(val,-1)

     % Delete TabList
  
       hh = cat( 2 , ud.Haxe , ud.Htext   , ...
                     ud.Htitle , ud.Hlist, ud.Hslider , ...
                     ud.edit{1}(:)' , ud.cntr{1}(:)'        );

       for h = hh
          try, delete(h), end
       end

       if  val == 0
       % call from outside, not by DeleteFcn
         set(ud.Hframe,'deletefcn','');
         delete(ud.Hframe);
       end

       return

  end

  %-----------------------------------------------------

  ok = all( isnumeric(val) );
  if ok
     ok = all(isfinite(val));
     if ok
        ok = all( ( mod(val,1) == 0 )  &  ( val >= 0 ) );
     end
  end

  if ~ok
   msg = [ msg0  'Value for ListNumbers to Delete must be Numeric.'];
   return
  end


  val = val(:);

  Nmax = get(ud.Hframe,'max');

  val( find( val < 1  | val > Nmax ) )  = [];

  if isempty(val)

     return

  end

    ntab    = ud.ntab + 1;

    nr      = ( 1 : Nmax )';
    nr(val) = [];

    str    = cell(1,ntab);
    str{1} = ( 1 : Nmax-length(val) )';
     
    for ii = 2 : ntab
      str1    = get(ud.Hlist(ii),'string'); 
      str{ii} = str1(nr,:);   
    end
   
    ud.userdata = ud.userdata(nr);

    Nmax =  Nmax - length(val);

    set( ud.Hframe , 'max'      , Nmax , ...
                     'userdata' , ud  )

    tab_list( ud.Hframe , 'set_prop' , 'string' , str );

        

%***********************************************************************
case { 'undo'  'redo' 'undobutton'  'redobutton' }

 warndlg('This Function is under Construction!','Sorry')



%***********************************************************************
case 'value'

  val = VarIn{1};
  if ~isnumeric(val)
   msg = [ msg0  'Value must be Numeric.'];
   return
  end

  val = val(:);

  Nmax = get(ud.Hframe,'max');

  val = getval( ud.Hlist(1) , Nmax , val );

  if isempty(val)

     set(ud.Hframe,'value',val);
     set(ud.Hlist ,'value',val);

  else

     top = val(1);  % get(ud.Hlist(1),'listboxtop');

     top_max = Nmax - ud.nrow + 1;
     top_max = top_max + ( 1 - top_max ) * ( top_max <= 0 );

     top = top + ( top_max - top ) * ( top > top_max );
 
     set(ud.Hframe,'value',val)
     set(ud.Hlist ,'value',val,'listboxtop',top)

     set(ud.Hslider,'value',1+get(ud.Hslider,'max')-top);
  
  end

  setenable( ud.Hframe );


%***********************************************************************
case 'val2str'


 nr   = VarIn{1}(:);
 val0 = VarIn{2}(:);

 if ~isnumeric(nr)
  msg = [' Input for TabularColumns must be Numeric.'  nl ];
 end

 if ~isnumeric(val0) & ~iscell(val0) & ~iscellstr(val0)  &  ~ischar(val0)  
  msg = [' Input for Values/Strings must be an numeric, char or CellArray.'  nl ];
 end

 if iscell(val0) & ~iscellstr(val0)
   try
     ok = isnumeric(cat(1,val0{:}));
   catch
     ok = 0;
   end
   if ~ok
    ok = 1;
    for ii = 1 : length(val0)
     ok = ( ok & ( isnumeric(val0{ii}) | ischar(val0{ii}) ) );
    end
    if ~ok
     msg = [' CellElements in CellArray of Values must be Numeric or Char.'  nl ];
    end
   end 
 end


 if ~isempty(msg)
  msg = [ msg0   msg ];
  return
 end

 if ischar(val0)
   val0 = cellstr(val0);
 end
 if isnumeric(val0)
   val0 = num2cell(val0);
 end


 if isempty(nr); nr = 1; end


 if length(nr) == 1
  nr = nr(1)+zeros(size(val0,1),1);
 end
 if length(val0) == 1
   val1 = cell(length(nr),1); val1(:) = val0;
   val0 = val1;
 end
 
 if length(nr) ~= length(val0)
  msg = [ msg0  ' Length of TabularColumns and Values/Strings  must be the same.' ];   
  return
 end



 out = cell(length(nr),5);  % { Ok  Value  String  Msg1  Msg2}

 out(:,1  ) = { 0};
 out(:,2  ) = {[]};
 out(:,3:5) = {''};

 nr = round(nr);

 for ii = 1:length(nr)

  if ~( ( 1 <= nr(ii) )    &   ( nr(ii) <= ud.ntab ) )
   % Out of TabularColumnRange
   
      out{ii,1} = NaN;

  else

    
    ud_ed = get(ud.edit{1}(nr(ii)+1),'userdata');
 
    form = ud_ed{1};

    if strcmp(form,'char')
     
     jj = 2 *  isnumeric( val0{ii} )  + ...
          3 *     ischar( val0{ii} )  ;

     if jj
       out(ii,[1 jj]) = [ {1}  val0(ii) ];
     end

    else

     if ~isempty(val0{ii})

       if ischar(val0{ii}) 
          val0{ii} = val0{ii}(1,:);
       end

      [val,str,msg1,msg2] = val2str(val0{ii},form);
  
       if isempty(msg1) & isempty(msg2)
        out(ii,1:3) = { 1  val  str };
       end 

        out(ii,4:5) = { msg1  msg2 };

     end


   end

  end

 end
 % ii
   

%***********************************************************************
case 'button'

% New Button to ControlBar


  ok = isempty(VarIn);

  if ~ok

    if ( mod(size(VarIn,1),2) ~= 0 )
      msg = [ 'Additional Inputs must contain UIControl-Property-Value-Pairs.' ];
    else 
      VarIn = reshape(VarIn,2,size(VarIn,1)/2);
      if ~iscellstr(VarIn(1,:))
        msg = [ 'Additional Inputs must contain UIControl-Property-Value-Pairs.' nl ...
                'Properties must be Strings.' ];
      end
    end

  end

  if ~isempty(msg)
     msg = [ msg0  msg ];
     return
  end
  
  h0 = ud.cntr{1}(1);

  h = uicontrol( 'parent'          , get(h0,'parent')          , ...
                 'style'           , 'pushbutton'              , ...
                 'backgroundcolor' , get(h0,'backgroundcolor') , ...
                 'foregroundcolor' , get(h0,'foregroundcolor') , ...
                 'fontname'        , get(h0,'fontname')        , ...
                 'fontunits'       , get(h0,'fontunits')       , ...
                 'fontsize'        , get(h0,'fontsize')              );

  if ~isempty(VarIn)

    try

      set(h,VarIn{:});
           
    catch

       msg = [ msg0  'Invalid Inputs.' nl lasterr ];
       delete(h) 
       return

    end

  end      

   ud.cntr{1} = cat(1,ud.cntr{1},h);
   ud.cntr{2} = cat(1,ud.cntr{2},{'on'});

   set(ud.Hframe,'userdata',ud);

   msg = tab_list(ud.Hframe,'position',ud.position);
 
   if ~isempty(msg)
     msg = [ msg0 'Error call TAB_LIST.' nl  msg ];
     delete(h)
     return
   end

   out = h;


%***********************************************************************
case 'visible'

   
  sets = VarIn{1};

  ok = ( ischar(sets) & ~isempty(sets) & ...
         ( prod(size(sets)) == size(sets,2) ) );
  if ok
     ok = any(strcmp(sets,{'on'  'off'}));
  end

  if ~ok

    msg = [ msg0  'Input for Visible must be ''on'' or ''off''.'];
    return

  end

  if strcmp(sets,'off')

       hh = cat( 2 , ud.Hframe , ud.Haxe , ud.Htext , ...
                     ud.Htitle , ud.Hlist, ud.Hslider  );

       set(hh,'visible','off');
  
       set( ud.edit{1} , 'visible' , 'off' );
       set( ud.cntr{1} , 'visible' , 'off' );


  else

       set(ud.Hframe,'visible','on');

       tab_list( ud.Hframe , 'ntab' , ud.ntab );

  end


%***********************************************************************
case 'enable'

   
  sets = VarIn{1};

  ok = ( ischar(sets) & ~isempty(sets) & ...
         ( prod(size(sets)) == size(sets,2) ) );
  if ok
     ok = any(strcmp(sets,{'on'  'off'}));
  end

  if ~ok

    msg = [ msg0  'Input for Enable must be ''on'' or ''off''.'];
    return

  end


       hh = cat( 2 , ud.Htitle , ud.Hlist, ud.Hslider  );

       set(hh,'enable',sets);
  
       set( ud.edit{1} , 'enable' , sets );
       set( ud.cntr{1} , 'enable' , sets );


   if strcmp( sets , 'off' )
      return
   end

  %--------------------------------------
  % Check Enability of Buttons

   setenable( ud.Hframe );


end
% switch action



%***********************************************************************
% CallBack 

if callback_ok  &  ~isempty(ud.callback)

 msg = '';

 try
   feval( ud.callback{:} , callback_in{:} );
 catch
   msg = lasterr;
 end
 
 if ~isempty(msg)

    msg = [ msg0 ' Error in UserCallBackString. ' nl ...
                  msg nl ];

    fprintf([ msg  nl ]);

 end

end



if Nout == 0
 out = [];
end


%***********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function val = getval(h,N,val)

if nargin < 3
   val = get(h,'value');
end

if ~isempty(val)

    val = val(:);

    val( find( val > N ) ) = [];

    if ~isempty(val)

       val( find( val < 1 ) ) = [];

       if prod(size(val)) > 1
          
          val = sort(val);

          val( find( diff(val) == 0 ) + 1 ) = [];

       end

       return

    end

end
   
% Value is EMPTY, Check with Min & Max of ListBoxHandle

if isempty( get(h,'string') )  | ...
   ( ( get(h,'max') - get(h,'min') ) > 1 ) 

   return

end

val = 1;

%***********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function setenable( hf );

% Sets Enability of Edit-, Insert-, Delete-, Reset-Button

ud = get( hf , 'userdata' );

h = ud.cntr{1}([1 3 6 7]);

Nmax = get(ud.Hframe,'max');


if ( Nmax == 0 ) | ...
   strcmp( get(ud.Hlist(1),'enable') , 'off' )

   set( h , 'enable' , 'off' );

   return

end

%--------------------------------------

Nmax = get(ud.Hframe,'max');

sets = { 'on'  'off' };

%--------------------------------------
% Insert-Button

   nr  = get( ud.edit{1}(1) , 'userdata' );
   nr  = nr{3};

   off = isempty(nr);
   if ~off 
       off = ( nr > Nmax );
   end

   set( h(2) , 'enable' , sets{ 1 + off  } );

%--------------------------------------
% Edit-, Del-Button

   off = isempty( get(ud.Hlist(1),'value') );

   set( h([1 3]) , 'enable' , sets{ 1 + off  } );

%--------------------------------------
% Reset-Button

   set( h(4) , 'enable' , 'on' );



%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [val,str,msg1,msg2]=val2str(val,form)

% VAL2STR  Converts between Valus and Formated Strings
%
% [ Value , String , Message1 , Message2 ] = VAL2STR( Argument , Format )
%
% Argument could be a Value or a String, Format gives the Format for the
%  String.
%
%  Valid Formats are:
%                                                   '#'   String
%  '%##.##'   for using by SPRINTF
%
%  'geo#'     for Geographic Coordinates  Longitude: 1     ##°##.##' <E/W>
%                                         Latitude : 2     ##°##.##' <N/S>
%
%  'time#'    for TimeConversation, first Value Day: 1    DD HH:MM
%                                              Hour: 2       HH:MM
%                                              Hour: 3       HH:MM:SS
%                                            Minute: 4          MM:SS
%
% the Outputs are the String in the given Format 
%  and the correspondending Value 
%    (this Value may be different (rounded)  from the InputValue,
%      depending on the Format )
%
% require functions STR2DAY, DAY2DHMS
%


str  = '';
msg1 = '';
msg2 = '';


if isempty(val) 
 val = [];
 str = '';
 return
end

if isempty(form)
  if ischar(val)
   str = val;
   val = [];
  end
  return
end


if ischar(val)
 if all(val==' ')
  str = val;
  val = [];
  return
 end
else
 if ~isnumeric(val)
  msg1 = ' VAL2STR: Inputs must be Numeric or String (Char).';
  val = [];
  str = '';
  return
 end
end




form = lower(form(:)');

ini = sum( cumsum( strcmp( form(end) , {'1' '2' '3' '4'} ) ) == 0 ) + 1;

ini = ini - 4 * ( ini == 5 );

is_geo  = ( ~isempty( strmatch('geo' ,form) )           );
is_time = ( ~isempty( strmatch('time',form) ) & ~is_geo );



ww = warnstat;

warning('off');


if ischar(val)

 %*******************************************************
 % STRING --> VALUE
 %-------------------------------------------------------

 

 if     is_geo

  ws =  ( ( any(lower(val)=='w')  &  ( ini == 1 ) ) | ...
          ( any(lower(val)=='s')  &  ( ini == 2 ) )       );

  val=str2day(val,2)*24*(1-2*ws);

 elseif is_time

  fak = [ 1  24  24  24*60 ];
  Ini = [ 1   2   2   3    ];
  val = str2day(val,Ini(ini)) * fak(ini);

 else
 
  ok = 1;
  eval( [ 'val = ['  val  '];' ] , 'ok=0;' )

  if ~ok
   msg1 = [ ' VAL2STR : ' lasterr ];
   val = [];
  end

 end

end
% isstr(val) STRING --> VALUE



if isempty(msg1)

%*************************************************
% VALUE --> STRING --> VALUE  using given FORMAT
%-------------------------------------------------

 if ~isempty( findstr(form,'%') )

  str = sprintf(form,val);

  ok = 1;

  eval( [ 'val = ['  str  '];' ] , 'ok=0;' )

  if ~ok
   msg2 = [' VAL2STR: ' lasterr ];
  end

 
 elseif is_geo  |  is_time

    fak = [ 1+23*is_geo   24   24  24*60 ];
   
    [DD,hh,mm,ss] = day2dhms(val/fak(ini) );

   MM = round(mm+ss/60);
   HH = hh + MM/60*( abs(MM) == 60 );
   MM = MM*(~( abs(MM) == 60 ));

   if is_geo

        ws = [ val < 0 ];
        WS = char( (69+18*ws)*(ini==1) + (78+5*ws)*(ini==2) );

        hh = hh+24*DD;

        gf = sprintf('%.0f',4);

       str = sprintf(['%'  gf  '.0f' char(176) ' %5.2f'' ' WS ],abs([hh mm+ss/60]));  

       val = hh + (mm+ss/60)/60;

   elseif is_time

      switch ini
       case 1
         str = sprintf('%3.0f  %2.2d:%2.2d',[DD HH MM]);  
         val = DD + HH/24 + MM/24/60 ;
       case 2  
         str = sprintf('%2.2d:%2.2d',[HH+DD*24 MM]);
         val = HH+DD*24 + MM/60 ; 
       case 3
          ss = round(ss);
          mm = mm + ss/60*[ abs(ss) == 60 ];
          ss = ss*(~( abs(ss) == 60 ));

         str = sprintf('%2.2d:%2.2d:%2.2d',[hh+DD*24  mm ss]);
         val = hh*60+DD*24*60 + mm + ss/60 ;
       case 4
          ss = round(ss);
          mm = mm + ss/60*[ abs(ss) == 60 ];
          ss = ss*(~( abs(ss) == 60 ));

         str = sprintf('%2.2d:%2.2d',[hh*60+DD*24*60+mm ss]);
         val = hh*60+DD*24*60 + mm + ss/60 ;
      end
      % ini
   end
   % form
  end
  % '%'

end


warning(ww);



%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function time=str2day(ss,ini)

% STR2DAY  converts string with any TimeNotation into decimal days
%
% dd = str2day( TimeString , INI ) 
%
%   takes the First Number in the TimeString as
%
%    Day  when INI == 1,
%    Hour      INI == 2,
%    Min       INI == 3, 
%    Sec       INI == 4.
%
% If there is any Error during evaluating the string, 
%  or the string contains an empty value, STR2DAY returns NaN.
%


if nargin < 2
 ini = 1;
end

if isempty(ss) | ~ischar(ss)
 ss = ' ';
end


ss = [ ' '  ss(1,:)  ' ' ];


%------------------------------------------
% Accept only [ 43 45 46 48 .. 57 ]
% "+-.0123456789"
 
ii = find( ss < 43 | ss > 57 | ss == 47 | ss == 44);

ss(ii) = 32+0*ii;


%------------------------------------------
% Accept "." only if Number before or after


ind = [ 2 : length(ss)-1 ];

ii = find( [ss(ind  ) == 46 ]  & ... 
          ~[ss(ind-1) >= 48 ]  & ...
          ~[ss(ind+1) >= 48 ]        );

ss(ii+1) = 32+0*ii;


%------------------------------------------
% Accept "+-" only if any Number follows

jj = find( [ ss == 43 ]  |  [ ss == 45 ] );

ii = find( jj > max(find(ss>=48)) );

ss(jj(ii)) = 32+0*jj(ii);




val = eval( [ '[' char(ss) ']' ] , '[]' );


if isempty(val)

 time = nan;

else

  time = [0 0 0 0];

  ende = ini-1+length(val);
  ende = ende+(4-ende)*[ende>4];

  time( ini : ende ) = val((ini:ende)-ini+1);

if 0
  fak = -1+2*[ time(ini) >=0 ]; 
  if ini < ende
   time(ini+1:ende) = time(ini+1:ende)*fak;
  end
end

  time = sum( time .* [ 1  1/24  1/24/60  1/24/3600 ] );

end


%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [dd,hh,mm,ss] = day2dhms(day)

% DAY2DHMS   converts  dezimal day  into  Day Hour Mine Sec 
%
%  [Day,Hour,Min,Sec] = day2dhms(day)
%

dd = fix( day );
hh = fix( 24 * (day-fix(day)) );
mm = fix( 60 * ( 24 * (day-dd) - hh ));
ss  = round( 24*3600*((day-dd)-hh/24-mm/(24*60)) );

mm_ss = fix(ss/60);
   ss = ss - mm_ss * 60;
   mm = mm + mm_ss;

hh_mm = fix(mm/60);
   mm = mm - hh_mm * 60;
   hh = hh + hh_mm;

dd_hh = fix(hh/24);
   hh = hh - dd_hh * 24;
   dd = dd + dd_hh;


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



