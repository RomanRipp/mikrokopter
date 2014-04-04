function  [Msg,varargout] = sel_list(h,action,varargin);

% SEL_LIST   Creates a SelectionList
%
% [ Msg, ... ] = SEL_LIST( Handle , Action , ... );
%
%  Note:  The first Output contains ErrorMessages for invalid Inputs.
%         If all is ok, Msg is empty.
%
%-------------------------------------------------------------------------
% 1. Create a new SelectList in a Figure, Action 'NEW'
%
%   
% [ Msg, FrameHandle, TextHandle, PixelPosition ] = ...
%    SEL_LIST( FigureHandle, 'New', Property1, PropertyValue1, ... );
%
% Creates a new SelectList in the Figure, defined by 
%   FigureHandle, and returns the Handle of the FrameBox, TitleText, and 
%   the Position in Pixels. The Units for the Position of the ListBox 
%   in the Figure are set to 'pixels'.
% 
% Additional Inputs are Properties and their Values:
%
%  'RowNumber'   , Number of Rows, determines Hight of List
%  'ColNumber'   , Number of Columns, determines Width of List
%  'Marker'      , { Marker  ReMarker }   Strings to Mark, ReMark Selection
%  'Position'    , [  Left  Bottom ]      PixelOffset
%                  [ -Right -Top   ]
%  'BorderWidth' , PixelWidth of UIControl-Border , default: 3
%  'SliderWidth' , SliderWidth in Pixel
%  'ListFont'    ,  ListFontName
%
%  'VisOff'      , { 'Edit' 'Mark' 'UnMark' } to hide this Objects
%
%  'CBFcn'       , CallBackFcn, called after UserActions in the SelList-Object,
%                  Syntax: feval( CBFcn{:} , [ 'Sel' Action ] , 1 , varargin{:} )
%
%  The Actions and varargin are:
%
%   - The selection of ListBoxEntries has changed by a left MouseClick in the ListBox:
%
%     'LIST'  { Mouse SelectionType  ListBoxValue  String(Value)  Mark(Value) }
%               
%   - The Marks of selected ListBoxEntries has changed,
%       by a left MouseClick on Mark- or UnMark-Button:
%
%   { 'MARK'  'UNMARK'  'REMARK' }  { Mark  String  ListBoxString  ListBoxValue }
%
%   - The selection of ListBoxEntries has changed by using the EditBox:
%
%     'EDIT'  { ListBoxValue  EditString }
%  
%
% More Inputs are FontProperties and their Values, by default 
%  the DefaultUIControlFontProperties of the Root will used.
% 
%
% Please:  don't change the Property 'Tag' of the TagFrameHandle,
%           this Property will used to identify the Handle as valid 
%           SelectList.
%
%
% The UserData of the FrameHandle contains all Properties in a StructArray.
%
%-----------------------------------------------------------------------------
% 2. Set a Title above the List, Action 'TITLE'
%
%  SEL_LIST( FrameHandle , 'Title' , Titel );
%  
%-----------------------------------------------------------------------------
% 3. Set Strings to the List, Action 'SET'
%
%  SEL_LIST( FrameHandle , 'Set' , String , Mark , Value )
%
%  String: CellString or CharacterArray
%  Mark:   Vector with same Length of String,
%                 Marked ==  1
%               UnMarked ==  0
%               ReMarked == -1
%  Value: Index to replace Strings in List, if this Input is missing, 
%          all old Strings will removed from List.
%
%-----------------------------------------------------------------------------
% 4. Set the List to specified Values, Action 'VALUE'
%
%  SEL_LIST( FrameHandle , 'Value' , String , Value )
% 
%-----------------------------------------------------------------------------
% 5. Mark, UnMark, Remark ListEntries, Action 'MARK'  'UNMARK'  'REMARK' 
%
%   SEL_LIST( FrameHandle ,   'Mark' , Value )
%                           'UnMark'
%                           'ReMark'
%
% Marks the ListEntries, specified by the Vector Value, if this input is missing,
%  the actual Values from the List will used.
%
%-----------------------------------------------------------------------------
% 6. Returns the State, Action 'GET'
%
% [ Msg , Mark , String , ListString , ListValue ] = SEL_LIST( FrameHandle , 'Get' );
%
%
%-----------------------------------------------------------------------------
% 7. Reset the List, Action 'RESET'
%
%  SEL_LIST( FrameHandle , 'Reset'  )
% 
%-----------------------------------------------------------------------------
% 8. Set the Enability of Mark-, UnMark-Button, EditField
%
%  SEL_LIST( FrameHandle , 'Enable' , sets )
%
%   where sets is 'on' or 'off', following Input can be a CellStringArray with
%     any of: 'Mark'  'UnMark'  'Edit'  
%     to set the Enability only for the specified Objects
%
%-----------------------------------------------------------------------------
% 9. Set the Visibility of List or Mark-, UnMark-Button, EditField
%
%  SEL_LIST( FrameHandle , 'Visible' , sets )
%
%   where sets is 'on' or 'off', following Input can be a CellStringArray with
%     any of: 'List'  'Mark'  'UnMark'  'Edit'  
%     to set the Visibility only for the specified Objects
%
%-----------------------------------------------------------------------------
% 10. Resize the SelectList after the Figure was resized, Action 'RESIZE' 
%
%  [ Msg, PixelPosition ] = SEL_LIST( FrameHandle, 'Resize', Mode )
%
%  Resize the SelectList in the Figure,  depending on Mode
%   call this from the ResizeFcnProperty of the Figure.  
%
%  Mode:  'normalized'  ( short:  'n' )
%
%   The SelectListPosition will set normalized to the Figure, with the normalized
%   Position when the SelectList was created. 
%
%  Mode:  'absolut'     ( short:  'a' )
%
%   The SelectListPosition will set absolut, 
%    defined by the original Position.
%
%-----------------------------------------------------------------------------
% 11. Delete the SelectList, Action 'DELETE'
%
%  MSG = SEL_LIST( FrameHandle , 'Delete' )
%
%



Nout = nargout - 1;

Nout = Nout * ( Nout > 0 );

varargout = cell(1,Nout);



Msg = '';

nl = char(10);


Msg0 = 'SEL_LIST: ';

Nin = nargin;

if Nin < 2
  Msg = [ Msg0  'Inputs H and Action are undefined.' ];
  return
end

%------------------------------------------------------
ok = ( isnumeric(h) & ( prod(size(h)) == 1 ) );
if ok
 ok = ishandle(h);
end

if ~ok
 Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
         ' First Input must be a Handle.' ];
end


%-----------------------------------------------------------
if ~ischar(action)
 Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
          'Action must be a String,'  ];
end

if ~isempty(Msg)
  Msg = [ Msg0   Msg ];
  return
end


%------------------------------------------------------------
action = upper(action);

Msg0 = [ 'SEL_LIST( ' action ' ): ' ];

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
          strcmp(get(h,'tag'),'SEL_FRAME')          );
 end

 if ~ok
   Msg = [ Msg0  'Handle must be a ' ...
           'valid Message-TagFrameUIControl.' ];
   return
 end

 fig = get(h,'parent');

  ud = get(h,'userdata');


end
%------------------------------------------------------------


  VarIn = varargin;
  VarIn = VarIn(:);
 
  Nin   = size(VarIn,1);

SelCall = ( Nin == 0 );  % SEL_LIST-Callback

CallIn  = cell(1,0);     % Input for CBFcn;

switch action

%*********************************************************
case 'NEW'

  
  if ( mod(Nin,2) ~= 0 )
    Msg = [ 'Additional Inputs must contain UIControl-Property-Value-Pairs.' ];
  else 
    Nin = Nin / 2;
    VarIn = reshape( VarIn , 2 , Nin )';
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

  scr_uni = get(0,'units');      set(0,'units','pixels')
  scr_si  = get(0,'ScreenSize'); set(0,'units',scr_uni);
  
  is_win = strcmp( upper(computer) , 'PCWIN' );

  if is_win
     ListFont = 'courier';
  else
     ListFont = { 'arrial'  get(0,'fixedwidthfontname') };
     ListFont = ListFont{ 1 + ( scr_si(4) >= 1050 ) } ;
  end 

  visbl = struct( 'Main' , { 'on' } , ...
                  'Edit' , { 'on' } , ...
                  'Mark' , { 'on' } , ...
                'UnMark' , { 'on' }       );


  enabl =  struct( 'Main' , { 'on' } , ...
                   'List' , { 'on' } , ...
                   'Edit' , { 'on' } , ...
                   'Mark' , { 'on' } , ...
                 'UnMark' , { 'on' }       );


  Marker = { ' x '  
             '   '   
             ' ~ '  };  % { 1  0  -1 }

  ud = struct( 'RowNumber'     , { 12  }         , ... % Nr of Row's
               'ColNumber'     , { 16  }         , ... % Nr of Row's
               'Marker'        , { Marker }      , ... % String to Mark Selection
               'Position'      , {[ 10 10 ]}     , ... % Position 
               'BorderWidth'   , {  3 }          , ... % BorderPixelWidth of UIControl
               'SliderWidth'   , { 15+3*is_win } , ... % SliderWidth in Pixel
               'ListFont'      , { ListFont }    , ... % FontName for List
               'CBFcn'         , { {} }          , ... % CallBackFcn
               'CharFit'       , { zeros(2,2,2) }  , ... % [ Width High ] of [ ListFont ; Font ]
               'LineOffset'    , {  2-[0 2]*is_win  } , ... % ListBox, [ Width hight ]
               'FrameHandle'   , { NaN }         , ... % Handle of Frame arround
               'TextHandle'    , { NaN }         , ... % Handle of TextObject, using for width
               'MarkHandle'    , { NaN }         , ... % Handle of Buttons
               'UnMarkHandle'  , { NaN }         , ... % Handle of Buttons
               'ListHandle'    , { NaN }         , ... % Handle of TextObject, using for width
               'EditHandle'    , { NaN }         , ... % Handle of TextObject, using for widt
               'Visibility'    , { visbl      }  , ...
                'Enability'    , { enabl      }  , ...
               'ResizeMode'    , { 'absolut'  }  , ... % Actual ResizeMode
               'PixelPosition' , { zeros(1,4) }  , ... % Actual PixelPosition [ Left Bottom Width High ]
               'NormPosition'  , { zeros(1,4) }  , ... % Default normalized Position
               'FigPos0'       , { figpos     }  );    % Default FigurePosition


  %-------------------------------------------------------------------
  % Look, if the first 8 Properties are given in Inputs

  fields = fieldnames(ud);
  fields = fields([1 2 3 4 5 6 7 8]);


  is_ud = zeros(Nin,1);

  for ii = 1 : Nin

     jj = find( strcmp( lower(VarIn{ii,1}) , lower(fields) ) );

     is_ud(ii) = ~isempty(jj);

     if is_ud(ii)

       val = VarIn{ii,2};
       msg = '';

       switch fields{jj}

         %-----------------------------------------------
         case  { 'ColNumber' 'RowNumber'  'Position'  'BorderWidth' }

            if ~isnumeric(val) | isempty(val)

               msg = [ 'Value for '  fields{jj} ' must be numeric and not empty.' ];

            else

            val_min = 1;

              if strcmp(fields{jj},'Position')
                 val_min = -inf;
              end

              ok1 =  ( all( val >= val_min )  &  all(isfinite(val)) );
              ok2 =  ( all( mod(val,1) == 0 )  |  strcmp(fields{jj},'Position') );
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


         %-----------------------------------------------
         case  'Marker'

            if ischar(val)
               val = cellstr(val);
            end

            ok = ( iscellstr(val) & ~isempty(val) );
            if ok

               val = val(:);

               sv  = size(val,1);
               sm  = size(ud.Marker,1);
                n  = min(sm,sv);

               for ii = 1 : n
                   ok1 = ( ischar(val{ii}) & ~isempty(val{ii}) & ...
                           ( prod(size(val{ii})) == size(val{ii},2) ) );
                   ok = ( ok & ok1 );
               end

            end

            if ok

               val = cat( 2 , val(1:n) , ud.Marker(n+1:sm) ); 

            else

               msg = [  msg  nl(1:(end*(~isempty(msg)))) ...
                        'Value for '  fields{jj} ' must contains Strings.' ];

            end
   
         %-----------------------------------------------
         case  'Listfont' 

            ok = ( ischar(val) & ~isempty(val) & ...
                   ( prod(size(val)) == size(val,2) ) );
            if ~ok
               msg = [  msg  nl(1:(end*(~isempty(msg)))) ...
                        'Value for '  fields{jj} ' must contains a String.' ];
            end

         %-----------------------------------------------
         case  'CBFcn'

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

       end

        if isempty(msg)
             ud     = setfield(ud,fields{jj},val);
        else
             Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) msg ];
        end

     end
     % is_ud(ii)

  end
  % ii

  if ~isempty(Msg)
    Msg = [ Msg0  Msg ];
    return
  end

  
  VarIn(find(is_ud),:) = [];


  %-------------------------------------------------------------------
  % Special Search for 'VisOff'

  is_vis = any( strcmp( lower(VarIn(:,1)) , 'visoff' ) );
 
  if is_vis

     jj = find( strcmp( lower(VarIn(:,1)) , 'visoff' ) );

     typ = VarIn{jj(end),2};

     VarIn(jj,:) = [];

     if ischar(typ)
        typ = cellstr(typ);
     end

     typs = { 'List' 'Edit' 'Mark' 'UnMark' };
     typs = typs( 1+strcmp(action,'VISIBLE') : end );

     ok = iscellstr(typ);
     if ok

        typ  = typ(:);
        
        for ii = 1 : size( typ , 1 )

            ok1 = ( ischar(typ{ii}) & ~isempty(typ{ii}) & ...
                    ( prod(size(typ{ii})) == size(typ{ii},2) ) );
            if ok1
               ok1 = any( strcmp( lower(typ{ii}) , lower(typs) ) );
            end

            ok = ( ok & ok1 );

        end
     end

     if ~ok

       Msg = [ Msg0  'Value for VisOff must contains any of '  ...
                     strhcat(typs,', ') '.' ];
       return
        
     end
  
  end

  %-------------------------------------------------------------------


  Nin = size(VarIn,1);


  %-------------------------------------------------------------------
  % Check UIControlProperties in the other Inputs

  VarIn(:,3) = { [] };  % DefaultValues

  
  % Set DefaultProperties 

  for ii = 1 : Nin

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

   FontName = { ud.ListFont  get(ht,'FontName') };

   nf = size(FontName,2);

   ext = zeros(2,4,nf);

   str = 'H';
  
   for ii = 1 : 2

     set( ht , 'string' , str(ones(1,ii),ones(1,ii)) );

     for jj = 1 : nf

        set( ht , 'FontName' , FontName{jj} );

        ext(ii,:,jj) = ceil( get(ht,'extent') );

     end

   end

   delete(ht);
   delete(axe);

   m = ext(2,[3 4],:) - ext(1,[3 4],:);
   n = ext(1,[3 4],:) - m;

  
   ud.CharFit = cat(1,m,n);



   %--------------------------------------------------------
   % Get ButtonImage

    bc0 = [ 0.95 1.0 1.0 ];

     cc = cell(1,2);

     ff = { 'mark.xpm'  'unmark.xpm' };

     for ii = [ 1  2 ]

        try

          [Msg,cc{ii}] = readxpm( ff{ii} );

          if isempty(Msg)
         
              jj = find(isnan(sum(cc{ii},3)));
 
            si12 = size(cc{ii},1) * size(cc{ii},2);

            cc{ii}(jj+0*si12) = bc0(1);
            cc{ii}(jj+1*si12) = bc0(2);
            cc{ii}(jj+2*si12) = bc0(3);

          end

        end
        
     end

     if isempty(cc{1})  |  isempty(cc{2})

         cc = cell(1,2);

        str = { 'x'  'o' };
         bc = [ 0.75  1.00  0.75
                1.00  0.75  0.75  ];

     else

         str = { ''  '' };
         bc  = bc0( [ 1  1 ] , : );

     end


   %--------------------------------------------------------
   % Build UIControls

   % Frame arround

   ud.FrameHandle = uicontrol('parent'   , fig          , ...
                              'units'    , 'pixels'     , ...
                              'position' ,  [ 0  0  1  1 ] , ...
                              'style'    , 'frame'      , ...
                              'selected' , 'off'        , ...
                              'enable'   , 'on'         , ...
                              'visible'  , 'on'         , ...
                              'tag'      , 'SEL_FRAME'  );

   form = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 );

     HF = sprintf(form,ud.FrameHandle);


   DelCB = [ 'sel_list('  HF  ',''Delete'',-1);' ];

   set(ud.FrameHandle,'DeleteFcn',DelCB);

   %--------------------------------------------------------

   ud.TextHandle = uicontrol('parent'   , fig     , ...
                             'units'    , 'pixels', ...
                             'position' ,  [ 0  0  1  1 ] , ...
                             'style'    , 'text'       , ...
                             'string'   , ''           , ...
                             'min'      , 0            , ...
                             'max'      , 1            , ...
                             'selected' , 'off'        , ...
                             'enable'   , 'on'         , ...
                             'visible'  , 'on'         , ...
                  'horizontalalignment' , 'center'     , ...
                             'tag'      , 'SEL_TEXT'  );


   %--------------------------------------------------------
   % EditHandle below
     
   CB = [ 'sel_list('  HF  ',''Edit'');'  ];

   ToolTip = [ 'Enter Index-Vector to select,'                    nl ...
               'Begin/End with ":" to start/end with First/Last.' nl ...
               'Use "end" for Last, ":" or "all" to select All.'          ];

   ud.EditHandle = uicontrol('parent'   , fig     , ...
                             'units'    , 'pixels', ...
                             'position' ,  [ 0  0  1  1 ] , ...
                             'style'    , 'edit'       , ...
                             'string'   , ''           , ...
                             'min'      , 0            , ...
                             'max'      , 1            , ...
                             'callback' , CB           , ...
                             'selected' , 'off'        , ...
                             'enable'   , 'off'        , ...
                             'visible'  , visbl.Edit   , ...
                        'ToolTipString' , ToolTip      , ...
                      'BackGroundColor' , bc0          , ...
                  'horizontalalignment' , 'center'     , ...
                             'tag'      , 'SEL_EDIT'  );


   %--------------------------------------------------------
   % MarkButton on Top

   CB = [ 'sel_list('  HF  ',''Mark'');'  ];

   ToolTip = 'Mark Selection';

   ud.MarkHandle = uicontrol('parent'   , fig     , ...
                             'units'    , 'pixels', ...
                             'position' ,  [ 0  0  1  1 ] , ...
                             'style'    , 'pushbutton' , ...
                             'string'   , str{1}       , ...
                             'cdata'    ,  cc{1}       , ...
                             'callback' , CB           , ...
                             'selected' , 'off'        , ...
                             'enable'   , 'off'        , ...
                             'visible'  , visbl.Mark   , ...
                        'ToolTipString' , ToolTip      , ...
                      'BackGroundColor' ,  bc(1,:)     , ...
                  'horizontalalignment' , 'center'     , ...
                             'tag'      , 'SEL_MARK'  );


   %--------------------------------------------------------
   % UnMarkButton on Top

   CB = [ 'sel_list('  HF  ',''UnMark'');'  ];

   ToolTip = 'UnMark Selection';

   ud.UnMarkHandle = uicontrol('parent'   , fig     , ...
                               'units'    , 'pixels', ...
                               'position' ,  [ 0  0  1  1 ] , ...
                               'style'    , 'pushbutton' , ...
                               'string'   , str{2}       , ...
                               'cdata'    ,  cc{2}       , ...
                               'callback' , CB           , ...
                               'selected' , 'off'        , ...
                               'enable'   , 'off'        , ...
                               'visible'  , visbl.UnMark , ...
                          'ToolTipString' , ToolTip      , ...
                        'BackGroundColor' ,  bc(2,:)     , ...
                    'horizontalalignment' , 'center'     , ...
                               'tag'      , 'SEL_UNMARK'  );
                      

   %--------------------------------------------------------
   % ListBox

   udl = struct( 'String' , {  cell(0,1) } , ...
                 'Marked' , { zeros(0,1) }       );

   CB = [ 'sel_list('  HF  ',''LIST'');'  ];

   ud.ListHandle = uicontrol('parent'   , fig     , ...
                             'units'    , 'pixels', ...
                             'position' ,  [ 0  0  1  1 ] , ...
                             'style'    , 'listbox'    , ...
                             'FontName' , ud.ListFont  , ...
                             'value'    ,  []          , ...
                             'string'   , {''}         , ...
                             'min'      , 0            , ...
                             'max'      , 2            , ...
                             'callback' , CB           , ...
                             'selected' , 'off'        , ...
                             'enable'   , 'off'        , ...
                             'visible'  , 'on'         , ...
                            'userdata'  , udl          , ...
                  'horizontalalignment' , 'left'       , ...
                             'tag'      , 'SEL_LIST'  );


   %--------------------------------------------------------

   % Set DefaultProperties back
   for ii = 1 : Nin
     set(0,['DefaultUIControl' VarIn{ii,1}],VarIn{ii,3});
   end


   set( ud.FrameHandle , 'userdata' , ud );

   %-------------------------------------------
   % Set correct Size of UIControls

   try

     if is_vis
          
       [Msg,pos] = sel_list(ud.FrameHandle,'Visible','off',typ);
       
     else

       [Msg,pos] = sel_list(ud.FrameHandle,'Resize','absolut');

     end

   catch
      Msg = lasterr;
   end

   if ~isempty(Msg)

     Msg = [ Msg0 'Error call SEL_LIST( Resize ).' nl Msg ];

     set(ud.FrameHandle,'DeleteFcn','');

     delete(ud.FrameHandle);
     delete(ud.TextHandle);
     delete(ud.MarkHandle);
     delete(ud.UnMarkHandle);
     delete(ud.ListHandle);
     delete(ud.EditHandle);

     return

   end

   ud = get( ud.FrameHandle , 'userdata' );
 
   ud.PixelPosition = pos;
   ud.NormPosition  = pos ./ figpos([ 3  4  3  4 ]);

   set( ud.FrameHandle , 'userdata' , ud );


   %-------------------------------------------

   out = { ud.FrameHandle ud.TextHandle pos }; 

   n = min(Nout,size(out,2));

   varargout(1:n) = out(1:n);


%*********************************************************
case 'RESIZE'

  if SelCall

    mode = 'a';   % Normalized

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


  % Check Visibility

  is_ed = strcmp( ud.Visibility.Edit   , 'on' );
  is_mk = strcmp( ud.Visibility.Mark   , 'on' );
  is_um = strcmp( ud.Visibility.UnMark , 'on' );
    
  % CharacterHight

  ch = ud.CharFit(1,2,:);
 
  sc = cat( 1 , size( get(  ud.MarkHandle,'cdata') ) , ...
                size( get(  ud.UnMarkHandle,'cdata') )       );

  sc = max(sc(:));

  sc = max( sc , ch(2) );

  wwc = sc + 2*ud.BorderWidth;  % ButtonWidth
  hhc = sc + 2*ud.BorderWidth;  % ButtonHight

  hht = ch(2);                     % TextHight

  hhe = ch(2) + 2*ud.BorderWidth;   % EditHight

  %--------------------------------------------------------------------------
  if strcmp(mode,'n')
  % Normalized

     pos = ud.NormPosition .* figpos([3 4 3 4]);

     pos([1 3]) =  ceil( pos([1 3]) );
     pos([2 4]) = floor( pos([2 4]) );
 
     wwl = pos(3) - 2*ud.BorderWidth;
     hhl = pos(4) - 2*ud.BorderWidt - is_ed*hhe - hhc;

     wwl = max( wwl , 1 );
     hhl = max( hhl , 1 );
    
  %--------------------------------------------------------------------------
  else
  % Absolut

    %  'Position'  , [  Left  Bottom  ] , default: [ 10  10 ]
    %                [  Left  Bottom  ]
    %                [ -Right -Top    ]

   hhl = ceil( ud.RowNumber * ( ch(1) + ud.LineOffset(2) ) + ...
               2*ud.BorderWidth );    % ListBoxHight    

   wwl =  ud.ColNumber * ud.CharFit(1,1,1) + ud.CharFit(2,1,1); 

   wwl = ceil( wwl +  ud.LineOffset(1) + ...
               ud.SliderWidth + 2*ud.BorderWidth );  % ListBoxWidth

     pos = zeros(1,4);

     pos(3) = wwl + 2*ud.BorderWidth;
     pos(4) = is_ed*hhe + hhl + hhc + 2*ud.BorderWidth;

     po = ud.Position;

     pos([1 2]) = ( 1 + po ) .* ( po >= 0 );

     pos([1 2]) = pos([1 2]) + ...
                   ( figpos([3 4]) - pos([3 4]) + 1 + po([1 2]) ) .* ...
                   ( po([1 2]) < 0 );

     pos([1 3]) =  ceil( pos([1 3]) );
     pos([2 4]) = floor( pos([2 4]) );
 
             
  end
  %--------------------------------------------------------------------------

  wwt = pos(3) - 2*ud.BorderWidth - (is_mk+is_um)*wwc;
  wwt = max( wwt , 1 );

  pos(3) = max( pos(3) , wwt + (is_mk+is_um)*wwc + 2*ud.BorderWidth );
  pos(4) = max( pos(4) , is_ed*hhe + hhl + hhc + 2*ud.BorderWidth );

  ud.PixelPosition = pos;
  ud.ResizeMode    = mode;

  %---------------------------------------------------------
  % Visibility

  sets = { 'off'  'on' };

  is_main = strcmp( ud.Visibility.Main , 'on' );

  %---------------------------------------------------------

  set( ud.FrameHandle , 'units'     , 'pixels' , ...
                        'position'  ,  pos     , ...
                        'visible'   , sets{1+is_main}, ...
                        'userdata'  ,  ud             )



  %---------------------------------------------------------
  % Edit

  pos1 =  cat( 2 , pos([1 2])+ud.BorderWidth , wwl , hhe );

  is_vis = ( is_main  &  strcmp( ud.Visibility.Edit , 'on' ) );

  set( ud.EditHandle , 'units'     , 'pixels' , ...
                       'position'  ,  pos1    , ...
                       'visible'   , sets{1+is_vis}       )

  %---------------------------------------------------------
  % List

  pos1 = cat( 2 , pos(1)+ud.BorderWidth , pos1(2)+pos1(4)*is_ed , wwl , hhl );
 
  set( ud.ListHandle , 'units'     , 'pixels' , ...
                       'position'  ,  pos1    , ...
                       'visible'   , sets{1+is_main}  )


  %---------------------------------------------------------
  % Text

  pos1 = cat( 2 , pos(1)+ud.BorderWidth , ...
                 pos1(2)+pos1(4)+ud.BorderWidth , wwt , hht );
 
  set( ud.TextHandle , 'units'     , 'pixels' , ...
                       'position'  ,  pos1    , ...
                       'visible'   , sets{1+is_main}       )


  %---------------------------------------------------------
  % Buttons

  pos2 = [ 0  0  wwc  hhc ];

  pos2([1 2]) = pos([1 2]) + pos([3 4]) - ud.BorderWidth - pos2([3 4]);
 

  is_vis = ( is_main  &  strcmp( ud.Visibility.Mark , 'on' ) );

  set( ud.MarkHandle , 'units'     , 'pixels' , ...
                       'position'  ,  pos2-[is_um*wwc 0 0 0] , ...
                       'visible'   , sets{1+is_vis}    )


  is_vis = ( is_main  &  strcmp( ud.Visibility.UnMark , 'on' ) );

  set( ud.UnMarkHandle , 'units'     , 'pixels' , ...
                         'position'  ,  pos2    , ...
                         'visible'   , sets{1+is_vis}  )

  
   %-------------------------------------------

   out = { ud.PixelPosition };  % !!!

   n = min(Nout,size(out,2));

   varargout(1:n) = out(1:n);


%***********************************************************************
case 'LIST'

   udl = get( ud.ListHandle , 'userdata' );

   if isempty( udl.String )
      return
   end

   val = get( ud.ListHandle , 'value' );
   sel = get( fig , 'SelectionType' );

   %-------------------------------------------

   out = { val  udl.String(val)  udl.Marked(val)  sel }; 

   n = min(Nout,size(out,2));

   varargout(1:n) = out(1:n);


   CallIn = { sel  val  udl.String(val)  udl.Marked(val) };


%***********************************************************************
case 'SET'

  if SelCall
    return
  end

  %---------------------------------------------------------------------
  str = VarIn{1};

  if isempty(str)

     sel_list( ud.FrameHandle , 'Reset' );

     return

  end

  if ischar(str)
    
    str = cellstr(str);

  end

  str = str(:);

  if ~iscellstr(str)

     Msg = [ Msg0 'First Input must be a Character- or CellStringArray.'];

     return

  end

  %---------------------------------------------------------------------
  if Nin < 2

    val = zeros(size(str,1),1);

  else

    val = VarIn{2};
    ok = ( isnumeric(val)  &  ~isempty(val) );
    if ok
       val = val(:);
        ok =  ( ( size(val,1) == size(str,1) )  &  ...
                 all(isfinite(val))  );
        if ok
           val = sign(val);
        end
    end

    if ~ok

      Msg = [ Msg0 'Second Input must be a Numeric Vector' nl ...
                   ' with same Size as String.'];
 
      return

    end

  end

  %---------------------------------------------------------------------
  if Nin < 3
  
     ind = [];

  else

     ind = VarIn{3}(:);

     Msg = check_ind(ind);

     if ~isempty(Msg)
        Msg = [ Msg0  Msg ];
        return
     end

    if ~( size(ind,1) == size(str,1) )

      Msg = [ Msg0 'Third Input must have same Size as String.'];
 
      return

    end

  end
     
  %---------------------------------------------------------------------


  udl = get( ud.ListHandle , 'userdata' );

  if isempty(ind)

    udl.String = str;
    udl.Marked = val;

  else

    n = size( udl.String , 1 );
    m =  max( ind );

    s = max( n , m );

    str1 =  cell( s , 1 );
    num1 = zeros( s , 1 );
 
    if n > 0

       str1(1:n) = udl.String;
       num1(1:n) = udl.Marked;

    end

    str1(ind) = str;
    num1(ind) = val;

    udl.String = str1;
    udl.Marked = num1;

  end

  set( ud.ListHandle , 'userdata' , udl );

  Msg = sel_list( ud.FrameHandle , 'NewList' );


   %-------------------------------------------

   out = { udl.Marked  udl.String  get(ud.ListHandle,'string') }; 

   n = min(Nout,size(out,2));

   varargout(1:n) = out(1:n);

 
%***********************************************************************
case 'GET'

   udl = get( ud.ListHandle , 'userdata' );
   str = get( ud.ListHandle , 'string' );
   val = get( ud.ListHandle , 'value' );

   %-------------------------------------------

   out = { udl.Marked  udl.String  str  val  }; 

   n = min(Nout,size(out,2));

   varargout(1:n) = out(1:n);


%***********************************************************************
case 'RESET'


   udl = struct( 'String' , {  cell(0,1) } , ...
                 'Marked' , { zeros(0,1) }         );


   set( ud.ListHandle , 'value'    ,  []    , ...
                        'string'   , { '' } , ...
                        'userdata' , udl          );

    hh = [  ud.ListHandle  ud.EditHandle  ...
            ud.MarkHandle  ud.UnMarkHandle ];

    set( hh , 'enable' , 'off' );
    

%***********************************************************************
case 'TITLE'


   if SelCall
      return
   end

   str = VarIn{1};
           
   if iscellstr(str)
      str = char(str);
   end

   if ~( ischar(str) & ...
         ( prod(size(str)) == size(str,2) ) )

       Msg = [ Msg0   'Input must be a String.' ];

       return

   end

   try
      set( ud.TextHandle , 'string' , str , VarIn{2:end} );
   catch
      Msg = [ Msg0  'Invalid Inputs.' nl lasterr ];
   end


%***********************************************************************
case { 'MARK'  'UNMARK'  'REMARK' }

  if SelCall

     val = get( ud.ListHandle , 'value' );

  else

     val = VarIn{1};
  
     [ Msg , val ] = sel_list( ud.FrameHandle , 'CheckValue' , val );

     if ~isempty(Msg)
        Msg = [ Msg0  Msg ];
        return
     end

   end

   udl = get( ud.ListHandle , 'userdata' );

   if isempty( udl.String )
      return
   end

   udl.Marked(val) =  1 * strcmp( upper(action) ,   'MARK' ) + ...
                     -1 * strcmp( upper(action) , 'REMARK' );

   set( ud.ListHandle , 'userdata' , udl );


   sel_list( ud.FrameHandle , 'NewList' );


   str = get(ud.ListHandle,'string');

   %-------------------------------------------

   out = { udl.Marked  udl.String  str }; 

   n = min(Nout,size(out,2));

   varargout(1:n) = out(1:n);


   CallIn = { udl.Marked  udl.String  str  val };



%***********************************************************************
case 'NEWLIST'

   udl = get( ud.ListHandle , 'userdata' );

   if isempty( udl.String )

      sel_list( ud.FrameHandle , 'Reset' );

      return

   end


   n = size( udl.String , 1 );

   str = cell( n , 2 );
   str(:,2) = udl.String;

   mark = char(ud.Marker);

   str(:,1) = { char( 32 * ones(1,size(mark,2)) ) }; 


   for ii = [  1   0  -1 ]
   
      jj  = find( udl.Marked == ii );

      str(jj,1) = { mark(-ii+2,:) };

   end


   for ii = 1 : n

       str{ii,1} = cat( 2 , str{ii,1} , str{ii,2} );
     
   end

   val = get( ud.ListHandle , 'value' );
   val = min( val , n );

   set( ud.ListHandle , 'value' , val , ...
                        'string' , str(:,1) );


   sel_list( ud.FrameHandle , 'SetEnable' );


%***********************************************************************
case 'EDIT'

  if SelCall

    str = get(ud.EditHandle,'string');
 
  else

    str = VarIn{1};

  end

  str = lower(rmblank(str,2));

  if isempty(str)
     return
  end

  udl = get(ud.ListHandle,'userdata');

    n = size(udl.String,1);

  str = strrep(str,'end',sprintf('%.0f',n));


  if any( strcmp( str , { ':'  'all' } ) )

    val = ( 1 : n );
 
  else

    str = strrep(str,'end',sprintf('%.0f',n));

    if strcmp( str(1) , ':' )
       str = cat( 2 , '1' , str );
    end

    if strcmp( str(end) , ':' )
       str = cat( 2 , str , sprintf('%.0f',n) );
    end

    [ Msg , val ] = check_fcn([ '[' str ']' ]);

    if ~isempty(Msg)

       Msg = [ 'Invalid Syntax.' nl Msg ];

       if SelCall
         set( ud.EditHandle , 'string' , '' );
         warndlg(Msg,'Invalid Input.','warn')
       end

       return

    end

  end


  [ Msg , val ] = sel_list( ud.FrameHandle , 'Value' , val );

  if ~isempty(Msg)

     if SelCall
       set( ud.EditHandle , 'string' , '' );
       warndlg(Msg,'Invalid Input.','warn')
     end

     return

  end


   %-------------------------------------------

   out = { val str }; 

   n = min(Nout,size(out,2));

   varargout(1:n) = out(1:n);


   CallIn = { val  str };


%***********************************************************************
case 'VALUE'

  if SelCall
     return
  end

  val = VarIn{1};

  [ Msg , val ] = sel_list( ud.FrameHandle , 'CheckValue' , val );

  if ~isempty(Msg)

    Msg = [ Msg0  Msg ];
    return

  end

  top = max( 1 , get( ud.ListHandle , 'listboxtop' ) );

  if ~isempty(val)

    ns  = size( get(ud.ListHandle,'string') , 1 );

    tmax = ns - ( ud.RowNumber - 1 ) * ( ns > ud.RowNumber );
    vmin = min(val);

    top = top + ( vmin - top ) * ( top+ud.RowNumber-1 < vmin );
    top = min( top , tmax );

   end

   set( ud.ListHandle , 'value'      , val , ...
                        'listboxtop' , top       );


   %-------------------------------------------

   out = { val }; 

   n = min(Nout,size(out,2));

   varargout(1:n) = out(1:n);


%***********************************************************************
case 'CHECKVALUE'

  if SelCall
     return
  end

  val = VarIn{1};

  Msg = check_ind(val);

  if ~isempty(Msg)
     return
  end

  if isempty(val)
     return
  end


  udl = get(ud.ListHandle,'userdata');

  lim = size(udl.String,1);


  val0 = val;

  val = sort(val);
  val = val( 1 : sum(val<=lim) );

  if isempty(val)  &  ~isempty(val0)

     Msg = 'No Index in Range.';

  end

   %-------------------------------------------

   out = { val }; 

   n = min(Nout,size(out,2));

   varargout(1:n) = out(1:n);


%***********************************************************************
case { 'VISIBLE'  'ENABLE' }

   
  ok = ~SelCall;

  if ok

    sets = VarIn{1};

    ok = ( ischar(sets) & ~isempty(sets) & ...
           ( prod(size(sets)) == size(sets,2) ) );
    if ok
       ok = any(strcmp(sets,{'on'  'off'}));
    end

  end

  if ~ok

    Msg = [ Msg0  'Input for must be ''on'' or ''off''.'];
    return

  end

  if strcmp( action , 'VISIBLE' )
     field = 'Visibility';
  else
     field = 'Enability';
  end

  %-------------------------------------------
  if Nin >= 2

     typ = VarIn{2};

     if ischar(typ)
        typ = cellstr(typ);
     end

     typs = { 'List' 'Edit' 'Mark' 'UnMark' };
     typs = typs( 1+strcmp(action,'VISIBLE') : end );

     ok = iscellstr(typ);
     if ok

        typ  = typ(:);
        
        for ii = 1 : size( typ , 1 )

            ok1 = ( ischar(typ{ii}) & ~isempty(typ{ii}) & ...
                    ( prod(size(typ{ii})) == size(typ{ii},2) ) );
            if ok1
               ok1 = any( strcmp( lower(typ{ii}) , lower(typs) ) );
            end

            ok = ( ok & ok1 );

        end
     end

     if ~ok

       Msg = [ Msg0  'Second Input must contain any of '  ...
                     strhcat(typs,', ') '.' ];
       return
        
     end


     val = getfield( ud , field );
     
     for ii = 1 : size( typ , 1 )

         jj = find( strcmp( lower(typ{ii}) , lower(typs) ) );
         
         val = setfield( val , typs{jj} , sets );

     end

     ud = setfield( ud , field , val );

  %-------------------------------------------
  else

    val      = getfield( ud , field );

    val.Main = sets;

    ud = setfield( ud , field , val );
         
  end


  set( ud.FrameHandle , 'userdata' , ud );


  %-------------------------------------------
  if strcmp( action , 'VISIBLE' )

    [Msg , pos ] = sel_list( ud.FrameHandle , 'Resize' , ud.ResizeMode );

     %-------------------------------------------

     out = { pos }; 

     n = min(Nout,size(out,2));

     varargout(1:n) = out(1:n);

  %-------------------------------------------
  else

    sel_list( ud.FrameHandle , 'SetEnable' );
  
  end


%***********************************************************************
case 'SETENABLE'


  sets = { 'off'  'on' };


  is_main = strcmp( ud.Enability.Main , 'on' );

  udl = get( ud.ListHandle , 'userdata' );

  is_list = ~isempty(udl.String);


  field = { 'List' 'Edit' 'Mark' 'UnMark' };

  for ff = field(:)'

      enabl = getfield( ud.Enability , ff{1} );

      is_enabl = ( is_main  &  is_list  &  strcmp( enabl , 'on' ) );

      h = getfield( ud , [ ff{1} 'Handle' ] );

      set( h , 'enable' , sets{ 1 + is_enabl } );

  end


%***********************************************************************
case 'DELETE'

   
    hh = [ ud.TextHandle ; ud.MarkHandle ; ud.UnMarkHandle ; ...
           ud.ListHandle ; ud.EditHandle ];

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


%***********************************************************************

if SelCall & ~isempty(ud.CBFcn) & ~isempty(CallIn)
 
       try
         feval( ud.CBFcn{:} , [ 'Sel' action ] , 1 , CallIn{:} );
       catch
         msg = lasterr;
       end
 
       if ~isempty(Msg)

          msg = [ Msg0 ' Error in UserCallBackString. ' nl ...
                        msg nl ];

          fprintf([ msg  nl ]);

       end

end


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function Msg = check_ind(ind)

% CHECK_IND  checks Index for valid Values ( Integer > 0 )

 Msg = '';

 ok = isnumeric(ind);
 if ok
    if ~isempty(ind)
       ind = ind(:);
        ok = all( ( ind >= 1 )  &  ( mod(ind,1) == 0 ) ...
                    &  isfinite(ind) ); 
    end
 end

 if ~ok
    Msg = ['Index must contains positive Integers.' ];
 end

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,out] = check_fcn(FUNCTION)


 FUNCTION = deblank(FUNCTION);
 if ~strcmp(FUNCTION(end),';');
   FUNCTION = [ FUNCTION  ';'  ];
 end
  
 eval([ 'out = ' FUNCTION ],'Msg = lasterr;');

 if exist('Msg') ~= 1

   Msg = '';

 end

 if exist('out') ~= 1

   out = [];

 end


