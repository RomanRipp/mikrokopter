function  [Msg,varargout] = msg_list(h,action,varargin);

% MSG_LIST   Creates a Message-ListBox and stores Messages into
%
% [ Msg, ... ] = MSG_LIST( Handle , Action , ... );
%
%  Note:  The first Output contains ErrorMessages for invalid Inputs.
%         If all is ok, Msg is empty.
%
%-------------------------------------------------------------------------
% 1. Create a new MessageListBox in a Figure, Action 'NEW'
%
%   
% [ Msg, ListHandle, PixelPosition ] = MSG_LIST( ...
%        FigureHandle, 'New', Property1, PropertyValue1, ... );
%
% Creates a new MessageListBox in the top of the Figure, defined by 
%   FigureHandle, and returns the Handle of the ListBox and 
%   the Position in Pixels. The Units for the Position of the ListBox 
%   in the Figure are set to 'pixels'.
% 
% Additional Inputs are Properties and their Values:
%
%  'RecentNumber'  , [ Number of recent ListEntries ] , default: 50 
%  'RowNumber'     , [ Number of displayed Rows     ] , default: 5
%  'PixelOffset'   , [  Left  Right  Top    ]         , default: [0 0 0]
%                    [  Left  Right -Bottom ]
%  'AppendString'  , [ String  ]                      , default: ' ... '
%  'NumberForm'    , [ FormatString for ListNumber ]  , default: '%4.3d:  '
%
% More Inputs are FontProperties and their Values, by default 
%  the DefaultUIControlFontProperties of the Figure will used.
% 
% Note: The High of the MessageListBox will defined by the FontSize, 
%        FontWeigh and the RowNumber.
%       Negative Values for PixelOffset-Top will used as Offset from Bottom !!!
%
% Please:  don't change the FontProperties of the ListBoxHandle and
%          don't change the Property 'Tag' of the ListBoxHandle,
%           this Property will used to identify the Handle as valid 
%           MessageListBox.
%
% For the using of AppendString and NumberForm see below under
%   Action MESSAGE.
%
% The UserData of the ListBoxHandle contains all Properties in a StructArray.
%
%----------------------------------------------------------------------------
% 2. Store a new Messages in the MessageListBox , Action 'MESSAGE', Mode 'new'
%
% [ Msg, ListNumber, NewString ] = MSG_LIST( ListHandle, 'Message', Text );
%  
% [ Msg, ListNumber, NewString ] = MSG_LIST( ListHandle, 'Message', Text, 'new' );
%
% Creates a new ListEntry from Text, and returns the ListNumber of this Entry,
%   and the new String, added to the ListBox.
%
% The NewString is build like:
%
%   NewString = { [ sprintf(NumberForm,ListNumber)  Text ] };  % CellString !!!
% 
% If Text contains NewLineCharacters ( char(10) ), a MultiLine-String will added.
%  In this case, the Output NewString  contains multiple Rows.
%
% The ListBoxtop is set automaticly, that the new Message is full visible.
% Only the last <RecentNumber> Messages will shown in the ListBox, 
%   that saves Memory.
%
% [ Msg, ListNumber, NewString ] = MSG_LIST( ListHandle, 'Message', Text, 'wait' );
%
% Append the AppendString as Delimiter for a following Text to append (see below).
%
%   NewString = { [ sprintf(NumberForm,ListNumber)  Text AppendSAtring ] };
%
% [ ... ] = MSG_LIST( ListHandle, 'Message', Text, 'wait', AppendString );
%
%  use this new AppendString, instead of the default,
%   when the MessageListBox was created.
%
%
%-----------------------------------------------------------------------------
% 3. Append a Messages to the last one, Action 'MESSAGE' , Mode 'append'
%
% [ Msg, ListNumber, NewString ] = MSG_LIST( ListHandle, 'Message', Text, 'append' );
% 
% Append the Text to the last ListEntry, using the AppendString as Delimiter between.
% The NewString is build like:
%
%   NewString = { [ OldString  AppendString  Text ] };   % CellString !!!
%
%
% [ Msg, ... ] = MSG_LIST( ListHandle, 'Message', Text, 'append', AppendString );
%   
%  use this new AppendString, instead of the default,
%   when the MessageListBox was created.
%
%              
%-----------------------------------------------------------------------------
% 4. Reset the MessageListBox , Action 'RESET' 
% 
%  Msg = MSG_LIST( ListHandle, 'Reset' );
%
% Deletes all ListBoxEntries.
%
%
%-----------------------------------------------------------------------------
% 5. Resize the MessageListBox after the Figure was resized, Action 'RESIZE' 
%
%  [ Msg, PixelPosition, RowNumber ] = MSG_LIST( ListHandle, 'Resize', Mode )
%
%  Resize the MessageListBox in the Figure,  depending on Mode
%   call this from the ResizeFcnPropertie of the Figure.  
%
%  Mode:  'normalized'  ( short:  'n' )
%
%   The ListBoxPosition will set normalized to the Figure, with the normalized
%   Position when the ListBox was created. 
%   In this case the RowNumber may be changed.
%
%  Mode:  'absolut'     ( short:  'a' )
%
%   The ListBoxPosition will set absolut, 
%    defined by the original RowNumber and the PixelOffset.
%
%
%-----------------------------------------------------------------------------
% 6. Example, execute the commands Step by Step
%
% % Create a Figure
% fig = figure('units'   , 'pixels' , ...
%              'position', [100 300 200 200 ] ); 
%
%   % Create a MessageListBox 
% [ Msg , HList ] = msg_list( fig , 'New' , ...
%                                   'RowNumber'    , 5  , ...
%                                   'RecentNumber' , 10 , ... 
%                                   'PixelOffset'  , [ 60  10  5 ] , ...
%                                   'FontUnits'    , 'points' , ...
%                                   'FontSize'     , 8                );
%
% if ~isempty(Msg)
%   % An Error occured
%   fprintf([ Msg  char(10) ]);
%   return
% end
%
% % New Message
% [ Msg, nr, str ] = msg_list( HList, 'Message', 'Read Data' );
%
% % Append a Message
% [ Msg, nr, str ] = msg_list( HList, 'Message', 'ok', 'append' );
%
% % MultiLine Message
% [ Msg, nr, str ] = msg_list( HList, 'Message', ...
%                         [ 'Read Data from file:'  char(10)  'test.nc' ] );
%
% % Append a Message, NewLine as AppendString
% [ Msg, nr, str ] = msg_list( HList, 'Message', 'Error', ...
%                                'append', char(10) );
%  
% % Some More new Messages
%  for ii = 1 : 9
%   [ Msg, nr, str ] = msg_list( HList, 'Message', ...
%                             sprintf('Message %.0f',ii) );
%  end    
%    
%  % 11 Messages are now in the ListBox.
%  % You see, that the 1. Message is not visible, because the
%  %  RecentNumber was set to 10.
%
%  % Resize the Figure
%   set(fig,'position',[100 300 300 300 ])
%
%  % Resize the ListBox normalized, the RowNumber is changing
%    [Msg, pos, nrow] = msg_list( HList, 'Resize', 'normalized' );
% 
%  % Resize the ListBox absolut
%    [Msg, pos, nrow] = msg_list( HList, 'Resize', 'absolut' );
% 


Nout = nargout-1;

varargout    = cell(Nout,1);
varargout(:) = { [] };


Msg = '';
out = [];

nl = char(10);


Msg0 = ' MSG_LIST: ';

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

%------------------------------------------------------------
% Check Handle

if strcmp(action,'NEW')

 if ~strcmp(get(h,'type'),'figure')
   Msg = [ Msg0  'Handle before NEW-Action must be a Figure.' ];
   return
 end
 fig = h;
   h = [];

%------------------------------------------------------------
else

 ok = strcmp(get(h,'type'),'uicontrol');
 if ok
   ok = ( strcmp(get(h,'style'),'frame')      &  ...
          strcmp(get(h,'tag'),'MESSAGE_FRAME')          );
 end

 if ~ok
   Msg = [ Msg0  'Handle before ' action '-Action must be a ' ...
           'valid Message-ListBoxUIControl.' ];
   return
 end

 fig = get(h,'parent');

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

  is_win = strcmp(upper(computer),'PCWIN');
  
  ud = struct( 'RecentNumber'  , { 50 } , ...  % Number of stored Entries 
               'RowNumber'     , {  5 } , ...  % Number of Rows shown, determine High
               'PixelOffset'   , { zeros(1,3) } , ... % from [ Left Right Top ] of Figure
               'AppendString'  , { ' ... '    } , ... % String between if 'append'
               'NumberForm'    , { '%4.3d:  ' } , ... % Format for New EntryNumber, using sprintf
               'NumberBlank'   , { ''         } , ... % Blanks with same Width like EntryNumber
               'ActualMode'    , { ''         } , ... % Mode of Actual Entry
               'ActualNumber'  , {  0 }         , ... % Actual Number of ListEntries
               'LineNumber'    , {  0 }         , ... % LineNumber of RecentListEntries 
               'LineWidth'     , {  0 }         , ... % PixelLineWidth of RecentListEntries 
               'PixelHigh'     , {  0 }         , ... % PixelHigh of Font
               'BorderWidth'   , { 3-0*is_win } , ... % BorderPixelWidth of Listbox
               'SliderWidth'   , {15+3*is_win } , ... % SliderPixelWidth of Listbox
               'LineOffset'    , { 2-[2 0]*is_win } , ... % PixelOffset per ListboxLine
               'RowNumber0'    , {  5 }         , ... % Default Number of Rows
               'FrameHandle'   , { NaN }        , ... % Handle of Frame arround
               'ListHandle'    , { NaN }        , ... % Handle of ListBox
               'AxeHandle'     , { NaN }        , ... % Handle of Axe, contains TextObject
               'TextHandle'    , { NaN }        , ... % Handle of TextObject, using for width
               'PixelPosition' , { zeros(1,4) } , ... % Actual PixelPosition [ Left Bottom Width High ]
               'NormPosition'  , { zeros(1,4) } , ... % Default normalized Position
               'FigPos0'       , { figpos     }  );   % Default FigurePosition


  %-------------------------------------------------------------------
  % Look, if the first 4 Properties are given in Inputs

  fields = fieldnames(ud);
  fields = fields([1 2 3 4 5]);

  nv = size(VarIn,1);

  is_ud = zeros(nv,1);

  for ii = 1 : nv
     jj = find(strcmp(VarIn{ii,1},fields));
     if ~isempty(jj)

       val = VarIn{ii,2};
       msg = '';

       if strcmp(fields{jj},'AppendString')
         if ~ischar(val)  |  ( size(val,1) > 1 )  | ( ndims(val) > 2 )
           msg = [ 'Value for '  fields{jj} ' must be a String.' ];
         end
       elseif   strcmp(fields{jj},'NumberForm')
         if ~ischar(val)  |  ( size(val,1) > 1 )  | ( ndims(val) > 2 )
           msg = [ 'Value for '  fields{jj} ' must be a String, ' ...
                   'using as Format for SPRINTF.'  ];
         else
           ok = 1;
           try
             sprintf(val,0);
           catch
             ok=0;
           end

           if ~ok
             msg = [ 'Value for '  fields{jj} ' must be a valid Format for SPRINTF.' ];
           end             
         end
       else
         if ~isnumeric(val) | isempty(val)
            msg = [ 'Value for '  fields{jj} ' must be numeric and not empty.' ];
         else
            valmin = 1-strcmp(fields{jj},'PixelOffset');
            if ~all( mod(val,1) == 0 )  &  ...
                all( val >= valmin )
                msg = [ 'Value for '  fields{jj} ' must contain Integers >= ' ...
                         sprintf('%.0f',valmin) '.' ];
            end
            si0 = size(getfield(ud,fields{jj}));
            if ~isequal(size(val),si0)
               str = sprintf(' %.0f by',si0);
               str = [ '[' str(1:end-2) ']' ];
               msg = [  msg  nl(1:(end*(~isempty(msg)))) ...
                        'Value for '  fields{jj} ' must be of Size ' str '.' ];
            end
         end
       end
       % AppendString

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

      VarIn{ii,3} = get(fig,['DefaultUIControl' VarIn{ii,1}]);

    catch

       msg = lasterr;

    end

    if isempty(msg)
       set(fig,['DefaultUIControl' VarIn{ii,1}],VarIn{ii,2});
    else
       Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
               'Invalid UIControlProperty: ' VarIn{ii,1} ];
       break
    end

  end

  if ~isempty(Msg)
     % Set DefaultProperties back
    for jj = 1 : ii-1
      set(fig,['DefaultUIControl' VarIn{jj,1}],VarIn{jj,3});
    end

    Msg = [ Msg0  Msg ];
    return

  end

  % Define Axes, containing TextObject

  ud.AxeHandle = axes('parent'  , fig , ...
                      'units'   , 'pixels' , ... 
                      'position',[ 0 0 1 1 ] , ...
                      'xlim'    , [ 0  1 ] , ...
                      'ylim'    , [ 0  1 ] , ...
                      'xtick'   , [] , ...
                      'ytick'   , [] , ...
                      'visible' , 'off' , ...
                      'tag'     , 'MESSAGE_AXES' , ...
                      'nextplot' , 'add' , ... 
                      'handlevisibility' , 'callback'   );

  ud.TextHandle = text('parent' , ud.AxeHandle , ...
          'units' , 'pixels' , ...
          'position',[ 1 1 0 ] , ...
          'string'  , '' , ... 
          'interpreter' , 'none' , ...
          'fontunits'  , get(fig,'DefaultUIControlFontUnits') , ... 
          'fontsize'   , get(fig,'DefaultUIControlFontSize')  , ... 
          'fontname'   , get(fig,'DefaultUIControlFontName')  , ... 
          'fontangle'  , get(fig,'DefaultUIControlFontAngle') , ... 
          'fontweight' , get(fig,'DefaultUIControlFontWeight') , ...
          'visible'    , 'off' , ...
          'tag'        , 'MESSAGE_TEXT'       );

   % Determine NumberBlank
           set( ud.TextHandle , 'string' , sprintf( ud.NumberForm , 0 ) );
    ext1 = get( ud.TextHandle , 'extent' );
           set( ud.TextHandle , 'string' , char(32*ones(1,10))  );
    ext2 = get( ud.TextHandle , 'extent' );
           set( ud.TextHandle , 'string' , ''  );
      Nblank = ceil( ext1(3) / ( ext2(3) / 10 ) ) ;

   ud.NumberBlank = char(32*ones(1,Nblank));

 
   
   ud.FrameHandle = uicontrol('parent'   , fig     , ...
                              'visible'  , 'off'   , ...
                              'selected' , 'off'   , ...
                              'units'    , 'pixels', ...
                              'position' ,  [ 0  0  1  1 ] , ...
                              'style'    , 'frame' , ...
                              'tag'      , 'MESSAGE_FRAME'  );

                      
  % Define List first as Text to determine right PixelSize per Line

   ud.ListHandle = uicontrol('parent'   , fig   , ...
                             'visible'  , 'off'  , ...
                             'selected' , 'off' , ...
                             'units'    , 'characters', ...
                             'position' , [0 0 1 ud.RowNumber ] , ...
                             'style'    , 'text' , ...
                             'min'      , 0 , ...
                             'max'      , 1 , ...
                             'string'   , { '' } , ...
                             'tag'      , 'MESSAGE_LIST'  );


%  [ characters ]  --- get(0,'defaultuicontrolfont...') --> [ pixels ]

   set(ud.ListHandle,'units','pixels');

   pos = get(ud.ListHandle,'position');  % PixelPosition

   ud.PixelHigh = pos(4) / ud.RowNumber;

   LineHigh = ud.PixelHigh + ud.LineOffset(2);

   % New Position of ListBix vincl. Border & LineOffset
   %  Frame with ud.BorderWidth arround !!!

   pos(4) = ud.RowNumber * LineHigh + 4 * ud.BorderWidth + ...
            ( ud.SliderWidth - LineHigh ) * ( ud.SliderWidth > LineHigh );

   is_top = ( ud.PixelOffset(3) >= 0 );

   pos(1) = 1 + ud.PixelOffset(1);
   pos(2) = ( figpos(4) - pos(4) + 1 ) * is_top - ud.PixelOffset(3);
   pos(3) =   figpos(3) - pos(1) + 1            - ud.PixelOffset(2);


   ud.PixelPosition = pos;
   ud.NormPosition  = pos ./ figpos([ 3  4  3  4 ]);

   ud.RowNumber0    = ud.RowNumber;


   % Note: Frame arround !!!

    pos([3 4]) = pos([3 4]) + ( 1+2*ud.BorderWidth - pos([3 4]) ) .* ...
                              ( 1+2*ud.BorderWidth > pos([3 4]) );

   set( ud.ListHandle , 'position' , ( pos + [ 1  1   -2  -2 ] * ud.BorderWidth ) , ...
                        'style'    , 'listbox' , ...
                        'visible'  , 'on'             )


   set( ud.AxeHandle   , 'position' , pos        );


   set( ud.FrameHandle , 'position' , pos  , ...
                         'visible'  , 'on' , ...
                         'userdata' , ud           );

   % Set DefaultProperties back
   for ii = 1 : size(VarIn,1)
     set(fig,['DefaultUIControl' VarIn{ii,1}],VarIn{ii,3});
   end

   varargout{1} = ud.FrameHandle;
   varargout{2} = pos;


%*********************************************************
case 'MESSAGE'

  if isempty(VarIn)
     Msg = [ Msg0  'Input Message is undefined.' ];
     return
  end

  if ~iscellstr(VarIn)
     Msg = [  Msg0 'Inputs must be Strings: Message , Mode.' ];
     return
  end

  txt = VarIn{1};

  if ~( ischar(txt) & ( prod(size(txt)) == size(txt,2) ) )
     Msg = [  Msg0 'Input Message must be a String.' ];
     return
  end

   ud = get(h,'userdata');

  add = { 'new'  ud.AppendString };

   nv = size(VarIn,1);

  VarIn( nv+1 : 3 ) = add( (nv+1:3) - 1 );
 
  mode = lower(VarIn{2}(1));

  if ~any(strcmp(mode,{'n' 'a' 'w'}))
     Msg = [  Msg0  'Input Mode  must be ''new'', ''wait'' or ''append''.' ];
     return
  end

  app = VarIn{3};  % AppendString

  if ~( ischar(app) & ( prod(size(app)) == size(app,2) ) )
     Msg = [  Msg0 'Input AppendString must be a String.' ];
     return
  end


  is_new = any(strcmp(mode,{'n' 'w'}));

  if ~is_new  &  ( ud.ActualNumber == 0 )
     Msg = [  Msg0  'Nothing to append on.' ];
     return
  end
     
     str0 = get(ud.ListHandle,'string');


    % New ActualNumber

    ud.ActualNumber = ud.ActualNumber+is_new;
 

    % New String

    if is_new
      str =  [ sprintf(ud.NumberForm,ud.ActualNumber)    txt  ];
      if strcmp(mode,'w')
         str = [ str app ];
      end
    else   
      str =  [ str0{end}  app(1:(end*(~strcmp(ud.ActualMode,'w'))))  txt  ]; 
    end

    ud.ActualMode = mode;

    % Check String for CR
    str(find(double(str)==13)) = [];


    % Build MultiLine Strings if NewLine

    % 1. "'"     --> "''"
    str = strrep( str , char(39) , char([39  39]) ); 

    % 2. NewLine --> "';'"
    str = strrep(str,char(10),[ char([39  59  39])  ud.NumberBlank ]);

    str = eval( char([  '{''' str  '''}' ]) );


    % Determine StringWidth in Pixel, using ud.TextHandle

    set(ud.TextHandle,'string','', ...
          'units'      , 'pixels'           , ...
          'fontunits'  , get(h,'FontUnits') , ... 
          'fontsize'   , get(h,'FontSize')  , ... 
          'fontname'   , get(h,'FontName')  , ... 
          'fontangle'  , get(h,'FontAngle') , ... 
          'fontweight' , get(h,'FontWeight'), ...
          'visible'    , 'off'      );

     ns = size(str,1);

    ext = zeros(ns,4);

     ok = zeros(ns,1);

    for ii = 1 : ns
       ok(ii) = ~all( double(str{ii}) == 32 );
       if ok(ii)   
         set(ud.TextHandle,'string',str{ii});
         ext(ii,:) = get(ud.TextHandle,'extent');
       end
    end

    set(ud.TextHandle,'string','');

      ok  = find(ok);
      str = str(ok);
      ext = max(ext(ok,3));

     str0 = [ str0(1:end-1+is_new) ; str ];


      % Note: Frame arround ==> "4*ud.BorderWidth"  !!!

      ww  = ud.PixelPosition(3) - ud.LineOffset(1) - 4*ud.BorderWidth - ud.SliderWidth;

      % Check, if 1. Entrie to remove
      rm_line = ( ( ud.ActualNumber > ud.RecentNumber )  &  is_new );

           n0 = min( ud.ActualNumber-is_new , ud.RecentNumber );
           n1 = n0 + is_new - rm_line;
                 
          ind = [ ( 1+rm_line : n0-1+is_new )  n0+(n0==0) ];

      ud.LineNumber     =  ud.LineNumber(ind);
      ud.LineNumber(n1) = (ud.LineNumber(n1)-1)*(~is_new) + size(str,1);

      ud.LineWidth      =  ud.LineWidth(ind);
      ud.LineWidth(n1)  = max( ud.LineWidth(n1)*(~is_new) , ext );


            ln = sum(ud.LineNumber);  % Number of Lines in ListBox
      
          wmax = max(ud.LineWidth);
          ww0  = ud.PixelPosition(3) - ud.LineOffset(1) - 4*ud.BorderWidth;  % WithOut Slider
 
       right_slider = ( ln > ud.RowNumber - ( wmax > ww0 ) );

        bott_slider = ( wmax > ww0-ud.SliderWidth*right_slider );

            lt = ln-ud.RowNumber+1+bott_slider;       % ListBoxTop
            lt = lt + ( 1 - lt ) * ( lt < 1 );


     set( ud.ListHandle , 'string'     , str0( end-ln+1 : end ) , ...
                          'listboxtop' , lt   );

     set( ud.FrameHandle , 'userdata'   , ud );


     varargout{1} = ud.ActualNumber;
     varargout{2} = str0( end-ud.LineNumber(end)+1 : end );

     drawnow  % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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


  ud = get(h,'userdata');

  % FigurePosition in Pixels

  figuni = get(fig,'units');
           set(fig,'units','pixels');
  figpos = get(fig,'position');
           set(fig,'units',figuni);


   LineHigh = ud.PixelHigh + ud.LineOffset(2);
  %--------------------------------------------------------------------------
  if strcmp(mode,'n')
  % Normalized

          pos = ud.NormPosition .* figpos([3 4 3 4]);


         % Note: Frame arround ==> "4*ud.BorderWidth"  !!!
        
         pos4 = pos(4) - 4 * ud.BorderWidth - ...
              ( ud.SliderWidth - LineHigh ) * ( ud.SliderWidth > LineHigh );

       ud.RowNumber = floor(pos4/LineHigh);  % New Actual RowNumber
           
         pos4 = ud.RowNumber * LineHigh + 4 * ud.BorderWidth + ...
                ( ud.SliderWidth - LineHigh ) * ( ud.SliderWidth > LineHigh );
       
       pos(2) = pos(2) + pos(4) - pos4;
       pos(4) = pos4;

  %--------------------------------------------------------------------------
  else
  % Absolut

      ud.RowNumber =  ud.RowNumber0;

     % Note: Frame arround ==> "4*ud.BorderWidth"  !!!

     pos(4) = ud.RowNumber * LineHigh + 4 * ud.BorderWidth + ...
              ( ud.SliderWidth - LineHigh ) * ( ud.SliderWidth > LineHigh );

     is_top = ( ud.PixelOffset(3) >= 0 );

     pos(1) = 1 + ud.PixelOffset(1);
     pos(2) = ( figpos(4) - pos(4) + 1 ) * is_top - ud.PixelOffset(3);
     pos(3) =   figpos(3) - pos(1) + 1            - ud.PixelOffset(2);

  end
  %--------------------------------------------------------------------------

      pos              = ceil(pos);

      ud.PixelPosition = pos;


    % Determine new ListBoxTop

            ln = sum(ud.LineNumber);  % Number of Lines in ListBox
      
          wmax = max(ud.LineWidth);
          ww0  = ud.PixelPosition(3) - ud.LineOffset(1) - 4*ud.BorderWidth;  % WithOut Slider
 
       right_slider = ( ln > ud.RowNumber - ( wmax > ww0 ) );

        bott_slider = ( wmax > ww0-ud.SliderWidth*right_slider );

           lt = ln-ud.RowNumber+1+bott_slider;       % ListBoxTop
           lt = lt + ( 1 - lt ) * ( lt < 1 );

         pos([3 4]) = pos([3 4]) + ( 1+2*ud.BorderWidth - pos([3 4]) ) .* ...
                                   ( 1+2*ud.BorderWidth > pos([3 4]) );

     set( ud.ListHandle , 'units'      , 'pixels' , ...
                          'position'   , ( pos + [ 1  1   -2  -2 ] * ud.BorderWidth ) , ...
                          'listboxtop' , lt        );
      
     set( ud.AxeHandle , 'units'      , 'pixels' , ...
                         'position'   ,  pos     , ...
                         'visible'    , 'off'          );

     set( ud.FrameHandle , 'units'      , 'pixels' , ...
                           'position'   ,  pos     , ...
                           'userdata'   , ud            );
 

    varargout{1} = ud.PixelPosition;  % Actual PixelPosition
    varargout{2} = ud.RowNumber;      % Actual RowNumber


%*********************************************************
case 'RESET'

   ud = get(h,'userdata');

   ud.ActualNumber = 0;
   ud.LineNumber   = 0;
   ud.LineWidth    = 0;

   set( ud.ListHandle , 'string'     , {''} , ...
                        'value'      ,   1  , ... 
                        'listboxtop' ,   1         )


%*********************************************************
case 'VISIBLE'

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
   
   ud = get(h,'userdata');

   set([ ud.FrameHandle ud.ListHandle ] , 'visible' , sets  );

   set([ ud.AxeHandle   ud.TextHandle ] , 'visible' , 'off' );


end

