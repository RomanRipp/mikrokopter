function [ Msg , varargout ]  = browser(action,varargin);
 
% BROWSER  12. Jul. 2006 Creates and Manages a BrowserFigure
%
% Syntax: [ Msg , ... ] = BROWSER( Action , FigureHandle , ... )
%         [ Msg , ... ] = BROWSER( Action ,   ItemHandle , ... )
%         [ Msg , ... ] = BROWSER( Action , FigureHandle , ItemTag , ... )
%         
%-----------------------------------------------------------------
% 1. Create New BrowserFigure 
% 
%  [ Msg , BrowserFig ] = BROWSER( 'New' , FigureName , LogoText );
%
%  Creates a new BrowserFigure and returns it's Handle.
%  
%    FigureName is the Name of this Figure.
%    LogoText   is a 3-Element CellString for the Logo, used by MSG_LOGO
%
%    Msg        contains ErrorMessages, empty  if all is ok.
%
%  additional Inputs:
%
%  [ Msg, BrowserFig, IH, IT, IS ] = BROWSER( 'New' , FigureName , LogoText , ...
%                                             ItemStructure );
%   
%   add the Items defined by the ItemStructure to the ContentsMenu
%    of the BrowserFigure, see below for more Infos about the ItemStructure.
%
%-----------------------------------------------------------------
% 2. Add Items to the ContentsMenu  or  existing ItemMenus
% 
%  [ Msg, IH, IT, IS ] = BROWSER( 'Add' , BrowserFigure , ItemStructure );
%
%  Adds Items, defined in the ItemStructure, to the ContentsMenu
%   of the BrowserFigure. 
%
%  The ItemStructure contains the Fields:
%  
%   FieldName            Value
%
%    'Tag'         String for the Tag of the Item
%    'Label'       String for the ItemLabel
%    'Titel'       String for the Header of the Text
%    'Text'        CharArray or CellStringArray for the Text
%    'File'        String for a FileName, contains the Text
%    'Note'        CharArray or CellStringArray for a FootNote
%    'Matlab'      String for an MatlabEvent
%    'Children'    ItemStructure for a SubItem
%
%  Note:
%
%    Any of the Fields 'Tag' , 'Label' or 'Titel'    must be defined 
%    and not empty!
%      If Tag   is empty, Label or Titel will used as Tag.
%      If Label is empty, Tag   or Titel will used as Label.
%      If Titel is empty, Label or Tag   will used as Titel.
% 
%    Any of the Fields 'Text', 'File', 'Matlab'  or 'Children' must 
%    be defined and not empty!
%
%    If the "Tag"   starts with '_', it will concatinated 
%           recursivly with the Tag  of the Parent.
%
%    If the "Titel" starts with '_', it will concatinated
%           with the Titel of the Parent, separated with ', '.
%
%    The Text of an Item, displayed in the Browser, will concatinated 
%      from the Titel, Text, File and the Text of the Children.
%
%    The TAB-Character (char(9)) in the Text of an Item (from Text and/or File)
%       will replaced with 8 Blanks.
%
%    The '\n'-Characters in "Text"  or  "Note"  will replaced 
%       with the NewLine-Character (char(10)).
%
%    In Text, File or Note can marked References to other Sections or Chapter
%     with: "@# ItemTag #@", where ItemTag is the String of the Tag of the 
%     Chapter or Section to refer.
%
%    Lines in File, starting with "@@", are Comments and will not added to the 
%     Text.
%  
%    If the String for "Text", "File" or "Matlab" starts with 
%       'MatlabHelp', the Text of the Item will generated using HELP, 
%       the following String is the Input for HELP, example:  
%        'MatlabHelpText  browser'  gets the Text from help('browser') and
%        displays the Help for this File (BROWSER.M)
%              
%    The String for "Matlab" has to start with:
%      - 'MatlabHelp', see above
%      - 'HelpWin'    opens HELPWIN, Inputs can follow 
%      - 'Eval'       calls EVAL, following Input required
%      - 'Image'      opens a Figure for an Image,
%                      following FileName required.
%              
%  
%  OutPuts:
%
%   IH   ItemHandle: Handles of ItemMenus, used to specify an Item
%   IT   ItemTag:    Tag's   of ItemMenus, used to specify an Item 
%                                             together with BrowserFigure
%   IS   Full ItemStructure in BrowserFigure
%
%
%  Use:
%
%  [ Msg, IH, IT, IS ] = BROWSER( 'Add' , ItemHandle , ItemStructure );
%
%    or ...
%
%  [ Msg, IH, IT, IS ] = BROWSER( 'Add' , BrowserFigure , ItemTag , ItemStructure );
%
%   to add an SubItem to an existing Item.
%
%
%  Add a DirectoryStructure:
%
%    If the Value for ItemStructure is an DirectoryName, specify a Directory
%     which contains an BROWSER.INI - File, the Structure specified by the
%     by the BROWSER.INI - File(s) will add to the ContentsMenu or an Item.
%   
%    Type  >> [ Msg , Example ] = browser('INI')
%
%      for more Informations about the StructureDefinition in BROWSER.INI,
%
%          >> [ Msg , Example ] = browser('INI',File)
%  
%      saves an example for an BROWSER.INI in File.
%        
%
%-----------------------------------------------------------------
% 3. Get the Text of an Item
%
%  [ Msg , CellText , V , Hist , Text ] = BROWSER( 'Text' , BrowserFigure )
%  [ Msg , CellText , V , Hist , Text ] = BROWSER( 'Text' , ItemHandle    )
%  [ Msg , CellText , V , Hist , Text ] = BROWSER( 'Text' , BrowserFigure , ItemTag )
%
%  returns the Text of the Item and all Children as
%     - CellStringArray CellText,
%     - CharArray       Text, 
%
%  the TitleHistory  Hist = { Nr  Titel },
%
%   and the Structure V with the Fields:
%
%  V.Nr         Nr. of the Item in the full ContentsStructure
%  V.Titel      Titel of Item
%  V.Text       Text  of Item
%  V.Note       FootNote
%  V.Matlab     Matlab
%  V.Reference  2-Column CellArray of References in Text: { String  Nr }
%  V.Image      FileName for Image
%  V.Children   Structure for Children
%
%-----------------------------------------------------------------
% 4. Activate an Item
%
%  [ Msg , ActiveHandle , ActiveTag , Text ] = BROWSER( 'Activate' , ItemHandle )
%  [ Msg , ActiveHandle , ActiveTag , Text ] = BROWSER( 'Activate' , ...
%                                                        BrowserFigure , ItemTag )
%
%  Activates the Item, specified by ItemHandle  or  BrowserFigure & ItemTag.
% 
%  the Handle, Tag  and Text of the active Item, displayed in the Browser, 
%    will returned.
%
% NOTE:  When an Item will Activated, the Children of the ContentsMenu will
%         rearranged. After this on some X-Server the Items are not visible.
%         Please use the <Refresh>-Menu to have a correct Visibility of the Items.
%
%-----------------------------------------------------------------
% 5. DisActivate the activated Item
%
%  [ Msg , ActiveHandle , ActiveTag ] = BROWSER( 'DisActivate' , BrowserFigure );
%
%   DisActivates the BrowserFigure and returns the Handle and Tag of the 
%     last active Item.
%
%-----------------------------------------------------------------
% 6. Create an ItemMenu in another Figure or Menu
%
%  [ Msg , MenuHandle ] = BROWSER( 'Menu' , BrowserFigure , Parent , Label )
%
%   Creates a copy of the ContentsMenu of the BrowserFigure
%     to the Figure or Menu, specified by FigureHandle.
%
%  The Handle of the new Menu is returned.
%  The new Menu gets the Label, specified by Label (optional).
%
%  [ Msg , MenuHandle ] = BROWSER( 'Menu' , ItemHandle              , ...
%                                   Parent , Label )
%  [ Msg , MenuHandle ] = BROWSER( 'Menu' , BrowserFigure , ItemTag , ...
%                                   Parent , Label )
%
%  Creates a copy of the ItemMenu, specified by  ItemHandle  or 
%    BrowserFigure & ItemTag.
%
%  Now a Item can activated from another Figure or Menu  using the new Menu.
%
%-----------------------------------------------------------------
% 7. Recurse trough the History
%
%  [ Msg , ActiveHandle , ActiveTag ] = BROWSER( 'Back'    , BrowserFigure );
%  [ Msg , ActiveHandle , ActiveTag ] = BROWSER( 'Forward' , BrowserFigure );
%     
%  jumps back or forward and returns the Handle and Tag of the active Item.
%
%-----------------------------------------------------------------
% 8. Reset
%
%  Msg = Browser( 'Reset' , BrowserFigure )
%
%   deletes all ContentsItems.
%
%-----------------------------------------------------------------
% 9. Close and Delete
%
% The command:  >> close(BrowserFigure)
%
%  sets the Visibility of the BrowserFigure to 'off'.
%
%         Use:  >> set(BrowserFigure,'visible','on')
%
%   or activate an Item to see the BrowserFigure again.
%
%         Use: >> delete(BrowserFigure)   to delete it. 
%  
%   All Children of Copies of ItemMenus in other Figures will deleted too.
%
%-----------------------------------------------------------------
% 10. Required M-Files
% 
%  MSG_LOGO (TEXTLOGO) , MSG_LIST
%
%  BRWSHELP, BRWSINI
%
%-----------------------------------------------------------------
% 11. Example
%
%  Type  >> [ Msg , ItemStructure ] = browser('Help')  
%
%   to see a example
%
%-----------------------------------------------------------------
%


Msg = '';
out = [];

Msg0 = 'BROWSER: ';

nl = char(10);

Nin  = nargin;
Nout = nargout - 1;

Nout = Nout * ( Nout > 0 );

varargout = cell(1,Nout);


if Nin < 1
  Msg = [ Msg0  'Not enough InputArguments.' ];
  return
end
 

%*****************************************************************
% Check Action-input

if ~( ischar(action)  &  ~isempty(action)  &  ...
      ( prod(size(action)) == size(action,2) ) ) 
  Msg = [ Msg0  'Action must be a nonempty String.' ];
  return
end

action = upper(action);

Msg0 = [ 'Browser( ' action ' ): ' ];


%*****************************************************************
% Catch HELP any INI

switch action

 %----------------------------------------------------------------
 case 'HELP'

  out = { [] };

  try 

   echo('brwshelp','on');

    [Msg,out{1}] = brwshelp;

   echo('brwshelp','off');

  catch

    Msg = [ Msg0 'Error using brwshelp.' nl lasterr ];

  end

    n = min(Nout,size(out,2));

    varargout(1:n) = out(1:n);

    return

 %----------------------------------------------------------------
 case 'INI'

   txt = help('brwsini');

   out = { txt };

   n = min(Nout,size(out,2));

   varargout(1:n) = out(1:n);

   txt = strrep(txt,'%','%%');
   txt = strrep(txt,'\','\\');

   if ( Nin >= 2 )
   % Write into FileName

     file = varargin{1};

     if ~( ischar(file)  &  ~isempty(file)  &  ...
         ( prod(size(file)) == size(file,2) ) ) 

          Msg = [ Msg0  'File must be a nonempty String.' ];

     else

        fid = fopen(file,'wt')

        if fid == -1
          
          Msg = [ Msg0  'Can''t open File:'  file  ];
 
        else

          txt = strrep(txt,'%','%%');
          txt = strrep(txt,'\','\\');

          fprintf(fid,txt);

          fclose(fid);

          return          

        end

     end

   end

   if Nout == 0
     fprintf(nl)
     fprintf(txt)
     fprintf(nl)
   end

   return

end
% switch

%*****************************************************************
% Check following Inputs

%-----------------------------------------------------------------
if Nin < 2
  Msg = [ Msg0  'Not enough InputArguments.' ];
  return
end

 
%-----------------------------------------------------------------

VarArg = varargin;

VarArg = VarArg(:);


%-----------------------------------------------------------------
if ~any( strcmp( action , { 'NEW' } ) );

  % Second Input must be a BrowserFigure  or  BrowserMenu

  fig = VarArg{1};

  VarArg = VarArg(2:end);

  [ ok , typ ] = checkin( fig , { 'figure' 'uimenu' } );

  if ~ok
    Msg = [ Msg0  'Second Input must be a Handle of '  ...
                  'a BrowserFigure or BrowserMenu.'        ];
    return
  end


  ud = get(fig,'userdata');

  %----------------------------------------------
  % Actions, requires BrowserFigure

  fig_act = { 'RESIZE'     'MINIMIZE'  ...
              'BACK'       'FORWARD'   ...
              'RESET'      'REFRESH'   ...
              'CLOSE'      'DELETE'    ...
              'MESSAGE'    'MSGBOX'    ...
              'PROBLEM'    'DISACTIVATE'   'SAVE'    };

  % Actions, requires BrowserMenu  or 
  %                   BrowserFigure and/or Tag 
  menu_act = { 'ADD'  'ACTIVATE'  'MENU'  'TEXT' };

  %----------------------------------------------
  if     any( strcmp( action , fig_act ) )
  % BrowserFigure required 

     ok = strcmp( typ , 'figure' );
     if ~ok

       Msg = [ Msg0  'Second Input must be a Handle of a BrowserFigure.' ];
       return

     end

  %----------------------------------------------
  elseif any( strcmp( action , menu_act ) )
  % BrowserFigure and/or Tag  or  BrowserMenu required 

     if strcmp( typ , 'uimenu' )

       % BrowserMenu

        mud = ud;
        fig = mud.Figure;
         ud = get( fig , 'userdata' );

     else 

       % BrowserFigure, check for Tag

       if isempty(VarArg)  &  ~strcmp(action,'TEXT')

          Msg = [ Msg0  'Not enough Input Arguments.' ];
          return

       elseif ~isempty(VarArg)

         tag    = VarArg{1};

         tag_ok = ( ischar(tag)  &  ~isempty(tag)  &  ...
                   ( prod(size(tag)) == size(tag,2) ) );

         if  ~tag_ok  &  any(strcmp( action , { 'ACTIVATE'  'TEXT' } ))
            Msg = [ Msg0  'Input Tag must be a nonempty String.' ];
            return
         end

         if tag_ok
            tag_ok = any( strcmp( tag , ud.Children.Tag ) );
            if  ~tag_ok  &  any(strcmp( action , { 'ACTIVATE'  'TEXT' } ))
               Msg = [ Msg0  'Invalid Tag.' ];
               return
            end
         end

         if tag_ok

            ii = find( strcmp( tag , ud.Children.Tag ) );

           if ( size(ii,1) > 1 )  &  ~strcmp( action , 'TEXT' )
              msg = [ Msg0 'Warning, Duplicate Items with Tag "'  tag '" !' ];
              browser('Message',fig,msg,'new');
              ii = ii(1);
           end

           mud = get( ud.Children.Handle(ii) , 'userdata' );
           typ = 'uimenu';
 
           VarArg = VarArg(2:end);

         end

       end

     end

  end
  %----------------------------------------------

end


switch upper(action)

%*******************************************************************
case 'NEW'


 if isempty(VarArg) | size(VarArg,1) < 2
   Msg = [ Msg0  'Inputs Name and LogoText are missing.' ];
   return
 end
   
 name = VarArg{1};

 if ~( ischar(name)  & ( prod(size(name)) == size(name,2) ) )
   Msg = [ Msg0  'Name must be a String.' ];
   return
 end


 %-----------------------------------------------------------------
 % Create Figure

 % Get Time first, used for Tag

        Time = floor(clock); % [ Year Month Day Hour Min Sec ]

        form = '%2.2d/%2.2d/%4.0f %2.0f:%2.2d:%2.2d';
         ind = [3 2 1 4 5 6];

         tag = sprintf(['Browser    ' form],round(Time(ind)));
      msgtag = sprintf(['BrowserMsg ' form],round(Time(ind)));


 
 fig = figure('menubar','none'        , ...
              'toolbar','none'        , ...
              'numbertitle' , 'off'   , ...
              'color'  , [ 1  1  1 ]  , ...
              'Name'   , name         , ...
              'tag'    , tag          , ...
              'visible'          , 'off'      , ...
              'createfcn'        , ''         , ...
              'integerhandle'    , 'on'       , ...
              'handlevisibility' , 'callback'        );


 %-----------------------------------------------------------------
 % UIMenu's

  Fig = epsstr(fig);

  bl1 = char( 32 * ones(1,2) );
  bl2 = char( 32 * ones(1,4) );

  HContents = uimenu( 'parent'   , fig        , ...
                      'label'    , 'Contents' , ...
                      'callback' , ''         , ...
                      'enable'   , 'on'       , ...
                      'visible'  , 'on'       , ...
                 'interruptible' , 'off'      , ...
                 'busyaction'    , 'cancel'          );
 

              uimenu( 'parent'   , fig    , ...
                      'label'    , bl1    , ...
                      'enable'   , 'off'  , ...
                      'visible'  , 'on'          );


  HBack     = uimenu( 'parent'   , fig      , ...
                      'label'    , '<<'     , ...
                      'callback' , ['browser(''Back'','  Fig ');' ] , ...
                      'enable'   , 'off'    , ...
                      'visible'  , 'on'     , ...
                 'interruptible' , 'off'    , ...
                 'busyaction'    , 'cancel'       );


  HForward  = uimenu( 'parent'   , fig      , ...
                      'label'    , '>>'     , ...
                      'callback' , ['browser(''Forward'','  Fig ');' ] , ...
                      'enable'   , 'off'    , ...
                      'visible'  , 'on'     , ...
                 'interruptible' , 'off'    , ...
                 'busyaction'    , 'cancel'       );

              uimenu( 'parent'   , fig   , ...
                      'label'    , bl2   , ...
                      'enable'   , 'off' , ...
                      'visible'  , 'on'        );

  HSave    =  uimenu( 'parent'   , fig      , ...
                      'label'    , 'Save'   , ...
                      'callback' , ''       , ...
                      'enable'   , 'on'     , ...
                      'visible'  , 'on'     , ...
                 'interruptible' , 'off'    , ...
                 'busyaction'    , 'cancel'       );

              uimenu( 'parent'   , fig   , ...
                      'label'    , bl2   , ...
                      'enable'   , 'off' , ...
                      'visible'  , 'on'        );

              uimenu( 'parent'   , fig         , ...
                      'label'    , 'Refresh'   , ...
                      'callback' , ['browser(''Refresh'','  Fig ');' ] , ...
                      'enable'   , 'on'     , ...
                      'visible'  , 'on'     , ...
                 'interruptible' , 'off'    , ...
                 'busyaction'    , 'cancel'       );

              uimenu( 'parent'   , fig   , ...
                      'label'    , bl2   , ...
                      'enable'   , 'off' , ...
                      'visible'  , 'on'        );

              uimenu( 'parent'   , fig      , ...
                      'label'    , '???'    , ...
                      'callback' , ['browser(''Problem'','  Fig ');' ] , ...
                      'enable'   , 'on'     , ...
                      'visible'  , 'on'     , ...
                 'interruptible' , 'off'    , ...
                 'busyaction'    , 'cancel'       );

              uimenu( 'parent'   , fig   , ...
                      'label'    , bl2   , ...
                      'enable'   , 'off' , ...
                      'visible'  , 'on'        );


  HClose    = uimenu( 'parent'   , fig      , ...
                      'label'    , 'Close'  , ...
                      'callback' , ['browser(''Close'','  Fig ');' ] , ...
                      'enable'   , 'on'     , ...
                      'visible'  , 'on'     , ...
                 'interruptible' , 'off'    , ...
                 'busyaction'    , 'cancel'       );


              uimenu( 'parent'   , HSave             , ...
                      'label'    , 'Text of Section' , ...
                      'callback' , ['browser(''Save'','  Fig ',''Section'');' ]     , ...
                      'enable'   , 'off'    , ...
                      'visible'  , 'on'     , ...
                 'interruptible' , 'off'    , ...
                 'busyaction'    , 'cancel'      );

              uimenu( 'parent'   , HSave             , ...
                      'label'    , 'Text of Chapter' , ...
                      'callback' , ['browser(''Save'','  Fig ',''Chapter'');' ]     , ...
                      'enable'   , 'off'    , ...
                      'visible'  , 'on'     , ...
                 'interruptible' , 'off'    , ...
                 'busyaction'    , 'cancel'      );

              uimenu( 'parent'   , HSave       , ...
                      'label'    , 'Full Text' , ...
                      'callback' , ['browser(''Save'','  Fig ',''All'');' ]     , ...
                      'enable'   , 'on'     , ...
                      'visible'  , 'on'     , ...
                 'interruptible' , 'off'    , ...
                 'busyaction'    , 'cancel'      );

 
 %-------------------------------------------
 % MSG_LOGO in Top !!!

 [Msg,HMsgFrame,HLogoFrame,HLogoText] = msg_logo(fig,'New',VarArg{2});

 if ~isempty(Msg)
   Msg = [ Msg0 'Error using MSG_LOGO( New ).' nl Msg ];
   delete(fig);
   return
 end

 %-------------------------------------------

 mud = get(HMsgFrame,'userdata');

 HMsgList = mud.ListHandle;

 ListPosition = get(HMsgList,'Position');

  BorderWidth = mud.BorderWidth;

 %-------------------------------------------

 lud = get(HLogoFrame,'userdata');

 Offset = lud.Offset.Msg(3);   % vertical Offset


 %-------------------------------------------
 % 2 TextControls for Chapter and Section

          str = 'Please select an Item in the Contents-Menu !';

 HSectionText = uicontrol('parent'          , fig        , ...
                          'units'           , 'pixels'   , ...
                          'style'           , 'text'     , ...
                          'string'          ,  str       , ...
                          'userdata'        ,  str       , ...
                          'foregroundcolor' , [ 0 0 0 ]  , ...
                          'backgroundcolor' , [ 1 1 1 ]  , ...
                          'FontName'        , get(HLogoText(2),'FontName'  ) , ...
                          'FontUnits'       , get(HLogoText(2),'FontUnits' ) , ...
                          'FontSize'        , get(HLogoText(2),'FontSize'  ) , ...
                          'FontWeight'      , get(HLogoText(2),'FontWeight') , ...
                          'horizontalalignment' , 'center' , ...
                          'tag'                 , 'BrowserSectionText' , ...
                          'visible'             , 'on'           );

 HChapterText = uicontrol('parent'          , fig           , ...
                          'units'           , 'pixels'      , ...
                          'style'           , 'text'        , ...
                          'string'          ,  name         , ...
                          'userdata'        ,  name         , ...
                          'foregroundcolor' , [ 0 0 0 ]     , ...
                          'backgroundcolor' , [ 1 1 1 ]     , ...
                          'FontName'        , get(HLogoText(1),'FontName'  ) , ...
                          'FontUnits'       , get(HLogoText(1),'FontUnits' ) , ...
                          'FontSize'        , get(HLogoText(1),'FontSize'  ) , ...
                          'FontWeight'      , get(HLogoText(1),'FontWeight') , ...
                          'horizontalalignment' , 'center' , ...
                          'tag'                 , 'BrowserChapterText' , ...
                          'visible'             , 'on'      );



 %-------------------------------------------
 % New MsgList at Bottom

 % Use same Inputs like MSG_LOGO

  MsgIn = { 'RowNumber'           mud.RowNumber             
            'RecentNumber'        100             
            'AppendString'        mud.AppendString         
            'PixelOffset'         ( Offset * [ 1  1 -1 ] )     
            'NumberForm'          mud.NumberForm         
            'FontSize'            get(HMsgList,'fontsize'  )              
            'FontName'            get(HMsgList,'fontname'  )       
            'FontWeight'          get(HMsgList,'fontweight')        
            'ForeGroundColor'     get(HMsgList,'foregroundcolor')     
            'BackGroundColor'     get(HMsgList,'backgroundcolor')   };


    MsgIn = permute(MsgIn,[2 1]);

    [Msg,HMessgFrame] = msg_list( fig , 'New' , MsgIn{:} );

    if ~isempty(Msg)
      Msg = [ Msg0 'Error using MSG_LIST( New ).' nl Msg ];
      delete(fig);
      return
    end

    mud = get( HMessgFrame , 'userdata' );    

    HMessgList =  mud.ListHandle;

    CB  = ['browser(''Minimize'','  Fig ');' ];
 
    mud = struct( 'Visible'     , { 1 }               , ...
                  'PixelHigh'   , { mud.PixelHigh   } , ...
                  'LineOffset'  , { mud.LineOffset  } , ...
                  'SliderWidth' , { mud.SliderWidth }        );

    set( HMessgList , 'callback'      , CB   , ...
                      'tooltipstring' , 'Click to Minimize'        );
    

    HMessgPush = uicontrol( 'parent'   , fig          , ...
                            'style'    , 'pushbutton' , ...
                            'units'    , 'pixels'     , ... 
                            'string'   , ''           , ...
                            'callback' , CB           , ...
                            'visible'  , 'off'        , ...
                            'userdata' , mud          , ...   
                            'backgroundcolor' , get(HMessgFrame,'backgroundcolor') , ... 
                            'foregroundcolor' , get(HMessgFrame,'foregroundcolor') , ... 
                            'tooltipstring'   , 'Click to Maximize'  , ...
                            'interruptible' , 'off'   , ...
                            'busyaction'    , 'cancel'      );
    

 %-------------------------------------------
 % ListBox for HelpText
 % DefaultHelpText created by function INTRO below
 
   HHelpFrame = uicontrol('parent'          , fig           , ...
                          'units'           , 'pixels'      , ...
                          'style'           , 'frame'       , ...
                          'foregroundcolor' , [ 0 0 0 ]     , ...
                          'backgroundcolor' , [ 1 1 1 ]     , ...
                          'visible'             , 'on'    , ...
                          'interruptible'       , 'off'   , ...
                          'busyaction'          , 'cancel'     );

  scr_uni = get(0,'units');      set(0,'units','pixels')
  scr_si  = get(0,'ScreenSize'); set(0,'units',scr_uni);

  is_win = strcmp(upper(computer),'PCWIN');
  if is_win
     FontName = 'courier';
  else
     FontName = { 'arrial'  get(0,'fixedwidthfontname') };
     FontName = FontName{ 1 + ( scr_si(4) >= 1050 ) } ;
  end 

    HHelpList = uicontrol('parent'          , fig           , ...
                          'units'           , 'pixels'      , ...
                          'style'           , 'listbox'     , ...
                          'string'          ,  intro        , ...
                          'userdata'        ,  intro        , ...
                          'foregroundcolor' , [ 0 0 0 ]     , ...
                          'backgroundcolor' , [ 1 1 1 ]     , ...
                          'FontName'        , FontName                       , ...
                          'FontUnits'       , get(HLogoText(2),'FontUnits' ) , ...
                          'FontSize'        , get(HLogoText(2),'FontSize'  ) , ...
                          'FontWeight'      , 'normal'    , ...
                          'horizontalalignment' , 'left'  , ...
                          'visible'             , 'on'    , ...
                          'interruptible'       , 'off'   , ...
                          'busyaction'          , 'cancel'     );
   

 %-------------------------------------------
 % LogoMsgList --> Frame

 set( HMsgList , 'style'           , 'frame'     , ...
                 'backgroundcolor' , [ 1  1  1 ]       );

 %-------------------------------------------
 % Minimum FigurePosition

      LogoPosition = get( HLogoFrame  , 'position' );
     MessgPosition = get( HMessgFrame , 'position' );


   MinWidth = sum( lud.Offset.Logo([1 2]) , 2 ) + ...
              sum( lud.Offset.Msg([1 2])  , 2 ) + ...
                3 * LogoPosition(3);

   MinHeight = 4 * Offset + 3*LogoPosition(4) + MessgPosition(4);
 

 %-------------------------------------------
 % Build UserData

 e0 = zeros(0,0); 
 e1 = zeros(0,1);        % Empty matrix: 0-by-1
 e2 = zeros(1,0);        % Empty matrix: 1-by-0

 ec =  cell(0,1);        % Empty cell array: 0-by-1

 es = char(zeros(1,0));  % Empty string: 1-by-0


 Menu   = struct( 'Contents' , { HContents } , ...
                  'Back'     , { HBack     } , ...
                  'Forward'  , { HForward  } , ...
                  'Save'     , { HSave     } , ...
                  'Close'    , { HClose    }         );
  
 Handle = struct( 'LogoFrame'   , { HLogoFrame   } , ...
                  'LogoText'    , { HLogoText    } , ...
                  'MsgFrame'    , { HMsgFrame    } , ...
                  'MsgList'     , { HMsgList     } , ...
                  'ChapterText' , { HChapterText } , ...
                  'SectionText' , { HSectionText } , ...
                  'MessgFrame'  , { HMessgFrame  } , ...
                  'MessgList'   , { HMessgList   } , ...
                  'MessgPush'   , { HMessgPush   } , ...
                  'HelpList'    , { HHelpList    } , ...
                  'HelpFrame'   , { HHelpFrame   }        );

 % All ContentsMenus's

 Children = struct( 'Handle' , { e1 } , ...     % All Items
                    'Tag'    , { ec } , ...     % their Tags
                    'Menu'   , { e1 }        ); % Menu's in other Figures


% Special Values for FileName, to get HelpText or launch HelpWin

 Matlab = struct( 'Help'  , {   'MatlabHelp'                } , ...
                  'Event' , { { 'HelpWin'  'Eval'  'Image' } }  );

 % 'Help'  used for 'Text' and 'File'  or  Matlab
 % 'Event' used for 'Matlab'


 Color = struct( 'Default'  , { get( 0 , 'DefaultUIMenuForeGroundColor' ) } , ...
                 'Selected' , { [ 0  0  0.5 ] }  );

 Character = struct( 'Comment'  , { '@@' }  , ...
                     'NewLine'  , { '\n' }  , ...
		     'TitelSep' , { '='  }  , ...
                     'NoteSep'  , { '-'  }  , ...
                     'RefSep'   , { ', ' }  , ...
                     'RefText'  , { { 'Section'  'Chapter' } }  , ...
                     'RefMark'  , { { '@#'  '#@' } }  );
		     
 % 'Comment'   used for 'File', hide Lines starting with this
 % 'NewLine'   used for 'Text' and 'Note'
 % 'Separator' SeperatorLine below Titel in BrowserText
  

 ud = struct( 'ID'          , { fig }         , ...
              'Tag'         , { tag }         , ...
              'Children'    , { Children }    , ...
              'MsgTag'      , { msgtag }      , ...
              'Handle'      , { Handle }      , ...
              'Menu'        , { Menu   }      , ...
              'Matlab'      , { Matlab }      , ...
              'BorderWidth' , { BorderWidth } , ...
              'Offset'      , { Offset      } , ...
              'MinPosition' , { [ MinWidth MinHeight ] } , ...
              'Color'       , { Color }        , ...
              'Character'   , { Character }              );
	      

  cud = struct( 'Children'   , { e1 } , ...  % ID's of existing Chapters
                'ActiveRoot' , { e0 } , ...  % ID   of actual activated Chapter
                'ActiveID'   , { e0 } , ...  % ID   of actual activated Menu
                'ActiveTag'  , { es } , ...  % Tag  of actual activated Menu
              'SelectedTag'  , { es } , ...  % Tag  of actual selected  Menu
                'History'    , { e2 } , ...  % ID's of selected Menu's (History)
                'Value'      , {  0 }   ...  % Value in History
                                            );

  set( HContents , 'userdata' , cud );


 %-------------------------------------------

 Fig = epsstr(fig);

 ResizeFcn = ['browser(''Resize'' ,' Fig ');' ];
  CloseFcn = ['browser(''Close''  ,' Fig ');' ]; 
 DeleteFcn = ['browser(''Delete'' ,' Fig ');' ]; 

 ch = get( fig , 'children' );

 set(fig,'userdata'        , ud        , ...
         'ResizeFcn'       , ResizeFcn , ...
         'CloseRequestFcn' , CloseFcn  , ...
         'DeleteFcn'       , DeleteFcn , ...
         'children'        , ch               );

 % Update Children too, better using PCWIN

 %-------------------------------------------
 % Set correct Size of UIControls

 Msg = browser('Resize',fig);

 if ~isempty(Msg)
   Msg = [ Msg0 'Error using BROWSER( Resize ).' nl Msg ];
   delete(fig);
   return
 end

  out    = cell(1,4);
  out(1) = { fig };

 %-------------------------------------------
 % Add Contents

 
 if size(VarArg,1) > 2

   [Msg,out{2},out{3},out{4}] = browser('Add',fig,VarArg{3:end});

   if ~isempty(Msg)
     Msg = [ Msg0 'Error using BROWSER( Add ).' nl Msg ];
   end

   if isempty(out{4})
     delete(fig);
     return
   end

 end

 %-------------------------------------------

   n = min(Nout,size(out,2));

   varargout(1:n) = out(1:n);



%*******************************************************************
case 'RESIZE'

  %------------------------------------------------------------------
  % FigurePosition in Pixels
  % Check with MinPosition

    figuni = get( fig , 'units'    );
    ResFcn = get( fig , 'ResizeFcn');

             set( fig , 'units'     , 'pixels', ... 
                        'resizefcn' , ''            );

    figpos = get( fig , 'position' );

    figtop = figpos(2)+figpos(4);

    figpos([3 4]) = figpos([3 4]) + ( ud.MinPosition - figpos([3 4]) ) .* ...
                                    ( ud.MinPosition > figpos([3 4]) );

    figpos(2) = figtop - figpos(4);  % UpperLeft Corner fixed

 
    set( fig , 'position'  , figpos );
    set( fig , 'units'     , figuni );
    set( fig , 'ResizeFcn' , ResFcn );


  %------------------------------------------------------------------
  % Resize Objects
      
  msg_logo(ud.Handle.LogoFrame  ,'Resize','absolut');

  msg_list(ud.Handle.MessgFrame ,'Resize','absolut');

  mud = get(ud.Handle.MessgList,'userdata');

  %------------------------------------------------------------------
  % Get new Positions

  LogoPosition = get(ud.Handle.LogoFrame  , 'position' );

  TextPosition = get(ud.Handle.LogoText,'position');
  TextPosition = cat(1,TextPosition{:});

  ListPosition = get(ud.Handle.MsgList,'position');

   MsgPosition = get(ud.Handle.MessgFrame , 'position' );


  if strcmp( get(ud.Handle.MessgPush,'visible') , 'on' )

    MsgPosition = get( ud.Handle.MessgPush , 'position' );

  end
 
  %--------------------------------------------------
  % Fit Chapter- and SectionText into MsgList

  pos    = ListPosition + [ 1  1 -2 -2] * ud.BorderWidth;

  HeadSpace = ( ListPosition(4) - 2*ud.BorderWidth - sum(TextPosition([1 2],4),1) ) / 3;

  pos(4) = TextPosition(2,4) + HeadSpace;

  set( ud.Handle.SectionText , 'units'    , 'pixels' , ...
                               'position' , pos            );

  pos(2) = pos(2) + pos(4);
  pos(4) = TextPosition(1,4) + HeadSpace;

  set( ud.Handle.ChapterText , 'units'    , 'pixels' , ...
                               'position' , pos            );


  %--------------------------------------------------
  % Fit HelpList between Logo and MessageList

  pos(1) = 1 + ud.Offset;
  pos(3) = figpos(3) - pos(1) + 1 - ud.Offset;
  pos(2) =  MsgPosition(2) + MsgPosition(4) + ud.Offset;
  pos(4) = LogoPosition(2) - ud.Offset - pos(2);

  set( ud.Handle.HelpFrame , 'units'    , 'pixels' , ...
                             'position' , pos            );
  
  set( ud.Handle.HelpList , 'units'    , 'pixels' , ...
                            'position' , pos + [ 1  1  -2 -2 ]*ud.BorderWidth  );

  %--------------------------------------------------
  % Fit MaximizePushButton

  pos    = get(ud.Handle.MessgFrame,'position');
  pos(4) = 2*ud.BorderWidth;

  set( ud.Handle.MessgPush , 'units'    , 'pixels' , ...
                             'position' , pos            );
 

%*******************************************************************
case 'MINIMIZE'

   mud = get(ud.Handle.MessgPush,'userdata');
 
   mud.Visible = 1 - mud.Visible;

   sets = { 'off'  'on' };    

   set( ud.Handle.MessgPush , 'userdata' , mud     , ...
                              'visible'  , sets{ 2-mud.Visible } );

   set( [ ud.Handle.MessgFrame ud.Handle.MessgList ]  , ...
                              'visible'  , sets{ 1+mud.Visible } );

   browser('Resize',fig);
 


%*******************************************************************
case 'ADD'

  % Add Contents
 
  if isempty(VarArg)
    Msg = [ Msg0 'Not enough InputArguments.' ];
    return
  end

  %-------------------------------------------------------------
  % Check for Directory

  if ( ( size(VarArg,1) == 1 )  &  ...
       ischar(VarArg{1})        &  ~isempty(VarArg{1})  &  ...
       ( prod(size(VarArg{1})) == size(VarArg{1},2) ) ) 

    d = dir(VarArg{1});

    if isempty(d)

      Msg = [ Msg0 'Invalid Input, single String must be a DirectoryName.' ];
  
      return

    end
 
      %------------------------------------------
      % Get ItemStructure

      [Msg,VarArg] = checkdir( ud.Handle.MessgFrame , VarArg{1} );

      if ~isempty(Msg)  &  isempty(VarArg)
        Msg = [ Msg0 'Invalid Input.' nl Msg ];
      end

      if isempty(VarArg)
        return
      end
  

  end
 
  %-------------------------------------------------------------
  % Check ItemStructure

  [Msg,V] = checkadd( ud.Handle.MessgFrame , VarArg , size(Msg0,2) , [] );
  
  if ~isempty(Msg)
    Msg = [ Msg0 'Invalid Inputs.' nl Msg ];
  end

  if isempty(V)
    return
  end


  %-------------------------------------------------------------
  % Add ItemStructure to Parent

  if strcmp( typ , 'figure' )
     parent = ud.Menu.Contents;
     root   = 0;
     nr     = zeros(1,0);
     hst    = zeros(1,0);
  else
     parent = mud.ID;
     root   = mud.Root + mud.ID * ( mud.Root == 0 );
     nr     = mud.Nr;
     hst    = mud.History;
  end


  out = cell(1,3);

  [ out{1} , out{2} , out{3} ] = makemenu(parent,root,fig,nr,hst,V);

   ud.Children.Handle = cat( 1 , ud.Children.Handle , out{1}(:) );
   ud.Children.Tag    = cat( 1 , ud.Children.Tag    , out{2}(:) );

  set( fig , 'userdata' , ud );
 

               n = min(Nout,size(out,2));

  varargout(1:n) = out(1:n);
    

%*******************************************************************
case 'ACTIVATE'

% Activates a Menu

   HC = ud.Menu.Contents;

   cud = get( HC , 'userdata' );

   %---------------------------------------------
   % Check Current Active Menu

   if isequal(mud.ID,cud.ActiveID) 

     % Selected Menu allready active

     % ... |  isequal(mud.ID,cud.ActiveRoot)  
     % Selected Root allready ActiveRoot
     
      out = { cud.ActiveID  cud.ActiveTag  get(ud.Handle.HelpList,'string') };

        n = min(Nout,size(out,2));

        varargout(1:n) = out(1:n);


      set( fig , 'visible' , 'on' )

      figure(fig)

      return

   end


   %---------------------------------------------
   % Check Current Active Chapter
   %  Activate new Chapter

   root = mud.Root + mud.ID * ( mud.Root == 0 );

   if  ~isequal( root , cud.ActiveRoot )   

      %------------------------------------------
      % DisActivate old Root

      if ~isempty( cud.ActiveRoot )
         browser( 'DisActivate' , fig );
      end

      cud = get( HC , 'userdata' );


      %------------------------------------------
      % Activate new Root

      rud = get( root , 'userdata' );

      if ~isempty( rud.Children )
      % Root contains Children  ==>
      %  Flip them to Top of ContentsMenu

        set( rud.ID , 'callback' , get(rud.Double,'callback') );

        set( rud.Double , 'visible' , 'off' );

        % Children of Root
        rch = cat( 1 , rud.Children , rud.Double  );

        set( rch , 'parent' , HC );

        % Other Chapter
        och = cud.Children( find( cud.Children ~= rud.ID ) );

        ch = cat( 1 , rud.ID , rud.Children  , och , rud.Double );

        set( HC , 'children' , ch(end:-1:1) );

        set( och(1:~isempty(och)) , 'separator' , 'on' );

      end

      set( ud.Handle.ChapterText , 'string' , rud.Titel );

   end


   %------------------------------------------------------------
   % Selected Tag

   cud.SelectedTag = mud.Tag;

   %------------------------------------------------------------
   % Root Selected, switch to 1. Children if Text and File empty

   if ( mud.Root == 0 )

      if isempty( mud.Text )  &  isempty( mud.File )  &  ...
         ~isempty(mud.Children)
         
         mud = get( mud.Children(1) , 'userdata' );

      end

   end 

   %------------------------------------------------------------

   cud.ActiveRoot = mud.Root + mud.ID * ( mud.Root == 0 );
   cud.ActiveID   = mud.ID;
   cud.ActiveTag  = mud.Tag;

   cud.Value   = cud.Value + ( size(cud.History,2) - cud.Value ) * ...
                             ( size(cud.History,2) < cud.Value );

   cud.History = cat( 2 , cud.History(1:cud.Value) , mud.ID );
   cud.Value   = cud.Value + 1;

   set( HC , 'userdata' , cud );

   sets = { 'off'  'on' };
   sets = sets{ 1 + ( cud.Value > 0 ) };

   set( ud.Menu.Back    , 'enable' , sets  );
   set( ud.Menu.Forward , 'enable' , 'off' );

   set( get(ud.Menu.Save,'children') , 'enable' , 'on' );

   %------------------------------------------------------------
   % Set all ContentsMenu's checked off

   set( ud.Children.Handle , 'foregroundcolor' , ud.Color.Default );

     if strcmp(computer,'PCWIN');
        set( ud.Children.Handle , 'checked' , 'off' );
     end

   %------------------------------------------------------------
   % Titel

   tt = mud.Titel ;

   % Check for Concatinate
   if strcmp( tt(1) , '_' )

     tt = tt( 2 : end );

     if mud.Root ~= 0
      
       pud = get( mud.Parent , 'userdata' );
       t0  = pud.Titel;
       t0  = t0( 1+strcmp(t0(1),'_') : end );
       tt  = cat( 2 , t0 , ', ' , tt );

     end

   end

   set( ud.Handle.SectionText , 'string' , tt );

   %------------------------------------------------------------
   % Check for MatlabEvent

    [MatlabOk,par] = getpar( ud.Matlab.Event , mud.Matlab );

    if MatlabOk

      try

        MatlabOk = event( ud.Handle.MessgFrame , ...
                          ud.Matlab.Event{MatlabOk} , par, ...
                          mud.Titel( 1+strcmp(mud.Titel(1),'_') : end ) , ud.MsgTag);

      catch
    
        MatlabOk = 0;

        Msg = [ Msg0  'Error using EVENT.' nl lasterr ]; 

        browser('MsgBox',fig,Msg,'Error','warn');

      end

    end


   %------------------------------------------------------------
   % Titel and Text

   % Mode == 1
   %  Don't use the ChapterNumber for Header
   %  Don't recurse trough Children if Chapter selected
   %

   pointer = get(fig,'pointer');
             set(fig,'pointer','watch');

   try

      txt = gettext( ud.Handle.MessgFrame , mud.ID , ...
                          cell(0,1) , 1 , ...
                          ud.Children.Handle , ud.Children.Tag , ...
                          ud.Character , ud.Color.Selected , ...
                          ud.Matlab , 0 );

   catch
    
     txt = { '' };
     ref = cell(0,2);

     Msg = [ Msg0  'Error using GETTEXT.' nl lasterr ]; 

     browser('MsgBox',fig,Msg,'Error','warn');

   end

   set(fig,'pointer',pointer);


     set( ud.Handle.HelpList , 'value'      ,  1  , ...
                               'listboxtop' ,  1  , ...
                               'string'     , txt        );

   %------------------------------------------------------------
   % All SubMenus are checked now, check the Parents

   set( mud.History , 'foregroundcolor' , ud.Color.Selected );

    if strcmp(computer,'PCWIN');
       set( mud.History , 'checked' , 'on' );
    end

   set(root,'checked','on');


   %------------------------------------------------------------

   out = { cud.ActiveID  cud.ActiveTag  txt };

     n = min(Nout,size(out,2));

   varargout(1:n) = out(1:n);


   %------------------------------------------------------------

   set( fig , 'visible' , 'on' );

   if 1   % ~( MatlabOk )

     figure(fig)

   end

   %------------------------------------------------------------



    
%*******************************************************************
case 'DISACTIVATE'

% DisActivates a Chapter

   HC = ud.Menu.Contents;

   cud = get( HC , 'userdata' );

   if isempty( cud.ActiveRoot )
      return
   end

  if strcmp(typ,'figure');
     root = cud.ActiveRoot;
  else
     root = mud.Root + mud.ID * ( mud.Root == 0 );
  end


  rud = get( root , 'userdata' );

  if ~isempty(rud.Children)
  % Root contains Children  ==>
  %  Flip them back from HC to Root

      % Children of Root
      rch = cat( 1 ,  rud.Double , rud.Children );
 
     set(  rch , 'parent'   , root );
     set( root , 'children' , rch(end:-1:1) , ...
                 'callback' , ''                  );

     set( rud.Double , 'visible' , 'on' );

  end

     set(   HC , 'children' , cud.Children(end:-1:1) );

     set( ud.Children.Handle , 'foregroundcolor' , ud.Color.Default );

     if strcmp(computer,'PCWIN');
        set( ud.Children.Handle , 'checked' , 'off' );
     end

     set( cud.Children , 'separator' , 'off' , ...
                         'checked'   , 'off'       );

     set( rud.Double , 'visible' , 'on' );

     set( ud.Handle.ChapterText , 'string'     , ...
     get( ud.Handle.ChapterText , 'userdata' )       );

     set( ud.Handle.SectionText , 'string'     , ...
     get( ud.Handle.SectionText , 'userdata' )       );

     set( ud.Handle.HelpList    , 'value'      ,  1  , ...
                                  'listboxtop' ,  1  , ...
                                  'string'     ,       ...
     get( ud.Handle.HelpList    , 'userdata' )             );


  out = { cud.ActiveID  cud.ActiveTag };

  cud.ActiveRoot = [];
  cud.ActiveID   = [];
  cud.ActiveTag  = '';

  set( HC , 'userdata' , cud );

  %------------------------------------------
  % Disable Save-Chapter and Save-Section

  ch = get(ud.Menu.Save,'children');
  set( ch([1 2]) , 'enable' , 'off' );

  %------------------------------------------

   n = min(Nout,size(out,2));

   varargout(1:n) = out(1:n);


%*******************************************************************
case 'BACK'

   HC = ud.Menu.Contents;

   cud = get( HC , 'userdata' );

   if cud.Value <= 0

      set( ud.Menu.Back , 'enable' , 'off' );

     out = { cud.ActiveID  cud.ActiveTag };
       n = min(Nout,size(out,2));

      varargout(1:n) = out(1:n);

      return
   end

   %------------------------------------------
   % Store current History

   hist = cud.History;
   val  = cud.Value - 1;

   %------------------------------------------
   % Activate previous Menu

   if val == 0
     browser( 'DisActivate' , fig    );
   else
     browser( 'Activate' , hist(val) );
   end

   %------------------------------------------
   % Reset current History

   cud = get( HC , 'userdata' );

   cud.History = hist;   
   cud.Value   = val;

   set( HC , 'userdata' , cud );

   sets = { 'off'  'on' };
   sets = sets{ 1 + ( val > 0 ) };

   set( ud.Menu.Back    , 'enable' , sets  );
   set( ud.Menu.Forward , 'enable' , 'on'  );
   

     out = { cud.ActiveID  cud.ActiveTag };

       n = min(Nout,size(out,2));

      varargout(1:n) = out(1:n);




%*******************************************************************
case 'FORWARD'

   HC = ud.Menu.Contents;

   cud = get( HC , 'userdata' );

   if cud.Value >= size(cud.History,2)
      set( ud.Menu.Forward , 'enable' , 'off' );

     out = { cud.ActiveID  cud.ActiveTag };
       n = min(Nout,size(out,2));

      varargout(1:n) = out(1:n);

      return

   end

   %------------------------------------------
   % Store current History

   hist = cud.History;
   val  = cud.Value + 1;

   %------------------------------------------
   % Activate following Menu

   browser('Activate',hist(val));

   %------------------------------------------
   % Reset current History

   cud = get( HC , 'userdata' );

   cud.History = hist;   
   cud.Value   = val;

   set( HC , 'userdata' , cud );

   sets = { 'off'  'on' };
   sets = sets{ 1 + ( val < size(hist,2) ) };

   set( ud.Menu.Forward , 'enable' , sets  );
   set( ud.Menu.Back    , 'enable' , 'on'  );


     out = { cud.ActiveID  cud.ActiveTag };

       n = min(Nout,size(out,2));

      varargout(1:n) = out(1:n);



%*******************************************************************
case 'MENU'

% Adds a ContentsMenu to another Figure

   if isempty(VarArg)
     Msg = [ Msg0 'Not enough InputArguments.' ];
     return
   end

  %---------------------------------------------------
  % Check FigureHandle

   parent = VarArg{1};

   VarArg = VarArg(2:end);

   ok = ( isnumeric(parent) & ( prod(size(parent)) == 1 ) );

   if ok
     ok = ishandle(parent);
     if ok
        ok = any( strcmp( get(parent,'type') , { 'figure' 'uimenu' } ) );
     end
   end

   if ~ok
     Msg = [ Msg0 'Parent must be a Handle for a Figure or UIMenu.' ];
     return
   end

  %---------------------------------------------------
  % Check for Label

   label = '';

   if ~isempty(VarArg)

      label = VarArg{1};
         ok = ( ischar(label)  &  ~isempty(label)  &  ...
                ( prod(size(label)) == size(label,2) ) );
     if ~ok
         Msg = [ Msg0  'Input Label must be a nonempty String.' ];
         return
     end

   end

  %---------------------------------------------------
  % Get current State and DisActivate

    HC  = ud.Menu.Contents;
    cud = get( HC , 'userdata' );

    browser( 'DisActivate' , fig );

  %---------------------------------------------------
  % Copy Menu

   if strcmp( typ , 'figure' );
   % Copy of full ContentsMenu
      id = HC;
   else
      id = mud.ID;
   end

   hm = copyobj( id , parent );

   if ~isempty(label)
      set( hm , 'label' , label );
   end

   NewMenu = hm;
   if ~isempty(get(hm,'children'))
     NewMenu = get(hm,'children');
   end

   ud.Children.Menu = cat( 1 , ud.Children.Menu, NewMenu );

   set( fig , 'userdata' , ud );

  %---------------------------------------------------
  % Reset State

   hist = cud.History;
   val  = cud.Value;

   if val > 0  
     browser( 'Activate' , hist(val) );
   end

   cud = get( HC , 'userdata' );

   cud.History = hist;
   cud.Value   = val;

   set( HC , 'userdata' , cud );


   sets = { 'off'  'on' };

   set( ud.Menu.Back    , 'enable' , sets{ 1 + ( val > 0            ) } );
   set( ud.Menu.Forward , 'enable' , sets{ 1 + ( val < size(hist,2) ) }  );


%*******************************************************************
case 'RESET'


  browser( 'DisActivate' , fig );

  % Delete all MessageBoxes
  % Delete all Children of ContentsMenu in other Figures
 
  browser( 'Delete' , fig );


 e0 = zeros(0,0); 
 e1 = zeros(0,1);        % Empty matrix: 0-by-1
 e2 = zeros(1,0);        % Empty matrix: 1-by-0

 ec =  cell(0,1);        % Empty cell array: 0-by-1

 es = char(zeros(1,0));  % Empty string: 1-by-0


 HC  = ud.Menu.Contents;
 cud = get( HC , 'userdata' );

 cud.Children   = e1;
 cud.ActiveRoot = e0;
 cud.ActiveID   = e0;
 cud.ActiveTag  = es;
 cud.History    = e2;
 cud.Value      =  0;

 ud.Children.Handle = e1;
 ud.Children.Tag    = ec;
 ud.Children.Menu   = e1;


  shh = get(0,'showhiddenhandles');
        set(0,'showhiddenhandles','on')

  delete( get(HC,'children') );

        set(0,'showhiddenhandles',shh)


   set( ud.Menu.Back    , 'enable' , 'off' );
   set( ud.Menu.Forward , 'enable' , 'off' );

  set( HC  , 'userdata' , cud );

  set( fig , 'userdata' ,  ud );

%*******************************************************************
case 'TEXT'

      out = cell(1,3);

    color = [];  % Don't set the Item to SelectedColor

    % Mode == 0
    %  Start with ChapterNumber
    %   Recurse


  if strcmp( typ , 'figure' )
    cud = get(  ud.Menu.Contents , 'userdata' );
    mud = get( cud.Children      , 'userdata' );
  end
  
  out{2} = struct( 'Nr'       , { [] } , ...
                   'Titel'    , { '' } , ...
                   'Text'     , { '' } , ...
                   'Note'     , { '' } , ...
                   'Reference', { cell(0,2) } , ...
                   'Image'    , { '' } , ...
                   'Contents' , { cell(0,2) } , ...
                   'Children' , { [] }           );

  out{2}.Children = out{2}(zeros(0,1));

  
  if iscell(mud)
     mud = cat(1,mud{:});
  end

  n = size(mud,1);

  out{1} = cell(n,1);
  out{2} = out{2}(ones(n,1));

  out{3} = cell(n,2);    % { Nr  Titel }

 try

   for ii = 1 : n

     [ out{1}{ii} , out{2}(ii) ] = ...
            gettext( ud.Handle.MessgFrame , mud(ii).ID , ...
                     cell(0,1) , 0 , ...
                     ud.Children.Handle , ud.Children.Tag , ...
                     ud.Character , color , ...
                     ud.Matlab , 0 );
   end

   
 catch
    
     Msg = [ Msg0  'Error using GETTEXT.' nl lasterr ]; 

 end
 

 if isempty(Msg)

   out{1} = cat(1,out{1}{:});

   MsgC = '';

   try
 
     [out{2},out{3},cont] = getcont(out{2});
 
   catch

     MsgC = lasterr;

   end

   if isempty(MsgC)

     out{1} = cat( 1 , cont , { '' } , out{1} );
     
   else

     Msg = [ Msg0  'Error using GETCONT.' nl MsgC ]; 
 
   end

 end


 if ( Nout >= 3 )  &  isempty(Msg) 

     out{4} = cell( prod(size(out{1})) , 2 );
     out{4}(:,1) = out{1}(:);
     out{4}(:,2) = { nl };
     out{4} = permute( out{4} , [ 2  1 ] );
     out{4} = cat( 2 , out{4}{:} );

 end


   n = min(Nout,size(out,2));

   varargout(1:n) = out(1:n);


%*******************************************************************
case 'SAVE'

   IsMsg = ( nargout == 0 );

   if isempty(VarArg)
      mode = 'All';
   else
      mode = VarArg{1};
   end

   if ~( ischar(mode)  &  ~isempty(mode)  &  ...
         ( prod(size(mode)) == size(mode,2) ) );

     Msg = [ Msg0  'Mode must be a String.' ];

     if IsMsg
        browser( 'MsgBox' , fig , Msg , 'Error' , 'warn' );
     end

     return

   end

   %------------------------------------------------
   if strcmp( mode , 'All' )

      TextIn = fig;

   %------------------------------------------------
   else

      HC = ud.Menu.Contents;

      cud = get( HC , 'userdata' );

      switch upper(mode)

        %----------------------------
        case 'SECTION'

          TextIn = cud.ActiveID;
 
        %----------------------------
        case 'CHAPTER'

          TextIn = cud.ActiveRoot;

        %----------------------------
        otherwise

          Msg = [ Msg0  'Mode must be any of ''Section'', ''Chapter'' or ''All''.' ];
 
          if IsMsg
             browser( 'MsgBox' ,fig ,  Msg , 'Error' , 'warn' );
          end

          return

      end

      %----------------------------
      if isempty(TextIn)

         Msg = [ Msg0  'No '  mode  ' activated.' ];
  
         if IsMsg
            browser( 'MsgBox' , fig , Msg , 'Error' , 'warn' );
         end
         
         return

      end

   end

   %------------------------------------------------
   % Get File

   name = cat( 2 , 'Save Text to ... ' );

   [file,pfad]= uiputfile( '*.txt' , name );

   if isequal(pfad,0) | isempty(file)

      return

   end

   file = cat( 2 , pfad , file );

   if ( exist(file,'dir') == 7 )

      Msg = [ Msg0  'Selected Filename is a Directory.' ];

      if IsMsg
         browser( 'MsgBox' , fig , Msg , 'Error' , 'warn' );
      end

      return
 
   end

   %------------------------------------------------
   % Open File for Check

   HMsg = ud.Handle.MessgFrame;

   msg_list( HMsg , 'Message' , ['Open File for check   '  file ] , 'new' );

   fid = fopen(file,'wt');

   if ( fid == -1 )

      msg_list( HMsg , 'Message' , 'Error' , 'append' );

      Msg = [ Msg0  'Can''t open File '  file  ' for writing.' ];

      if IsMsg
         browser( 'MsgBox' , fig , Msg , 'Error' , 'warn' );
      end

      return

   end

   fclose(fid);

   msg_list( HMsg , 'Message' , 'ok' , 'append' );

   %------------------------------------------------

   pointer = get(fig,'pointer');
             set(fig,'pointer','watch');

   %------------------------------------------------
   % Get Text

   try

     [ Msg , c , v , h , txt ] = browser( 'Text' , TextIn );

   catch
     
       Msg = lasterr;

   end

   if ~isempty(Msg)

      Msg = [ Msg0 'Error call BROWSER( Text ).' nl Msg ];

      if IsMsg
         browser( 'MsgBox' , fig , Msg , 'Error' , 'warn' );
      end

      set(fig,'pointer',pointer);

      return

   end 

   %------------------------------------------------
   % Open File for Check

   HMsg = ud.Handle.MessgFrame;

   msg_list( HMsg , 'Message' , ['Open File for writing   '  file ] , 'new' );

   fid = fopen(file,'wt');

   if ( fid == -1 )

      msg_list( HMsg , 'Message' , 'Error' , 'append' );

      Msg = [ Msg0  'Can''t open File '  file  ' for writing.' ];

      if IsMsg
         browser( 'MsgBox' , fig , Msg , 'Error' , 'warn' );
      end

      set(fig,'pointer',pointer);

      return

   end

   msg_list( HMsg , 'Message' , 'ok' , 'append' );

   %------------------------------------------------

   out = { file txt };

   %------------------------------------------------
   % Write Text

   txt = strrep( txt , '%' , '%%' );
   txt = strrep( txt , '\' , '\\' );

   msg_list( HMsg , 'Message' , ['Write Text into File   '  file ] , 'new' );

   Msg = '';

   try
      fprintf( fid , txt );
      fclose( fid );
   catch
      Msg = lasterr;
   end

   set(fig,'pointer',pointer);

   if ~isempty(Msg)

      msg_list( HMsg , 'Message' , 'Error' , 'append' );
      
      Msg = [ Msg0 'Error write Text into File '  nl file nl nl Msg ];

      if IsMsg
         browser( 'MsgBox' , fig , Msg , 'Error' , 'warn' );
      end

      return

   end

   msg_list( HMsg , 'Message' , 'ok' , 'append' );

   %------------------------------------------------

   n = min(Nout,size(out,2));

   varargout(1:n) = out(1:n);


%*******************************************************************
case 'REFRESH'


  ch = get( ud.Menu.Contents , 'Children' );

       set( ud.Menu.Contents , 'Children' , ch );


%*******************************************************************
case 'PROBLEM'

  txt = [ 'On some Systems the ContentsMenu''s are not visible' ...
           nl    ...
          ' after Activating a new Chapter!'  ...
           nl nl ... 
          'Use the <Refresh>-Menu to have a correct Visibility'  ...
           nl  ...
          ' of the ContentsMenu''s.'    nl         ];
 
 try

   [msg,c,cmap] = readxpm('menu_bug.xpm');

 catch

   msg = lasterr;

 end


 if ~isempty(msg)

   browser('MsgBox',fig,txt,'Sorry','warn');

   return

 end

 offs = 10;  % PixelOffset

 s1 = size(c,1);
 s2 = size(c,2);

 xl = [ 1  s2 ] + 0.5*[ -1  1 ];
 yl = [ 1  s1 ] + 0.5*[ -1  1 ];

 MsgFig = figure( 'menubar'     , 'none'      , ...
               'numbertitle' , 'off'       , ...
               'color'       , [ 1  1  1 ] , ...
               'colormap'    , cmap        , ...
               'name'        , 'Sorry'     , ...
               'visible'     , 'off'       , ...
               'nextplot'    , 'add'       , ...
               'tag'         , ud.MsgTag   , ... 
               'handlevisibility' , 'callback' , ...
               'createfcn'        , ''         , ...
               'resize'           , 'off'      , ...
               'closerequestfcn'  , 'closereq'     );


  hc  = uicontrol( 'parent' , MsgFig , ...
                   'units'  , 'pixels' , ...
                   'position' , [ 1  1  offs offs ] , ...
                   'style'    , 'pushbutton'        , ... 
                   'string'   , 'Ok'                , ...
                   'ForeGroundColor' , [ 0  0  0 ]  , ...
                   'BackGroundColor' , get(ud.Handle.LogoFrame,'backgroundcolor') , ...
                   'FontUnits'       , get(ud.Handle.HelpList,'FontUnits' ) , ...
                   'FontSize'        , get(ud.Handle.HelpList,'FontSize'  ) , ...
                   'FontWeight'      , 'normal'        , ...
                   'horizontalalignment' , 'center'    , ...
                   'visible'             , 'on'        , ...
                   'handlevisibility'    , 'on'        , ...
                   'callback'            , 'delete(get(gcbo,''parent''));'  );


 axe = axes( 'parent'   , MsgFig          , ...
             'units'    , 'normalized'    , ...
             'position' , [ 0  0  1  1 ]  , ...
             'color'    , 'none'  , ...
             'xlim'     , xl      , ...
             'ylim'     , yl      , ...
             'xtick'    , []      , ...
             'ytick'    , []      , ...
             'xdir'     , 'normal'  , ...
             'ydir'     , 'reverse' , ...
             'visible'  , 'off'   , ...
             'nextplot' , 'add'   , ...
             'dataaspectratio'  , [ 1  1  1 ] , ...
             'handlevisibility' , 'callback'      );

 ht = text( 'parent' , axe , ...
            'units'  , 'pixels' , ...
            'position' , [ 1  1  0 ] , ...
            'string'   , 'Ok'        , ...
            'color'    , [ 0  0  0 ] , ...
            'FontUnits'       , get(ud.Handle.HelpList,'FontUnits' ) , ...
            'FontSize'        , get(ud.Handle.HelpList,'FontSize'  ) , ...
            'FontWeight'      , 'normal'        , ...
            'horizontalalignment' , 'left'      , ...
            'verticalalignment'   , 'bottom'       , ...
            'clipping' , 'off'    , ...
            'visible'  , 'on'     , ...
            'interpreter'         , 'none'          );


  ext_ok  = get( ht , 'extent' );

  set( ht , 'string' , txt );

  ext_txt = get( ht , 'extent' );

 
  uicpos         = zeros(1,4);
  uicpos(2)      = offs + 1;
  uicpos([3 4])  = [ 3  2 ] .* ext_ok([3 4]);


  axe_c     =  max( s1 , ext_txt(4) ) / 2;  % Half AxeHeight

  axepos    = zeros(1,4);
  axepos(1) = offs + 1;
  axepos(2) = uicpos(2) + uicpos(4) + offs + floor( axe_c - s1/2 ) ;
  axepos(3) = s2; 
  axepos(4) = s1; 
  
  txtpos    = zeros(1,3);
  txtpos(1) = axepos(3) + offs + 1;
  txtpos(2) = ceil( axepos(4)/2 - ext_txt(4)/2 - 1/3*ext_ok(4) );
             
  figpos    = zeros(1,4);
  figpos(3) = axepos(1) + axepos(3) + offs + ext_txt(3) + offs;
  figpos(4) = axepos(2) + axepos(4) + 2*offs; 
  
  ppos      = get(fig,'position');

  figpos(1) = ppos(1) + ( ppos(3) - figpos(3) ) / 2;
  figpos(2) = ppos(2) + ( ppos(4) - figpos(4) ) / 2;

  figpos([1 2]) = figpos([1 2]) + ( 50 - figpos([1 2]) ) .* ...
                                  ( 50 > figpos([1 2]) );
 
  uicpos(1) = floor( ( figpos(3) - uicpos(3) ) / 2 );

  
  set( MsgFig , 'units' , 'pixels' , 'position' , figpos );
  set( axe    , 'units' , 'pixels' , 'position' , axepos );
  set( ht     , 'units' , 'pixels' , 'position' , txtpos );
  set( hc     , 'units' , 'pixels' , 'position' , uicpos );

  hi = image( 'parent' , axe , ...
              'xdata' , [ 1  s2 ] , ...
              'ydata' , [ 1  s1 ] , ...
              'cdata' , c         , ...
              'cdatamapping' , 'direct' , ...
              'visible'      , 'on'           );

  set( MsgFig , 'visible' , 'on' );



%*******************************************************************
case 'CLOSE'

  % Delete all MessageBoxes

  delmsgbox(ud.MsgTag);

  set( fig , 'visible' , 'off' );


%*******************************************************************
case 'DELETE'

  % Delete all MessageBoxes

  delmsgbox(ud.MsgTag);


  % Delete all Children of ContentsMenu in other Figures
  
    hm = ud.Children.Menu;

    hm = hm( find( ishandle(hm) ) );
    hm = hm( find( strcmp( get(hm,'type') , 'uimenu' ) ) );

    if ~isempty(hm)

      delete(hm);
 
    end

%*******************************************************************
case 'MESSAGE'
  
  Msg = msg_list( ud.Handle.MessgFrame , 'Message' , VarArg{:} );

  if ~isempty(Msg)
    Msg = [ Msg0 'Error call MSG_LIST.' nl Msg ];
    browser('MsgBox',fig,Msg,'Error','warn');
  end



%*******************************************************************
case 'MSGBOX'

  MsgFig = msgbox(VarArg{:});

  set( MsgFig , 'resize' , 'on' , ...
                'tag'    , ud.MsgTag )

  drawnow


end
% case


%*******************************************************************
% SubFunctions
%*******************************************************************


%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [ok,typ] = checkin(h,types)

% CHECKIN  Checks 2. Input for BrowserFigure or BrowserMenu
%


  typ = '';

  ok = ( isnumeric(h) & ( prod(size(h)) == 1 ) );

  if ok

     ok = ishandle(h);

     if ok

       typ = get(h,'type');

       ok = any( strcmp( typ , types ) );

       if ok
         ud = get(h,'userdata');  
         ok = isstruct(ud);
         if ok
           fields = fieldnames(ud);
           ok = ( any(strcmp(fields,'ID' ))  &  ...
                  any(strcmp(fields,'Tag'))          );
           if ok
             ok = ( isequal( getfield(ud,'ID' ) ,      h         )  &  ...
                    isequal( getfield(ud,'Tag') ,  get(h,'tag')  )         );
             if ok  &  strcmp( typ , 'uimenu' )
             % BrowserMenu  ==>  Check for BrowserFigure
               ok = any( strcmp( fields , 'Figure' ) );
               if ok
                  ok = checkin( ud.Figure , 'figure' );
               end
             end 
           end
         end
       end
     end
  end


%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [Msg,V] = checkdir(HMsg,pfad)

% CHECKDIR  Checks Input to add Directory to Contents
%            
%   works recursivly, returns ItemStructure
%

Msg = '';

nl  = char(10);

fs  = filesep;

form = '%3.3d.';   % NumberFormat
sep  = ':';        % Separator 

%-----------------------------------------------
% DefaultStructure

e2 = ones(1,0);
ec = char(e2);  %  Empty string: 1-by-0

V0  = struct( 'Tag'    , { ec } , ...
              'Label'  , { ec } , ...
              'Titel'  , { ec } , ...
              'Text'   , { cell(0,1) } , ...
              'File'   , { ec } , ...
              'Note'   , { ec } , ...
              'Matlab' , { ec } , ...
              'Children' , { [] }       );

V0.Children = V0(e2);

V = V0(e2);


%-----------------------------------------------
% Read Directory

pfad = cat( 2 , pfad , fs( 1 : (end*(~strcmp(pfad(end),fs))) ) );

d = dir(pfad);
txt = [ 'CHECKDIR, dir( '  pfad  ' )' ];

msg_list(HMsg,'Message',txt,'new');

d = dir(pfad);

if isempty(d)

  Msg = 'Can''t read Directory.';

  msg_list(HMsg,'Message',Msg,'append');

  Msg = [ Msg nl pfad ];

  return

end


%-----------------------------------------------
% Look for BROWSER.INI

name = cellstr(str2mat(d.name));

if ~any( strcmp( upper(name) , 'BROWSER.INI' ) )
  Msg = [ 'No  BROWSER.INI  in Directory:' nl pfad ];
  return  
end

ii = find( strcmp( upper(name) , 'BROWSER.INI' ) );

file = cat( 2 , pfad , d(ii).name );



%-----------------------------------------------
% Read  BROWSER.INI

txt = [ 'CHECKDIR, ReadFile( '  file ' ) ' ];

msg_list(HMsg,'Message',txt,'new');

[Msg,bb] = readfile(file,'');

if ~isempty(Msg)

  msg_list(HMsg,'Message',Msg,'append');

  Msg = [ 'Error read  BROWSER.INI:' file nl Msg ];

  return

end

%-----------------------------------------------
% Get Definitions 

bb = rmblank(char(bb),2);

s1 = size(bb,1);
s2 = size(bb,2);

one1 = ones(s1,1);

sp = size(sep,2);
  
% Counter: 000 001 002 ...
cc = 0;

vc = e2;  % CounterVector for Items

ok = 1;

field = fieldnames(V0);
field = field(:)';  

while  ok  |  ( cc == 1 )
  
  ok  = 0;
  okV = 0;

  V1 = V0;

  Nr = sprintf(form,cc);

  for ff = field;  

    % String of ItemField
    if strcmp(ff{1},'Children')
       def = cat( 2 , Nr , 'Directory' );
    else
       def = cat( 2 , Nr , ff{1} );
    end

    sd = size(def,2);

    if s2  > sd+sp 
        
      % STRMATCH, Lines matching ItemField

      ss =  sum( ( double(bb(:,1:sd)) == one1*double(def) ) , 2 );


      if any( ss == sd )

         ok = ( ok  |  1 );

         ii = find( ss == sd );
         ii = ii(end);           % Use the last Line which found

         jj = findstr( bb( ii , sd+1 : s2 ) , sep );

        if ~isempty(jj)

          % Value of the Definition
          str = rmblank( bb( ii , sd+jj(1)+sp : end ) , 2 );

          if ~isempty(str)
    
             okV = ( okV  |  1 );

             switch ff{1}

              %----------------------------------------------
              case 'Children'

                str0 = str;
                str1 = cat( 2 , pfad , str );

                ok0  = ~isempty(dir(str0));
                ok1  = ~isempty(dir(str1));

                if ok0 & ~ok1
                   str1 = str0;
                end
                            
                [MsgV,V1] = checkdir( HMsg , str1 );

                Msg = cat( 2 , Msg , nl(1:(end*(~isempty(Msg)))) , MsgV );
 
              %----------------------------------------------
              case 'File'

                  [p,n,e] = fileparts(str);

                   name = cat(2,n,e);  % Only FileName
 
                   str1 = cat(2,pfad,str);
                   str2 = cat(2,pfad,name);

                   str3 = which(name);

                   ok1 = ( exist(str1) == 2 );
                   ok2 = ( exist(str2) == 2 );

                   ok3 = ~isempty(str3); % File in  Matlab's SearchPath

                   ok13 = isequal( str1 , str3 );
                   ok23 = isequal( str2 , str3 );
  
                   if  ( ok13 | ok23 )           |  ...
                       ( ok3  & ~( ok1 | ok2 ) )
                   % File in Matlab's SearchPath

                      str = name;

                   elseif ok1

                      str = str1;

                   elseif ok2

                      str = str2;

                   end


                 V1 = setfield( V1 , ff{1} , str );
  
              %----------------------------------------------
              otherwise

                 V1 = setfield( V1 , ff{1} , str );

             end
             % switch

          end
          % ~isempty(str)

       end
       % ~isempty(jj)
      end                   
      % STRMATCH
    end
    % s2  > sd+sp
  end
  % ff

  if okV  &  ~isempty(V1)
    vc                 = cat( 2 , vc , cc ); 
    V(prod(size(V))+1) = V1;
  end

  cc = cc + 1;

end
% while


if ~isempty(V)

  if ( vc(1) == 0 ) 

     V(1).Children = V(2:end);

     V  = V(1);

  end

end


%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [Msg,V] = checkadd(HMsg,V0,nbl,nr)

% CHECKADD  Checks Input to add to Contents
%
%   works recursivly
%

Msg = '';

nl  = char(10);

if nargin < 5
  nr  = [];
end

blank = char( 32 * ones(1,nbl) );

if isempty(nr)
  txt = '';
else
  txt = [ ' for '  sprintf('%.0f.',nr)  ' ' ];
end

Msg0 = [ blank 'Invalid Structure'  txt  ': '  ];

 nl0 = [ nl char(32*ones(1,size(Msg0,2))) ];

%-----------------------------------------------
% DefaultStructure

e1 = ones(0,1);
e2 = ones(1,0);
ec = char(e2);  %  Empty string: 1-by-0

V1  = struct( 'Tag'    , { ec } , ...
              'Label'  , { ec } , ...
              'Titel'  , { ec } , ...
              'Text'   , { cell(0,1) } , ...
              'File'   , { ec } , ...
              'Note'   , { ec } , ...
              'Matlab' , { ec } , ...
              'Children' , { [] } , ...
              'History'  , { e2 }       );

V1.Children = V1(e1);

V = V1(e1);


%-----------------------------------------------
% Basic Check of Input

[Msg,V0] = structchk(V0);

if ~isempty(Msg)
  Msg = [ Msg0  Msg ];
  return
end


V0 = V0(:);


%-----------------------------------------------
% Check for allowed FieldNames

f0 = fieldnames(V0);
f1 = fieldnames(V1);

n0 = size(f0,1);
n1 = size(f1,1);

ok0 = zeros(n0,1);

ok1      = cell(2,n1);

ok1(1,:) = f1(:)';
ok1(2,:) = {0};
ok1      = struct(ok1{:});


for ii = 1 : n1
  is_field = strcmp( f1{ii} , f0 );
  if any( is_field )
    ok1 = setfield( ok1 , f1{ii} , 1 );
    ok0( find(  is_field ) ) = 1;
  end
end


%-----------------------------------------------
% Check for required FieldNames

ok = ( ok1.Tag  |  ok1.Label  |  ok1.Titel );
if ~ok
  Msg = [ Msg nl0(1:(end*(~isempty(Msg))))  ...
          'FieldNames ''Tag'', ''Label'' or ''Titel'' required.' ];
end

ok = ( ok1.Text  |  ok1.File  |  ok1.Children );
if ~ok
  Msg = [ Msg nl0(1:(end*(~isempty(Msg))))  ...
          'FieldNames ''Text'', ''File'' or ''Children'' required.' ];
end


if ~isempty(Msg)
  Msg = [ Msg0  Msg ];
  return
end


%-----------------------------------------------
% Fill up Structure, Check new Values

  n0 = size(V0,1);
 
   V = V1(ones(n0,1));

 okV = ones(n0,1);

MsgC = '';   % ErrorMessages for Children, to append later


for ii = 1 : n0

  MsgV = '';   % ErrorMessage

  NrTxt = sprintf('%.0f.',cat(2,nr,ii));

  Msg0 = [ blank 'Invalid Values for ' NrTxt ' : ' ];

   nl0 = [ nl char(32*ones(1,size(Msg0,2))) ];

  %----------------------------------------------------------------
  % Fill up V1 with Values from V0

  ok2      = cell(2,n1);

  ok2(1,:) = f1(:)';
  ok2(2,:) = { 1 };
  ok2      = struct(ok2{:});


  for jj = 1 : n1

    if getfield(ok1,f1{jj});

      val = getfield( V0 , {ii} , f1{jj} );

      if ~isempty(val)

        ok  = 0;

        switch f1{jj}

         %----------------------------------------------------------------
         case { 'Tag'  'Label'  'Titel'  'Matlab'  'File' }
           ok = ( ischar(val)  &  ( prod(size(val)) == size(val,2) ) );
           if ~ok
             MsgV = [ MsgV nl0(1:(end*(~isempty(MsgV)))) ...
                      'Values for '''  f1{jj}  ''' must be a String.' ];
           end  

         %----------------------------------------------------------------
         case { 'Text'   'Note' }

           ok = ( ischar(val) | iscellstr(val) );

           if ok

             if iscellstr(val)
               val = char(val(:));
             end

           else

               MsgV = [ MsgV nl0(1:(end*(~isempty(MsgV)))) ...
                        'Values for '''  f1{jj}  ...
                        ''' must be a CharArray or CellStringArray.' ];

           end

         %----------------------------------------------------------------
         case { 'Children' }

           [msg,val] = checkadd(HMsg,val,nbl,cat(2,nr,ii));

           ok = ( isempty(msg) | ~isempty(val) );
           if ~isempty(msg)
             MsgC = [ MsgC nl(1:(end*(~isempty(MsgC)))) msg ];
           end           

        end
        % switch

        if ok
      
          V = setfield( V , {ii} , f1{jj} , val );

        end

        ok2 = setfield( ok2 , f1{jj} , ok );


      elseif ~strcmp( f1{jj} , 'Children' )


        val = ec;  %  Empty string: 1-by-0
 

      end
      % ~isempty(val) 

    end
    % ok1

  end
  % jj

  %------------------------------------------------------
  % Check Tag and Title for Concatinate

  bad = 1;

  ct = struct( 'Tag'   , { 0 } , ...
               'Titel' , { 0 }       );

  fct = fieldnames(ct);

  for ff = fct(:)

     val = getfield( V(ii) , ff{1} );

     if getfield( ok2 , ff{1} )
     % Value ok, check for '_'
        if ~isempty(val)

           cc = strcmp(val(1),'_') * ( 1 + ( size(val,2) == 1 ) );

           ct = setfield( ct , ff{1} , cc );

        end
     end

  end
         
  %------------------------------------------------------
  % Check for required Values
           
  if ( ( ( isempty(V(ii).Tag  ) | ( ct.Tag == 2 ) ) & ok2.Tag   )  &  ...
       (   isempty(V(ii).Label)                     & ok2.Label )  &  ...
       ( ( isempty(V(ii).Titel) | ( ct.Tag == 2 ) ) & ok2.Titel )         ) 

     MsgV = [ MsgV nl0(1:(end*(~isempty(MsgV)))) ...
              'Any of ''Tag'', ''Label''  or  ''Titel'' must be nonempty.' ];

  else

    if isempty(V(ii).Titel) | ( ct.Titel == 2 )
      if isempty(V(ii).Label)
         V(ii).Titel = cat( 2 , V(ii).Titel , V(ii).Tag );
      else
         V(ii).Titel = cat( 2 , V(ii).Titel , V(ii).Label );
      end 
    end

    if isempty(V(ii).Label)
      if isempty(V(ii).Tag) | ( ct.Tag == 2 )
         V(ii).Label = V(ii).Titel;
      else
         V(ii).Label = V(ii).Tag;
      end 
    end

    if isempty(V(ii).Tag)  |  ( ct.Tag == 2 )
      if isempty(V(ii).Label)
         V(ii).Tag = cat( 2 , V(ii).Tag , V(ii).Titel );
      else
         V(ii).Tag = cat( 2 , V(ii).Tag , V(ii).Label );
      end 
    end

  end

  if ( ( isempty(V(ii).Text    ) & ok2.Text     )   &  ...
       ( isempty(V(ii).File    ) & ok2.File     )   &  ...
       ( isempty(V(ii).Matlab  ) & ok2.Matlab   )   &  ...
       ( isempty(V(ii).Children) & ok2.Children )          ) 

     MsgV = [ MsgV nl0(1:(end*(~isempty(MsgV)))) ...
              'Any of ''Text'', ''File'', ''Matlab''  or ' ...
              '''Children'' must be nonempty.' ];

  end

  %------------------------------------------------------
  % Append ErrorMessage

  okV(ii) = isempty(MsgV);
 
  if ~okV(ii)

    Msg = [ Msg nl(1:(end*(~isempty(Msg)))) Msg0 MsgV ];

  end
     
end
% ii


%------------------------------------------------------
% Append ChildrenMessage

Msg = [ Msg nl(1:(end*(~isempty(Msg))))  MsgC ];


V = V(find(okV));




%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [ Hdl , Tag , V1 ] = makemenu(parent,root,fig,nr,hst,V);

% MAKEMENU   Create Menus from Structure 
%              works recursivly
%
% Returns new Structure and All Handles and Tags
%

V = V(:);
n = size(V,1);

zz = zeros(0,1);

Hdl = zz;         % created Handles
Tag = cell(0,1);  %   Tags of created Handles
 

V1 = struct( 'ID'       , { []     } , ...
             'Parent'   , { parent } , ...
             'Root'     , { root   } , ...
             'Figure'   , { fig    } , ...
             'Nr'       , { nr } , ...
             'Tag'      , { '' } , ...
             'Double'   , { zz } , ...
             'Label'    , { '' } , ...
             'Titel'    , { '' } , ...
             'Text'     , { '' } , ...
             'File'     , { '' } , ...
             'Note'     , { '' } , ...
             'Matlab'   , { '' } , ...
             'Children' , { zz } , ...
             'History'  , { [] }               );

V1 = V1(ones(n,1));


pud = get( parent , 'userdata' );


for ii = 1 : n

  if strcmp( V(ii).Tag(1) , '_' )
    V(ii).Tag = V(ii).Tag(2:end);
    if root ~= 0
      V(ii).Tag = cat( 2 , pud.Tag , V(ii).Tag );
    end
  end

  V1(ii).Nr     = cat(2,nr,ii);
  V1(ii).Tag    = V(ii).Tag;
  V1(ii).Label  = V(ii).Label;
  V1(ii).Titel  = V(ii).Titel;
  V1(ii).Text   = V(ii).Text;
  V1(ii).File   = V(ii).File;
  V1(ii).Note   = V(ii).Note;
  V1(ii).Matlab = V(ii).Matlab;

  if root == 0
    label = V1(ii).Label;
  else
    label = [ '  '  V1(ii).Label ];
  end
  
  V1(ii).ID = uimenu('parent'    , parent       , ... 
                     'label'     , label        , ...
                     'tag'       , V1(ii).Tag   , ...
                     'separator' , 'off'       , ...
                     'checked'   , 'off'       , ...
                     'enable'    , 'on'        , ...
                     'visible'   , 'on'        , ...
                     'interruptible' , 'off'   , ...
                     'busyaction'    , 'cancel'          );

  CB = [ 'browser(''Activate'','  epsstr(V1(ii).ID) ');' ];

  V1(ii).History = cat( 2 , hst , V1(ii).ID );

  set( V1(ii).ID , 'userdata' , V1(ii) , ...
                   'callback' , CB           );


  if ~isempty(V(ii).Children)

       rootc = root + V1(ii).ID * ( root == 0 );

      [ h , tag , V1(ii).Children ] = makemenu( V1(ii).ID , rootc , fig , ...
                                       V1(ii).Nr , V1(ii).History , V(ii).Children );

       Hdl = cat( 1 , Hdl , V1(ii).ID , h );
       Tag = cat( 1 , Tag , { V(ii).Tag } , tag );

  else

       Hdl = cat( 1 , Hdl , V1(ii).ID );
       Tag = cat( 1 , Tag , { V(ii).Tag } );

  end

end
% ii


% Set new Objects to Children of Parent !!!

pud = get( parent , 'userdata' );

pud.Children = cat( 1 , pud.Children , cat(1,V1.ID) );

% Check for new Children of a Chapter or SubItem 
if ( root ~= 0 )

  if isempty(pud.Double)
  % Create a new Double

      pa = get(parent,'parent');
      ch = get(parent,'children');
      
      set( ch , 'parent' , pa );

      ziel = parent;
      vsbl = 'on';

      if root == parent
      % Double of a Chapter, Check if Activate
      %  pa == ContentsMenu
        cud = get(pa,'userdata');
        if isequal(cud.ActiveRoot,parent)
          ziel = pa;
          vsbl = 'off';
        end
      end
  
      pud.Double = copyobj( parent , ziel );

      set( ch , 'parent' , parent );

      set( pud.Double , 'label'     , pud.Label                , ...
                        'separator' , 'off'                    , ...
                        'callback'  , get(parent,'callback')   , ...
                        'userdata'  , []                       , ... 
                        'tag'       , cat(2,pud.Tag,'_DOUBLE') , ...
                        'visible'   , vsbl                         );

      set( parent     , 'callback'  , ''         );

      set( pud.Children(1) , 'separator' , 'on' );
       
  end

   if isequal( get(pud.Double,'parent') , parent )
     ch = cat( 1 , pud.Double , pud.Children );
     set( parent , 'children' , ch(end:-1:1) );
   end

end


set( parent , 'userdata' , pud );



%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [txt,V] = gettext( HMsg , id , txt , mode , ...
                            ch , ct , cini , color , Matlab , nh );

% GETTEXT  Get Text from Text and File for ContentsMenu
%           works recursivly
 
  Nout = nargout;

    nl = { '' };

   mud = get(id,'userdata');

  if isempty(txt)
     txt = cell(0,1);
  end


  %------------------------------------------------
  % ChildrenNumber 

  nc = prod(size(mud.Children));


  %------------------------------------------------
  % Check for Recurse
  %
  %  No Children  or ...
  %     Chapter Selected  &  ( Mode == 1 ) 
  %   ==>  don't recurse

  is_recurse =  ~(   ( nc == 0 )  |  ...
                   ( ( mud.Root == 0 )  &  ( mode == 1 ) ) );


  %------------------------------------------------
  % Header

  titel = mud.Titel;
  titel = titel( 1+strcmp(titel(1),'_') : end );

  NrStart = mode + 1; % !!!
  NrEnde  = prod(size(mud.Nr));

  head = [ sprintf('%.0f.',mud.Nr(NrStart:NrEnde))  ' '  titel ]; 

  % Number of HeaderCharacters
  nh   = max( nh , size(head,2) );


 
  %------------------------------------------------
  % Text  &  File

  txt1 = cell(0,1);
  ref1 = cell(0,2);

  for ff = { 'Text'  'File' }

     val = getfield(mud,ff{1});

     if ~isempty( val )

       %------------------------------------------------------------
       % Check for MatlabHelp

       msg0 = [ 'GETTEXT( '  head  ' )' ];
       
       [ok,hlp] = helptext( HMsg , Matlab.Help , val , msg0 );
      
       if ok
       
         if ~isempty(hlp)       
               txt1 = cat( 1 , txt1 , nl , hlp , nl );
         end

       elseif strcmp(ff{1},'File')
       %------------------------------------------------------------
       % Read File

          msg = [ 'GETTEXT( '  head  ' ), ReadFile( '  val  ' )'  ];
 
          msg_list(HMsg,'Message',msg,'new');

          [Msg,bb,rr] = readfile( val , '' , cini.RefMark );

          if ~isempty(Msg)

             msg_list(HMsg,'Message',Msg,'append');

             bb = cat( 1 , nl , {[ 'Error reading File: ' val ]} );

             try
               [m,Msg] = char2cell( Msg , '' );
               if isempty(m)
                 bb = cat( 1 , bb , nl , Msg , nl );
               end
             end 

          else

            ref1 = cat(1,ref1,rr);
 
          end


             % Remove CommentLines
	     n  = size(cini.Comment,2);
	     ok = char(bb);
	     if size(ok,2) >= n
	        ok = ( double(ok(:,1:n)) == ...
		       cini.Comment(ones(1,size(ok,1)),:) );
                ok = find( ( sum(ok,2) < n ) );
                bb = bb(ok);
             end
	     		         
             txt1 = cat( 1 , txt1 , bb );


       else
       %------------------------------------------------------------
       %  Text

            try

               [Msg,bb,rr] = char2cell(val,cini.NewLine,cini.RefMark);

            catch

                Msg = 'Error';

            end

            if isempty(Msg);

               txt1 = cat( 1 , txt1 , bb );
               ref1 = cat( 1 , ref1 , rr );


            else

                msg = [ 'GETTEXT( '  head  ' ), CHAR2CELL( Text )'  Msg  ];
 
                msg_list(HMsg,'Message',msg,'new');

            end

       end

     end
     % ~isempty( val )

  end
  % ff


  %------------------------------------------------
  % Note

   note = cell(0,1);

   if ~isempty(mud.Note)

       try

          [Msg,note,rr] = char2cell(mud.Note,cini.NewLine,cini.RefMark);

       catch

           Msg = 'Error';

       end

       if isempty(Msg);

          ref1 = cat( 1 , ref1 , rr );

           sep = cini.NoteSep( ones(1,size(char(note),2)) );
 
          note = cat( 1 , {sep} , note , nl );

       else

           msg = [ 'GETTEXT( '  head  ' ), CHAR2CELL( Note )'  Msg  ];
 
           msg_list(HMsg,'Message',msg,'new');

       end

   end


  %------------------------------------------------
  % Replace Reference
 
  ref = cell(0,2);

  if ~isempty(txt1)

    for ii = 1 : size(ref1,1)

       jj = find( strcmp( ct , ref1{ii,2} ) );

       if ~isempty(jj)

         rud = get( ch(jj) , 'userdata' );

         if iscell(rud)
           rud = cat(1,rud{:});
         end

         str = '';

         for rr = 1 : prod(size(rud))
          
           is_root = ( rud(rr).Root == 0 );

           n0 = mode * (~is_root) + 1; % !!!
           n1 = prod(size(rud(rr).Nr));

           str1 = cat( 2 , cini.RefText{1+is_root} , ' ' , ...
                           sprintf('%.0f.',rud(rr).Nr(n0:n1))  ); 
             
           if ~( rud(rr).Nr(1) == mud.Nr(1) )  &  ~( n0 == 1 )

              str1 = cat( 2 , cini.RefText{2} , ' ' , ...
                           sprintf('%.0f, ',rud(rr).Nr(1)) , str1 ); 

           end

           if mode == 1
              str2 = get(rud(rr).History,'userdata');
              if iscell(str2)
                 str2 = cat(1,str2{:});
              end
              str2 = cellstr(str2mat(str2.Label));
              str2 = cat( 2 , '<' , strhcat(str2,'>-<') , '>' );
              str1 = cat( 2 , str1 , ' (' , str2 , ')' );
           end

           str = cat( 2 , str , cini.RefSep(1:(end*(~isempty(str)))) , str1 );

           ref = cat( 1 , ref , { str1  rud(rr).Nr } );

         end
         % rr

         txt1 = strrep( txt1 , ref1{ii,1} , str );
         
       end
       % ~isempty(jj)

    end
    % ii

  end
  % ~isempty(txt1)


  %------------------------------------------------
  % Matlab

  ImageFile = '';

  if ~isempty(mud.Matlab)
  
    %------------------------------------------------------------
    % Check for MatlabHelp

    msg0 = [ 'GETTEXT( '  head  ' )' ];
       
    [ok,txt2] = helptext( HMsg , Matlab.Help , mud.Matlab , msg0 );


    if ~ok  &  ( isempty(txt1)  |  ( Nout > 1 )  ) 
    %------------------------------------------------------------
    % Check for Event
       
      [MatlabOk,par] = getpar( Matlab.Event , mud.Matlab );

      if MatlabOk

        typ = upper( Matlab.Event{MatlabOk} );

        if strcmp( typ , 'IMAGE' )
           ImageFile = par;
        end

        if isempty(txt1)

          switch typ
 
           case 'HELPWIN'

             if isempty(par)
               p0 = '';
               p1 = '';
             else
               p0 = [ ' for '  par         ];
               p1 = [ '('''    par   ''')' ];
             end
 
             txt2 = cat( 1 , ...
                    { [ ' To get the MatlabHelp' p0 ] }       , ...
                    { '  type in the MatlabCommandWindow:' }  , ...
                    nl                                        , ...
                    { [ '>> help'     p1 ] }                  , ...
                    { '   to display the MatlabHelp,' }       , ...
                    nl                                        , ...
                    { [ '>> helpwin'  p1 ] }                  , ...
                    { '   to open the MatlabHelpWindow,' } );

             if isempty(par)
               txt2 = cat( 1 ,  txt2 ,  nl , ... 
                       { '>> helpdesk' }    , ...               
                      { '   to open the WebBrowser with the MatlabOnlineHelp.' } );
             end
          

          otherwise

             if isempty(par)
               txt2 = { typ };
             else
               [m,txt2] = char2cell(par);
             end

             for ii = 1 : size(txt2,1)
                 txt2{ii} = cat( 2 , '>> ' , txt2{ii} );
             end

          end
          % switch 

        end
        % isempty(txt1)

      end
      % MatlabOk

    end
    % Check for Event
    
    if ~isempty(txt2)
   
      txt1 = cat( 1 , txt1 , nl , txt2 , nl );

    end

  end
  % Matlab
  

  %------------------------------------------------
  % Append to Header and Text1 to Text

  txt = cat( 1 , txt , { head } );

  if ~isempty(txt1)
    
    if ~isempty(note)
 
       txt1 = cat( 1 , txt1 , nl , note );

    end

     sep = cini.TitelSep( ones(1,nh) );

     txt = cat( 1 , txt , { sep } , txt1 , nl );

     nh  = 0;   % !!!!

  end
             



  %------------------------------------------------
  % StructureOutput

  if Nout > 1


     V = struct( 'Nr'       , { mud.Nr    } , ...
                 'Titel'    , { titel     } , ...
                 'Text'     , { txt1      } , ...
                 'Note'     , { mud.Note  } , ...
                 'Reference', { ref       } , ...
                 'Image'    , { ImageFile } , ...
                 'Contents' , { cell(0,2) } , ...
                 'Children' , { [] }              );


             nc = prod(size(mud.Children));

     V.Children = V(ones( nc , 1 ));

  end

  %------------------------------------------------

  if ~isempty(color)

    set( id , 'foregroundcolor' , color );

    if strcmp(computer,'PCWIN');
       set(id,'checked','on');
    end

  end


  if ~is_recurse
    
      return
  
  end

  %------------------------------------------------
  % Recurse trough Children

  if Nout < 2

    for id = mud.Children(:)'

       txt = gettext( HMsg , id , txt , mode , ...
                      ch , ct , cini , color , Matlab , nh );
      
    end

    if isempty(txt1)  &  ~isempty(note)
 
       txt = cat( 1 , txt , note , nl );

    end

    return

  end
  

  %------------------------------------------------

  for ii = 1 : nc

    [txt,V.Children(ii)] = gettext( HMsg , mud.Children(ii) , txt , mode , ...
                                    ch , ct , cini , color , Matlab , nh );

    if isempty(txt1)  &  ~isempty(note)
 
       txt = cat( 1 , txt , note , nl );

    end

  end



%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [V,c,nr,str]  = getcont(V)

% GETCONT  returns Table of Contents from TextStructure


Nout = nargout;

c   = cell(0,2);

nr  = cell(0,1);
str = cell(0,1);

if isempty(V)
  return
end

bl = ' ';

n = prod(size(V));

for ii = 1 : n

   c   = cat( 1 , c   , { V(ii).Nr  V(ii).Titel } );

   nr  = cat( 1 , nr  , { sprintf('%.0f.',V(ii).Nr) } );
   str = cat( 1 , str , { ...
             cat( 2 , bl(ones(1,2*(prod(size(V(ii).Nr))-1))) , V(ii).Titel ) } );

   if ~isempty(V(ii).Children)

      [ V(ii).Children , c1 , nr1 , str1 ] = getcont(V(ii).Children);

      c   = cat( 1 , c   ,   c1  );
      nr  = cat( 1 , nr  ,  nr1  );
      str = cat( 1 , str , str1  );
 
      V(ii).Contents = c1;

   end  

   if prod(size(V(ii).Nr)) == 1
   % Chapter

     nr  = cat( 1 , nr  , { '' } );
     str = cat( 1 , str , { '' } );

   end

end
% ii

if Nout < 4

   n  = size(nr,1);

   nr = cat( 2 , bl(ones(n,1)) , char(nr) , bl(ones(n,4)) , char(str) );

   sep = char( double('=') * ones(1,size(nr,2)) );

   nr = cat( 1 , { '' } , cellstr(nr) , { '' } , { sep } , { '' } );

end
 
%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [Msg,bb,ref] = readfile(file,mnl,ref01);

% READFILE  Read Text from File and Converts to CellStringArray

 MaxSize = 500000;  % MaximumFileSize [kB]

 Msg = '';

 bb  = cell(0,1);

 ref = cell(0,2);

 nl = char(10);

 % Marker for NewLine, using in char2cell
 if nargin < 2
   mnl = '';
 end

 if nargin < 3
   ref01 = { '' '' };
 end

 %----------------------------------------------

  fid = fopen(file,'r');

  if fid == -1  

     Msg = [ 'Can''t open File.' ];
     return

  end

 %----------------------------------------------
 % Check Size of File (max 200000 kB)

  d = dir( file );

  if isempty(d)
   d = dir( which(file) );  % which(file) gives the full Name
                           %  [ PathName FileName ]
  end

  if d.bytes > MaxSize
    Msg = [' File too large, Limit = '  ...
            sprintf('%.0f kB',MaxSize) '.' ];
    fclose(fid);
    return
  end


 %----------------------------------------------

  bb = fread(fid,'char');

  fclose(fid);


 %----------------------------------------------
 % Char to CellString

  try

    [Msg,bb,ref] = char2cell(bb,mnl,ref01);

    if ~isempty(Msg)
      Msg = [ Msg ' in File.' ];
    end

  catch

      Msg = ['READFILE: Error using CHAR2CELL.' ...
              nl lasterr ];  

  end


%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,txt] = helptext(HMsg,ini,val,Msg)

% HELPTEXT  returns HelpText, using HELP

ok  = 0;
txt = cell(0,1);

if ~( ischar(val)  |  iscellstr(val)  )  |  isempty(val) 
  return
end

if iscellstr(val)
  val = val(:);
end
    
if size(val,1) ~= 1
  return
end

  
[ok,par] = getpar( ini , char(val) );

if ok

   msg = [ Msg  ', Help( '  par  ' )'  ];
 
   msg_list(HMsg,'Message',msg,'new');

   try 

      txt = help(par);

      if isempty(txt)
         p = fileparts(par);
         if ~isempty(p)
             fs = filesep;
             ps = pathsep;
             mp = matlabpath;
             cp = seperate(mp,ps);
             ok = ( ~any(strcmp(cp,p)) & ~any( p == '@' ) );
             ok = ( ok & isempty(findstr([fs 'private'],p)) );
             if ok 
                matlabpath(cat(2,mp,ps,p));
                txt = help(par);
                matlabpath(mp);
             end
         end
      end

      if isempty(txt)
         txt = [ 'No Help available for '  par '!' ];
      end 

      [msg,txt] = char2cell(txt,'');

      if ~isempty(msg)
         txt = cell(0,1);
      end

   catch

      msg_list(HMsg,'Message','Error','append');

   end

end
% ok

%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = seperate(name,sep);

% SEPERATE  separates Name
%
% String = SEPERATE( Name , Seperator )
%
%


str = cell(1,0);

if isempty(name)
   return
end

n = size(name,2);

%---------------------------------------------
% Find Seperator in Name

is = ( double(name) == double(sep) );

if all(is)
   str    = cell(1,n-1);
   str(:) = { '' };
   return
end

%---------------------------------------------

i0 = ~is(1);
i1 = ~is(n);

is = cat( 2 , ones(1,i0) , is , ones(1,i1) );

is = find( is );

is = is(:);

ni = size(is,1) - 1;

if ni == 0 
   return
end
     
%---------------------------------------------
% [ Start  End ]

ind = ( 1 : ni ) ;

is  = cat( 2 , is(ind)+1 , is(ind+1)-1 ) - i0;

%---------------------------------------------

ni = size(is,1);

if ni == 0
   return
end

%---------------------------------------------
   
is = is(1:ni,:);

str = cell(1,ni);
  
for ii = 1 : ni
    str{ii} = name( is(ii,1) : is(ii,2) );
end

%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [Msg,bb,ref] = char2cell(bb,mnl,ref01);

% CHAR2CELL Converts CharArray to CellStringArray

 Msg = '';

 ref = cell(0,2);

 Nin  = nargin;
 Nout = nargout;
 
  if isempty(bb)
     bb = cell(0,1);
     return
  end

  % Marker for NewLine
  if Nin < 2
    mnl = '';
  end

  if ischar(bb)
     bb = double(bb);
  end

  if ( size(bb,1) > 1 )  &  ( size(bb,2) > 1 )
     bb = cat( 2 , bb , 10*ones(size(bb,1),1) );
     bb = permute( bb , [ 2 1 ] );
     bb = bb(:);
  end

  if ( size(bb,1) > 1 ) 
     bb = permute( bb , [ 2 1 ] );
  end

  %---------------------------------------------------
  % Check Characters

  ok = all( ( bb ==  9 ) |  ...
            ( bb == 10 ) |  ...
            ( bb == 13 ) |  ...
            (  28 <= bb  &   bb <= 126 ) | ...
            ( 160 <= bb  &   bb <= 255 )        );

  if ~ok
    Msg = ['Invalid Characters' ];
    return
  end
 

  %---------------------------------------------------
  % Remove CR
  bb( find( bb == 13 ) ) = [];


  bb = char(bb);


  %---------------------------------------------------
  % TAB  --> 8 Blanks 
  bb = strrep( bb , char(9) , char(32*ones(1,8)) ); 

  %---------------------------------------------------
  % mnl --> NewLine   % !!!!!
  if ~isempty(mnl)
    bb = strrep( bb , mnl , char(10) ); 
  end

  %---------------------------------------------------
  % Reference
  if ( Nin == 3 )  &  ( Nout == 3 )

     [MsgR,ref] = get_ref(bb,ref01{:});

  end

  %---------------------------------------------------
  % Form CellString

  % 1. "'"     --> "''"
  bb = strrep( bb , char(39) , char([39  39]) ); 

  % 2. NL --> "';'"
  bb = strrep( bb , char(10) , char([39  59  39]) );
  

  bb = eval([  '{'''  bb  '''}' ]);


%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,out] = get_ref(str,s1,s2,c1,c2)

% GET_REF  Returns a Reference from a String
%
% [Msg,Reference] = GET_REF( String , S1 , S2 , C1 , C2 );
%
%  S1:  StartMarker  Characters which marks the Begin of Reference
%  S2:    EndMarker  Characters which marks the End   of Reference
%
%  C1:  Characters in Reference to remove,
%        default: [ CR NL TAB ]
%  C2:  Characters closed to C1 to remove,
%        default: [ Space '-' ]
%
%  Reference = { FullString  ReferenceString }  2-Column Cell-Array
%
%

Msg = '';

out = cell(0,2);


nl = char(10);

Msg0 = 'GET_REF: ';

nm0 = size(Msg0,2);

nl0 = char([ 10 32*ones(1,nm0+0) ]);

%---------------------------------------------------------

if ~( ischar(str) &  ( prod(size(str)) == size(str,2) ) )
  Msg = '1. Input must be a String.';
end
 
if ~( ischar(s1) &  ( prod(size(s1)) == size(s1,2) ) )
  Msg = [ Msg nl0(1:(end*(~isempty(Msg)))) ...
           'S2 must be a String.' ];
end
 
if ~( ischar(s2) &  ( prod(size(s2)) == size(s2,2) ) )
  Msg = [ Msg nl0(1:(end*(~isempty(Msg)))) ...
           'S2 must be a String.' ];
end

%---------------------------------------------------------

if nargin < 4

  c1 = [ 13  10  9 ];

else

  if ischar(c1)
    c1 = double(c1);
  end
  ok = isnumeric(c1);
  if ok & ~isempty(c1)
    c1 = c1(:)';
    ok = all( ( mod(c1,1) == 0 )  & ...
              ( c1 >= 0 ) & isfinite(c1)  );
  end
  if ~ok
      Msg = [ Msg nl0(1:(end*(~isempty(Msg)))) ...
              'Input C1 must be a String or ASCII-Codes.'];
  end

end   

%---------------------------------------------------------

if nargin < 5

  c2 = [ double('-')  32  ];

else

  if ischar(c2)
    c2 = double(c2);
  end
  ok = isnumeric(c2);
  if ok & ~isempty(c2)
    c2 = c2(:)';
    ok = all( ( mod(c2,1) == 0 )  & ...
              ( c2 >= 0 ) & isfinite(c2)  );
  end
  if ~ok
      Msg = [ Msg nl0(1:(end*(~isempty(Msg)))) ...
              'Input C1 must be a String or ASCII-Codes.'];
  end

end   

%---------------------------------------------------------

if ~isempty(Msg)
  Msg = [ Msg0  Msg ];
end

if isempty(str)  |  isempty(s1)  |  isempty(s2)  |  ...
   ~isempty(Msg)
 
  return

end


%**********************************************************

n = size(str,2);


%--------------------------------
% Start of Reference at End of s1

i1 = findstr(str,s1);

if isempty(i1)
  return
end

i1 = i1 + size(s1,2) - 1;
i1( find( i1 > n ) ) = [];

if isempty(i1)
  return
end
 
%--------------------------------
% End of Reference at Begin of s2

i2 = findstr(str,s2);

if isempty(i2)
  return
end

i2( find( i2 < i1(1) ) ) = [];

if isempty(i2)
  return
end

%----------------------------

ok = zeros(1,n);  % Start/End
ii = zeros(1,n);  % Index

ok(i1) = 1;  % Start of Reference
ok(i2) = 2;  % End   of Reference

ii(i1) = i1;
ii(i2) = i2;

jj = find(ok);

ok = ok(jj);
ii = ii(jj);


%----------------------------
% Following Start-End

ok = find( abs(diff(ok)) == 1 );

if isempty(ok)
   return
end

%----------------------------
% Start's

ok = ok(1:2:end);

n = prod(size(ok));

out      = cell(n,2);
out(:,1) = { [ s1  s2 ] };
out(:,2) = { '' };


for jj = 1 : n

  if ( ii(ok(jj)+1) - ii(ok(jj)) ) > 1
  % Reference not empty

    % String in Reference
 
    sr = str( ii(ok(jj))+1 : ii(ok(jj)+1)-1 );


    % Full String

    out{jj,1} = [ s1  sr  s2 ];


    % Remove bad Characters from Reference

    sr = rmblank( sr , 2 );

    b1 = zeros( 1 , size(sr,2) );
    for r1 = c1
      b1 = ( b1  |  ( double(sr) == r1 ) );
    end

    b2 = zeros( 1 , size(sr,2) );
    for r2 = c2
      b2 = ( b2  |  ( double(sr) == r2 ) );
    end
  

    i0 = find( diff(cat(2,0,( b1 | b2 )  )) ==  1 );  % Start of Group
    i1 = find( diff(cat(2,  ( b1 | b2 ),0)) == -1 );  % End   of Group
 
     b = 2 * b1 + 1 * b2;
    cb = cumsum(b,2);

    lg = i1 - i0 + 1;              % Length of Group
    sg = cb(i1) - cb(i0) + b(i0);  % Sum    of Group         

    ib = find( sg > lg );          % Groups incl. b1

    lg = lg(ib);
    i0 = i0(ib);

    % Build IndexVector to remove

     b = grp2ind(i0,lg);

     sr(b) = [];
 
     out{jj,2} = sr;

  end
  % Reference not empty

end
% jj



%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ii = grp2ind(i0,l);

% GRP2IND  Built IndexVector from StartIndex and Length
%
% Index = GRP2IND( StartIndex , GroupLength )
%

if isempty(i0);
   ii = [];
   return
end

si = size(i0);

if ( sum( si > 1 ) > 1 )
   error('StartIndex must be a Vector.');
end

i0 = i0(:);
l  =  l(:);

if ~isequal(size(i0,1),size(l,1))
   error('Size of StartIndex and GroupLenght must be the same.');
end

n = size(l,1);

ii = ones(sum(l),1);
jj = cumsum( cat(1,1,l) , 1 );

ii(jj(1:n)) = i0;

if n > 1
   ii(jj(2:n)) = ii(jj(2:n))-(i0(1:n-1,1)+l(1:n-1)-1);
end

ii = cumsum(ii,1);

if n == 1
   return
end

jj = find( si > 1 );

if jj == 1
   return
end
  
perm = cat( 2 , (1:jj-1)+1 , 1 , (jj+1:size(si,2)) );

ii = permute(ii,perm);


%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,par] = getpar(ini,str);

% GETPAR(INI,String)  Get Parameter from String, starts with INI
%
%  INI can be a String, CharrArray or CellStringArray
%

Nout = nargout;

ok  = 0;
par = '';

if iscellstr(ini);
   ini = char(ini(:));
end

ini = rmblank(ini,2);
str = rmblank(str,2);

if isempty(str)
  return
end

n1 = size(ini,1);
n2 = size(ini,2);

ns = size(str,2);

ok1 = ( ischar(ini)  &  ( prod(size(ini)) == n1*n2 ) );
ok2 = ( ischar(str)  &  ( prod(size(str)) == ns )  );

if  ~ok1  |  ~ok2  
  
  error('Invalid Inputs.');

end

%----------------------------------------------

if isempty(ini);

   ok = 0;

   return

end

%----------------------------------------------

n = min(n2,ns);

sn = sum( ( ini(:,1:n) == str(ones(1,n1),1:n) ) , 2 );

ok = ( sn == n );

if any(ok)

   ok = find(ok);
   
else
 
   ok = zeros(n1,1);
   for ii = 1 : n1
     in = rmblank( ini(ii,:) , 2 );
     if ~isempty(in)
         n2     = size(in,2);
         n      = min( n2 , ns );
         sn     = sum( ( ini(ii,1:n) == str(1:n) ) , 2 );
         ok(ii) = n2 * ( sn == n );
     end
   end

   if ~any(ok)
     ok     = 0;
   else
    [n2,ok] = max(ok);
   end

end
  
if ~ok  |  ( Nout < 2 )
   return
end

%----------------------------------------------

if ns > n2
    
    par = rmblank( str(n2+1:ns) , 2 );
  
end



%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = event(HMsg,typ,par,titel,tag)

% EVENT  executes an MatlabEvent
%

nl = char(10);

ok  = 0;

Msg = '';

fig = get(HMsg,'parent');


switch upper(typ)

  case 'HELPWIN'

    if isempty(par)

         msg_list(HMsg,'Message','HELPWIN','new');

         try
           helpwin;
           ok = 1;
         end

    else

         msg_list(HMsg,'Message',[ 'HELPWIN( '  par  ' )'  ],'new');

         try
           helpwin(par);
           ok = 1;
         end

    end

    if ~ok

         msg_list(HMsg,'Message','Error','append');

         Msg = [ Msg0 'Error call HELPWIN('''  par  ''')' nl nl lasterr ];

    end

  %------------------------------------------------------------------
  case 'EVAL'

    if ~isempty(par)

         msg_list(HMsg,'Message',[ 'EVAL( '  par  ' )'  ],'new');

         try

           eval(par);
           ok = 1;

         catch

           msg_list(HMsg,'Message','Error','append');

           Msg = [ Msg0 'Error call EVAL('''  par  ''')' nl nl lasterr ];

        end

    end
  
  %------------------------------------------------------------------
  case 'IMAGE'

    if ~isempty(par)

         msg_list(HMsg,'Message',[ 'IMAGE( '  par  ' )'  ],'new');

         try

           Msg = showimage(fig,par,titel,tag);
           ok = 1;

         catch

           msg_list(HMsg,'Message','Error','append');
           Msg = [ Msg0 'Error call IMAGE('''  par  ''')' nl nl lasterr ];

        end

    end

end
% switch  


    if ~isempty(Msg)

        browser('MsgBox',fig,Msg,'Warning','warn');          

    end

 
%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function Msg = showimage(parent,file,name,tag)

Msg = '';

 nl = char(10);

%-----------------------------------------------------------

try

  [ c , cmap ] = imread(file);

catch

  Msg = [ 'Can''t read Image: '  file   nl nl lasterr ];
 
  return

end

%-----------------------------------------------------------

if isempty(cmap);
  cmap = ones(1,3);
end

%-----------------------------------------------------------

si = size(c);

si(2) = si(2) + 2; % !!!

xlim = [ 1  si(2) ] + 0.5 * [-1  1 ];
ylim = [ 1  si(1) ] + 0.5 * [-1  1 ];


puni = get(parent,'units');
       set(parent,'units','pixels');
ppos = get(parent,'position');
       set(parent,'units',puni);

scr_uni = get(0,'units');      set(0,'units','pixels')
scr_si  = get(0,'ScreenSize'); set(0,'units',scr_uni);


figpos        = zeros(1,4);
figpos([3 4]) = si([2 1]);

figpos(1) = scr_si(3)-40-figpos(3);
figpos(2) = ppos(2)+ppos(4)-60-figpos(4);
 
axepos        = ones(1,4);
axepos([3 4]) = si([2 1]);


%-----------------------------------------------------------

 fig = figure( 'paperorientation'   , 'portrait' , ...
               'position'           , figpos     , ...
               'color'              , [1 1 1]    , ...
               'numbertitle'        , 'off'      , ...
               'menubar'            , 'none'     , ...
               'toolbar'            , 'none'     , ...
               'name'               , name       , ...
               'colormap'           , cmap       , ...
               'visible'            , 'on'       , ...
               'handlevisibility'   , 'callback' , ...
               'createfcn'          , ''         , ...
               'resize'             , 'off'             );

%-----------------------------------------------------------

pap_si = get(fig,'papersize');
ppi    = get(0,'screenpixelsperinch');

pappos([3 4]) = figpos([3 4]) / ppi;
pappos([1 2]) = (pap_si-pappos([3 4])) / 2;

set(fig,'paperunits'    , 'inches'   , ...
        'paperposition' , pappos            );

%-----------------------------------------------------------


axe = axes( 'parent'   , fig            , ...
            'units'    , 'pixels'       , ...
            'position' , axepos         , ...
            'xlim'     , xlim           , ...
            'ylim'     , ylim           , ...
            'xtick'    , []             , ...
            'ytick'    , []             , ...
            'color'    , 'none'         , ...
            'xgrid'    , 'off'          , ...
            'ygrid'    , 'off'          , ...
            'xdir'     , 'normal'       , ...
            'ydir'     , 'reverse'      , ...
            'dataaspectratio' , [ 1  1  1 ]   , ...
            'visible'            , 'off'      , ...
            'handlevisibility'   , 'callback'         );


image( 'parent' , axe , ...
       'xdata'  , [ 1  si(2) ] , ...
       'ydata'  , [ 1  si(1) ] , ...
       'cdata'  , c            , ...
       'cdatamapping' , 'direct'  , ...
       'visible'      , 'on'               );

               
%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  delmsgbox( tag );

% DELMSGBOX  Deletes all MessageBoxes, specified by Tag


  shh = get(0,'showhiddenhandles');
        set(0,'showhiddenhandles','on');

  delete( findobj( 0 , 'tag' , tag ) );

        set(0,'showhiddenhandles',shh);


%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
 

%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [Msg,V] = structchk(V,txt);

% STRUCTCHK( V )   Checks if V is an StructArray,
%
%  or CellArray with Field-Value-Pairs and build it.
%
%  [ Msg , Structure ] = STRUCTCHK( V )
%
%  returns the Error's in Msg
%    or the Stucture build from V
%

Msg = '';

 nl = char(10);

if nargin == 2
   txt = [ ' after '  txt  '-Action' ];
else
   txt = '';
end


V = V(:);

if isstruct(V)
  return
end

NV = size(V,1);

if     NV == 0

   Msg = ['Empty Input'  txt  '.'];

elseif NV == 1
%---------------------------------------------
% Single Input must be StructArray

  if iscell(V)
   V = V{1};
  end

  if ~isstruct(V)
    Msg = ['Single Input'  txt  ' must be a StructArray.' ];
    return
  end

  V = V(:);

else
%---------------------------------------------
% FieldName-Value-Pairs ==>  build StructArray

  if mod(NV,2) ~= 0
    Msg = ['Number of Multiple Input'   txt  ' must be even' ...
           ' (Field-Value-Pairs).' ];
    return
  end

  ok = iscellstr(V(1:2:end-1));  % FieldNames
  msg = '';

  if ok

    try 
      V = struct(V{:});
    catch
      msg = [ nl lasterr];
    end

    ok = isempty(msg);

  end

  if ~ok
    Msg = ['Multiple Input'  txt  ' must be Field-Value-Pairs.'  ...
            nl(1:(end*(~isempty(msg))))  msg ];
    return
  end

  V = V(:);

end

%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
 n = 10;
end
if Nin < 2
 del = char((32*ones(1,3)));
end


if isempty(str)
 str = '';
 return
end


if ~( ischar(str)  |  iscellstr(str) )
   error('StringArray must be a CharArray or CellStringArray.');
end

if iscellstr(str)
  str = char(str);
end

str = double(str);

    jj    = find( sum( ( str == 32 ) , 2 ) == size(str,2) );
str(jj,:) = [];
 
str = cellstr(char(str));


str = str(:);

str( : , 2 ) = { del };

str(n:n:size(str,1),2) = {nl};

str(    size(str,1),2) = {''};


str = str';

str = cat(2,str{:});


%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = rmblank(str,dim,cc)

% RMBLANK  Remove Blanks, NewLines at Begin and End of CharacterArrays
%
% String = RMBLANK( CharArray )
%
% CharArray  2-dimensional CharacterArray
%
% further Options:
%
% String = RMBLANK( CharArray , DIM , CHAR )
%
%  
%  DIM  specifies Dimension to work, 
%       default: 2
%
%  CHAR specifies BlankCharacters to remove
%       default:  [ 32  13  10  9 ];  % [ Space CR LF TAB ]
%

  
msg = '';
 nl = char(10);

Nin = nargin;

if Nin < 1
  error('Not enough Input Arguments.')
else
  str0 = str;
  if ischar(str)
    str = double(str);
  end
  ok = isnumeric(str);
  if ok
    ok = all( ( mod(str(:),1) == 0 )  & ...
              ( str(:) >= 0 ) & isfinite(str(:))  );
  end
  if ~ok
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CharArray must be a String or ASCII-Codes.'];
  end
  if size(str,1)*size(str,2) ~= prod(size(str))
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CharArray must be 2-dimensional.'];
  end     
end

if Nin < 2
  dim = 2;
else
  if ~isnumeric(dim)
    msg = [ msg nl(1:(end*(~isempty(msg)))) ...
            'Input DIM must be numeric.' ];
  elseif ~isempty(dim)
    dim = dim(:)';
    if ~all( ( dim == 1 ) |  ( dim == 2 ) )
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
             'Values for Input DIM must Integers larger ZERO.' ];
    end
  end 
end

if Nin < 3
  cc = [ 32  13  10  9 ];  % [ Space CR LF TAB ]
else
  if ischar(cc)
    cc = double(cc);
  end
  ok = isnumeric(cc);
  if ok & ~isempty(cc)
    cc = cc(:)';
    ok = all( ( mod(cc,1) == 0 )  & ...
              ( cc >= 0 ) & isfinite(cc)  );
  end
  if ~ok
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CHAR must be a String or ASCII-Codes.'];
  end
end

if ~isempty(msg)
  error(msg)
end


if isempty(str)
 str = '';
 return
end

if isempty(dim) | isempty(cc)
  str = str0;
  return
end


     jj  = find(str == 0 );
 str(jj) = cc(1);

  blank  = 0*str;

  for ii = cc
    blank = ( blank | ( str == ii ) );
  end

  si = size(str);

  for d = dim

    bad = ( sum(blank,3-d) == si(3-d) );
    jj  = find( bad );
    if ~isempty(jj) 

         p  = [ 3-d  d ];
        str = permute(str,p);

         jj = jj(:)';
         nb = size(jj,2);
        jj1 = find( jj ==   ( 1 : nb ) );       % Blank at Begin
        jj2 = find( jj == ( ( 1 : nb ) + ...    % Blank at End
                            ( si(d) - nb ) ) );
        str(:,jj([jj1 jj2])) = [];

        str = permute(str,p);

    end
    
  end

  str = char(str);



%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function txt = intro

%  INTRO   returns DefaultHelpString


txt = { ...
''
'           __-----_.                        ______'
'          /  \      \           o  O  O   _(      )__'
'         /    |  |   \_---_   o._.      _(           )_'
'        |     |            \   | |""""(_   Let''s see... )'
'        |     |             |@ | |    (_               _)'
'         \___/   ___       /   | |      (__          _)'
'           \____(____\___/     | |         (________)'
'           |__|                | |          |'
'           /   \-_             | |         |'''
'         /      \_ "__ _       !_!--v---v--"'
'        /         "|  |>)      |""""""""|'
'       |          _|  | ._--""||        |'
'       _\_____________|_|_____||________|_'
'      /                                   \'
'     /_____________________________________\'
'     /                                     \'
'    /_______________________________________\'
'    /                                       \'
'   /_________________________________________\'
'        {                               }'
'        <_______________________________|'
'        |                               >'
'        {_______________________________|               ________'
'        <                               }              / SNOOPY \'
'        |_______________________________|             /__________\'
'\|/       \\/             \||//           |//                       \|/    |/'
'kharas@expert.cc.purdue.edu'
};


%*******************************************************************
% The End
%*******************************************************************

