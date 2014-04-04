function  [Msg,fig] = mkgui(Config,name,bg)

% MKGUI  Creates a GUI-Figure from Structure
%
% [ Msg , fig ] = MKGUI( Config , FigureName , BackGroundColor )
%
% For an example to build the ConfigureStructure and
%  call MKGUI look into the Directory:  "example/"
%
%  type:   >> [Msg,fig] = EXAMPLE_gui
%
%
% requires: MSG_LOGO, TAG_FRAME
%
% see also: STDGUI, TAG_FRAME, TAB_LIST, SEL_LIST, MSG_LIST, MSG_LOGO, MKMENU
%
% 

Msg = '';
fig = '';

nl = char(10);

Msg0 = 'MKGUI: ';


if nargin < 1
   Msg = [ Msg0 'Input Config is missing.' ];
   return
end

if nargin < 2
   name = '';
end

if nargin < 3
   bg = [ 1  1  1 ];
end

%*******************************************

 form = 'GUI %2.2d/%2.2d/%4.0f %2.0f:%2.2d:%2.2d';
        
   cl = clock;
  ind = [3 2 1 4 5 6];

  tag = sprintf(form,round(cl(ind))) ;

%*******************************************
% PrePositioning of Figure
 
[fs,scr_si] = bestfont;

figpos = [ NaN  NaN  ceil( 2/3 * scr_si([3 4])) ];
figpos(1) = 50;
figpos(2) = scr_si(4)-60-figpos(4);


fig = figure( 'paperunits'   , 'inches'     , ...
              'paperorientation','portrait' , ...
              'units'        , 'pixels'     , ...
              'position'     , figpos       , ...
              'color'        , [1 1 1]      , ...
              'menubar'      , 'none'       , ...
              'toolbar'      , 'none'       , ...
              'numbertitle'  , 'off'        , ...
              'name'         ,  name        , ...
              'colormap'     , [ 1  1  1 ]  , ...
              'createfcn'    , ''           , ...
              'tag'          , tag          , ...
              'resize'       , 'off'        , ...
              'visible'      , 'off'        , ...
              'integerhandle'    , 'on'     , ...
              'handlevisibility' , 'callback'           );

%******************************************************************

 fud = struct( 'Children' , { ...
               struct( 'Menu'     , { [] } , ...
                       'Frame'    , { [] }        ) } );


%******************************************************************

field = fieldnames(Config);

is_menu    = any(strcmp(field,'Menu'));
is_logo    = any(strcmp(field,'Logo'));
is_control = any(strcmp(field,'Control'));


if is_menu
   is_menu = ~isempty(Config.Menu);
end

if is_logo
   is_logo = ~isempty(Config.Logo);
end

if is_control
   is_control = ~isempty(Config.Control);
end


%******************************************************************
% UIMenu's

if is_menu

  try
    [Msg,HM] = mkmenu(fig,'Menu',Config.Menu,fig);
  catch
    Msg = lasterr;
  end

  if ~isempty(Msg)
    Msg = [ Msg0  'Error call MKMENU.' nl Msg ];
    return
  end

  fud.Children.Menu = HM;

end

%******************************************************************

if ~( is_logo | is_control )

  set( fig , 'userdata' , fud );
 
  return

end



%******************************************************************
% Frame arround

try
  [Msg,Frame,hb,pos,bh] = tag_frame(fig,'new','TagNumber' , 0 , ...
                                 'Position'        , [ 0  0 -0 -0 ] , ...
                                 'BackGroundColor' ,   bg           , ...
                                 'ForeGroundColor' , [ 0  0  0 ]          );
catch
  Msg = lasterr;
end


if ~isempty(Msg)
  Msg = [ Msg0  'Error call TAG_FRAME( New ).' nl Msg ];
  return
end

fud.Children.Frame = Frame;

set( fig , 'userdata' , fud );


%******************************************************************
% FrameChildren

 fch = struct( 'Logo'    , { [] } , ...
               'Control' , { [] }      );

%******************************************************************
% Logo

pos0 = [ 0  -1 ];

if is_logo

  try
    [Msg,hl,hf,ht] = tag_frame(Frame,'add',0,'msg_logo',Config.Logo{:});
  catch
    Msg = lasterr;
  end

  if ~isempty(Msg)
    Msg = [ Msg0  'Error call TAG_FRAME( Add  "msg_logo" ).' nl Msg ];
    return
  end

  pos = get(hl,'position');

  pos0(2) = -1 * ( figpos(4) - pos(2) + 1 );

  fch.Logo = struct( 'Frame'   , { hf } , ...
                     'Text'    , { ht } , ...
                     'Message' , { hl }       );

end

%******************************************************************
% UIControls

if is_control

  cnf = make_cnf(fig);

  try
    [Msg,fch.Control,pos] = mkuic(Frame,0,pos0,Config.Control,cnf);
  catch
    Msg = lasterr;
  end

  if ~isempty(Msg)
    Msg = [ Msg0  'Error call MKUIC.' nl Msg ];
    return
  end

end


%******************************************************************
% Fit FigurePosition

 ud = get(Frame,'userdata');

 pos = pos + 2 * ud.TAG_FRAME.BorderWidth;

 figpos(2) = figpos(2) + (figpos(4)-pos(2)) + 30*(~is_menu);

 figpos([3 4]) = pos;

 orient = { 'portrait'  'landscape' };
 orient = orient{ 1 + ( figpos(3) >= figpos(4) ) };

 set( fig , 'paperorientation' , orient , ...
            'position'         , figpos       );

 [Msg,pos] = tag_frame(Frame,'resize','absolut');

 if ~isempty(Msg) , return, end

 wygiwys(fig);

%******************************************************************


 ud.Children = fch;


 set( Frame , 'userdata' , ud );



%******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,HM] = mkmenu(par,grp,cfg,fig);

%  FieldValues of Config for UIMenu-Definition:
%
% UIMenu:
%    { Label  SeparatorOn   CallbackFcn   UserData }
%    { Label  SeparatorOn   SubStructure  UserData }
%
%   SeperatorOn = [ 0  |   1 ];
%
%   UserData are optional
%
% UIContextMenu:
%    {  NaN   SubStructure  CallbackFcn   UserData }
%


%-----------------------------------------------------------

Msg = '';
HM  = []; 

nl = char(10);

if isempty(cfg)
  return
end


%************************************************************


field = fieldnames(cfg);
field = field(:)';

    n = size(field,2);

%------------------------------------------------------------
% HandleStructure

HM      = cell( 2 , n+1 );
HM(1,:) = cat( 2 , field , {'Children'} );
HM(2,:) = { {[]} };

HM      = struct( HM{:} );

HC      = zeros( n , 1 );

%------------------------------------------------------------

HP = epsstr(par);

sets = { 'off'  'on' };

for ii = 1 : n

  cc = getfield( cfg , field{ii} );

  nc = prod(size(cc));

  CB    = '';
  usd   = [];

  if nc < 3
     cc = cat( 2 , cc , { '' } );  % CallBackFcn
  end
  if nc >= 4
     usd = cc{4};
  end


  if ischar( cc{3} )  &  ~isempty(cc{3})
       CB = [ cc{3} '('  HP  ','''  grp  ''',''' field{ii} ''',1,'  ...
              sprintf('%.0f',ii)  ');'  ];
  end

  tag = cat( 2 , grp , '.' , field{ii} );

  if ischar(cc{1})

    %----------------------------------
    % UIMenu

    enb = sets{ 1 + (~isempty(rmblank(cc{1},2))) };
    sep = sets{ 1 + isequal(cc{2},1) };
 
    HC(ii) = uimenu( ...
       'parent'          , par      , ...
       'callback'        , CB       , ...
       'userdata'        , usd      , ...
       'label'           , cc{1}    , ...
       'tag'             , tag      , ...
       'enable'          , enb      , ...
       'separator'       , sep      , ...
       'checked'         , 'off'    , ...
       'Interruptible'   , 'off'    , ...
       'BusyAction'      , 'cancel'  );

     ch = cc{3};

   else

     %----------------------------------
     % UIContextMenu

     HC(ii) = uicontextmenu( ...
       'parent'          , fig      , ...
       'callback'        , CB       , ...
       'userdata'        , usd      , ...
       'tag'             , tag      , ...
       'Interruptible'   , 'off'    , ...
       'BusyAction'      , 'cancel'  );
 
     ch = cc{2};

   end 

  if isstruct(ch)
  % Children

    try
       [Msg,ch] = mkmenu(HC(ii),field{ii},ch,fig);
    catch
        Msg = lasterr;
    end

    if ~isempty(Msg) 
       Msg = [ tag  ': Error call MKMENU( '  field{ii}  ' ).' ...
               nl Msg ];
       return
    end

  end  

  ud = struct( 'Root'     , { get(par,'parent') } , ...
               'Parent'   , { par } , ...
               'Children' , { ch  } , ...
               'UserData' , { usd }       );

  set( HC(ii) , 'userdata' , ud );

  HM = setfield( HM , field{ii} , HC(ii) );

end

HM.Children = HC;


%******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


function [Msg,HC,pos_out] = mkuic(par,nr,pos0,cfg,cnf);


Msg = [];
HC  = [];

pos_out = [0  0];

if isempty(cfg)
  return
end

%************************************************************
% ParentHandle for CallBack

ud = get( par , 'userdata' );

if nr == 0

   HB = par;

else

   ud = get( par , 'userdata' );

   HB = ud.TAG_FRAME.ButtonHandle(nr);

end


HP = epsstr(HB);



%************************************************************

% Fit for single Characters

is_win = strcmp( upper(computer) , 'PCWIN' );


hd  = cnf.Control.HSpace;   % Horizontal Distance between
vd  = cnf.Control.VSpace;   % Vertical Distance between


fcp = get(par,'ForeGroundColor');  % ForeGroundColor

bcp = get(par,'BackGroundColor');  % TextLabelBackGroundColor

% bcp = 'r';     % Use this to check Extension of TextBox !!!

ct = cnf.Control;  % Config for TextLabel above

ct.SingleHeight = 1 * ct.PixelFit(1,2) + ct.PixelFit(2,2);
ct.FontWeight   = 'bold';

ListBoxExt = '_ListBoxFrame';  % Extension for Frame arround ListBox

%************************************************************

field0 = fieldnames(cfg);
field0 = field0(:)';

    ng = size(field0,2);

%------------------------------------------------------------
% HandleStructure

HC      = cell( 2 , ( ng ) );
HC(1,:) = field0;
HC(2,:) = { {[]} };

HC      = struct( HC{:} );

%************************************************************

pos0   = pos0 + [ hd  -vd ];

wmax_g = 0;   % Max. OffsetWidth


%************************************************************
%------------------------------------------------------------
for iig = 1 : ng
% Group

  cfg1 = getfield( cfg , field0{iig} );

%**********************************************************
if isstruct( cfg1)


  field01 = fieldnames( cfg1 );  % Original

  field1 = field01(:)';          % to expand with fields of TAG_FRAME's

      nc = size(field01,1);

  cc = struct2cell( cfg1 );
  cc = cat(1,cc{:});

  %------------------------------------------------------------
  % Expand HandleStructure
  % 
  % Check for Text only

  text_only = 1;

  for iic = nc : -1 : 1

     field2 = field1(iic);

     typ = lower(cc{iic,2}.Type);

     text_only = ( text_only & strcmp(typ,'text') );

     switch typ 

       case 'tag_frame' 
       % Add TagButton-Fields

         field2 = fieldnames(cc{iic,2}.Option);

         field2 = field2(:)';

       case 'listbox'
       % Add Frame

         field2  = cat( 2 , field1(iic) , ...
                        { cat( 2 , field01{iic} , ListBoxExt ) } ); 

     end


     field1 = cat( 2 , field1(1:iic-1) , field2 , ...
                       field1(iic+1:end) );


  end

  %------------------------------------------------------------
  % Structure of Group

  H1      = cell( 2 , size(field1,2) );
  H1(1,:) = field1;
  H1(2,:) = { {[]} };

  H1      = struct( H1{:} );


  pp = cat(1,cc{:,1});    % Positions
  hp = real(pp);          % HorizNr.

  pos_h = pos0;   % Start for Horizontal

  vmin_h = 0;     % Max. OffsetHight

  nh = max(hp);

  %----------------------------------------------
  % Go Horizontal 
  for iih = 1 : nh

    ih = find( hp == iih );  % Same horizontal Position

    vp = imag( pp(ih) );    % VertNr.

    nv = max(vp);
       
    pos_v = pos_h;  % Start for Vertical

    wmax_v = 0;  % Max. OffsetWidth
 
    %--------------------------------------------
    % Go Vertical
    for iiv = 0 : nv

      iv = find( vp == iiv );

      nn = size(iv,1);

      vmin_n = 0;  % Max. OffsetHight
      
      
      %--------------------------------------------
      % Controls at same Position
  
      for iin = 1 : nn

          iic = ih(iv(iin));  % Index in cc

          ccn = cc{iic,2};

          typ = lower( ccn.Type );

        %--------------------------------------------------------
        % Check Width / Position
 
          if isnan(real(ccn.Width))

             ccn.Width = size( char(ccn.String) , 2 ) + i * imag(ccn.Width);

          end

        %--------------------------------------------------------
        % Check Color

          % BackGroundColor

          if isnan(ccn.Color)
             ccn.Color = bcp;
          end

          % ForeGroundColor

          fc0 = fcp;  

          if isequal( ccn.Color , fc0 )

             fc0 = 1 - fc0;
 
          end

        %--------------------------------------------------------
        % Font's  &  Positioning

          switch typ
           case 'tag_frame'
              c = cnf.Tag;
           case { 'tab_list'  'sel_list'  'listbox' }
              c = cnf.List;
           case 'edit'
              c = cnf.Edit;
           case 'text'
              c = cnf.Text;
           case 'button'
              c = cnf.Button;
              typ = 'pushbutton';
           case 'popupmenu'
              c = cnf.Control;
              c.SliderWidth = cnf.List.SliderWidth;
           otherwise
              c = cnf.Control;
          end

        %--------------------------------------------------------
        % Check for Different Width
     
        c = setstyle(c,ccn.Style,cnf);


        %--------------------------------------------------------

         pos_n = pos_v;

        
        %--------------------------------------------------------
        switch typ

          %*******************************************************************
          case 'tag_frame'

          % Note: TagFrames only 1 in Horizontal !!!!


            % CallBack
            if ~isempty(ccn.CBFcn)

              CB = { ccn.CBFcn HB  field0{iig} };

            else
 
              CB = cell(0,1);

            end


             fieldt = fieldnames(ccn.Option);
                 nt = size(fieldt,1)-1;  
                 nb = prod(size(ccn.String));

               ntag = max(nt,nb);
 
%  'Position'    , [  Left  Bottom  -Right -Top  ] , default: [ 5  5 -5 -5]
%                  [  Left  Bottom   Width  High ] , ... and all Combinations ...
%                  [ -Right -Top     Width  High ]



             % First in Group  or Single in Group ==> LeftOffset == TagFrameOffset
             pos_n(1) = pos_n(1) + ( c.HSpace - hd ) * ( ( iig == 1 ) | ( nc == 1 ) );

             % First Group & Single in Group  ==> TopOffset == TagFrameOffset

             pos_n(2) = pos_n(2) + ( -c.VSpace + vd ) * ...
                                   ( ( iig == 1 )  &  ( nc == 1 ) );

               pos = [ pos_n  -c.HSpace  2*ceil(3*cnf.Text.SingleHeight+2*c.BorderWidth) ];


             VarIn = { 'TagNumber'       , max(nt,nb)    , ...
                       'CBFcn'           , CB            , ...
                       'Position'        , pos           , ...
                       'BorderWidth'     , c.BorderWidth , ...
                       'ForeGroundColor' , fcp           , ...
                       'BackGroundColor' , ccn.Color     , ...
                       'FontName'        , c.FontName    , ...
                       'FontUnits'       , c.FontUnits   , ...
                       'FontSize'        , c.FontSize    , ...
                       'FontWeight'      , c.FontWeight  , ...
                       'Interruptible'   , 'off'         , ...
                       'BusyAction'      , 'cancel'             };
                                            
             [Msg,h,hb,pos,bh] = tag_frame(par,'add',nr,'tag_frame',VarIn{:});

             if ~isempty(Msg) , return, end


             % Set ButtonStrings
             for ib = 1 : nb 
                 set(hb(ib),'string',ccn.String{ib});
             end

             posb = zeros( nt+1 , 2 );  %  [ Width  Hight ]

             hh = cat(1,h,hb);

             % Set ButtonChildren
             for it = 0 : nt 
         
                [Msg,H2,posb(it+1,:)] = mkuic( h , it , [0 0] , ...
                                      getfield(ccn.Option,fieldt{it+1}), cnf );

                 if ~isempty(Msg) , return, end


                 H1 = setfield( H1 , fieldt{it+1} , hh(it+1) );

                 ud = get( hh(it+1) , 'userdata' );

                 ud.Children = H2;


                 CB = '';

                 if ~isempty(ccn.CBFcn)  &  ( it > 0 )
      
                   CB = [ ccn.CBFcn '('  HP ','''  field0{iig}   ''','''  ...
                                               fieldt{it+1} ''',1);' ];

                 end

                 set( hh(it+1) , 'userdata' , ud , ...
                                 'callback' , [ get(hh(it+1),'callback')  CB ] );  

                 
             end  

             %------------------------------------------------------------ 

               ud = get( h , 'userdata' );

             %------------------------------------------------------------ 
             % Set Button to Children

               for it = 1 : nt

                  ud.Children = setfield( ud.Children , fieldt{it+1} , ...
                                          ud.TAG_FRAME.ButtonHandle(it)    );

               end           

             %------------------------------------------------------------ 
             % Adjust Frame in Position in Hight !!!

               pos    = max( posb , [] , 1 ) + 2*ud.TAG_FRAME.BorderWidth;
               pos(2) = pos(2) + bh;

               ud.TAG_FRAME.Position(4) = pos(2);

             %------------------------------------------------------------ 

               set( h , 'userdata' , ud );

             %------------------------------------------------------------ 
             % Adjust OriginalPosition in UserData of Parent-Frame !!!

               pud = get(par,'userdata');
               
               if nr == 0

                  jj = find( pud.TAG_FRAME.FrameChildren.Handle == h );

                  pud.TAG_FRAME.FrameChildren.Position{jj}(4) = pos(2);
 
                  set(par,'userdata',pud);

               else

                   hud = get( pud.TAG_FRAME.HideHandle(nr) , 'userdata' );

                    jj = find( hud.Handle == h );

                   hud.Position{jj}(4) = pos(2);

                   set( pud.TAG_FRAME.HideHandle(nr) , 'userdata' , hud )

               end
 

             %------------------------------------------------------------ 

               % Last Group & Single in Group  ==> Offset == TagFrameOffset
               
               last = ( ( iig == ng )  &  ( nc == 1 ) );

               % wmax_v will added to pos_h, move offs same like pos_n above !!!
               % Last Group & Single in Group  ==> LeftOffset == TagFrameOffset

               c.HSpace = (~last) * hd + ...
                            last  * ( 2*c.HSpace - hd );


               c.VSpace = (~last) * vd + ...
                            last  * c.VSpace ;

     
               pos_n(2) = pos_n(2) - pos(2);


               pos = pos([1 2 1 2]);  % Need pos(3) below !!!


          %*******************************************************************
          case 'sel_list'


            % CallBack
            if ~isempty(ccn.CBFcn)

              CB = { ccn.CBFcn HB  field0{iig} };

            else
 
              CB = cell(0,1);

            end

            hht  = ceil( ct.SingleHeight  * (~isempty(ccn.Text)) );
            pos = pos_n + [ 0  -hht ];

            [ Msg , h ] = tag_frame(par,'add',nr,'sel_list' , ...
                      'fontname'    , c.FontName   , ...
                      'fontunits'   , c.FontUnits  , ...
                      'fontsize'    , c.FontSize   , ...
                      'fontweight'  , c.FontWeight , ...
                       ccn.Option{:}               , ...
                      'ListFont'    , c.FontName   , ...
                      'position'    , pos          , ...
                      'CBFcn'       , CB           , ...
                      'foregroundcolor' , [0 0 0]  , ...
                      'backgroundcolor' , [1 1 1]  , ...
                      'Interruptible'   , 'off'    , ...
                      'BusyAction'      , 'cancel'        );

            if ~isempty(Msg) , return, end


            H1 = setfield( H1 , field01{iic} , h );


            %------------------------------------------------------
            % TextLabel above
            if ~isempty(ccn.Text) 

                pos        = get(h,'position');
                pos([1 2]) = pos_n;
                pos(3)     = floor(pos(3));
                pos(4)     = hht;
 
               [Msg,htext] = tag_frame(par,'add',nr,'uicontrol' , ...
                           'position'    ,  pos          , ...
                           'style'       , 'text'        , ...
                           'string'      , ccn.Text      , ...
                           'fontname'    , ct.FontName   , ...
                           'fontunits'   , ct.FontUnits  , ...
                           'fontsize'    , ct.FontSize   , ...
                           'fontweight'  , ct.FontWeight , ... 
                           'cdata'       , []            , ...
                           'foregroundcolor'     , fcp   , ...
                           'backgroundcolor'     , bcp   , ...
                           'horizontalalignment' , 'center'         , ...
                           'callback'  , ''     );

              if ~isempty(Msg) , return, end

            end

            pos = get(h,'position');

            pos_n(2) = pos_n(2) - ceil( pos(4) + hht );

               
          %*******************************************************************
          case 'tab_list'


            hht  = ceil( ct.SingleHeight  * (~isempty(ccn.Text)) );

            pos = pos_n + [ 0  -hht ];

            [ Msg , h ] = tag_frame(par,'add',nr,'tab_list' , ...
                      'fontname'    , c.FontName   , ...
                      'fontunits'   , c.FontUnits  , ...
                      'fontsize'    , c.FontSize   , ...
                      'fontweight'  , c.FontWeight , ...
                       ccn.Option{:}               , ...
                      'position'    , pos          , ...
                      'foregroundcolor' , [0 0 0]  , ...
                      'backgroundcolor' , [1 1 1]  , ...
                      'Interruptible'   , 'off'    , ...
                      'BusyAction'      , 'cancel'               );

            if ~isempty(Msg) , return, end


            H1 = setfield( H1 , field01{iic} , h );


            ud = get( h , 'userdata' );
            set( ud.cntr{1} , 'backgroundcolor' , ccn.Color );


            % CallBack
            if ~isempty(ccn.CBFcn)

              CB = { ccn.CBFcn HB  field0{iig} };

              Msg = tab_list( h , 'callback'    , CB , ...
                                  'setquestfcn' , CB , ...
                                  'delquestfcn' , CB      );
             
              if ~isempty(Msg) , return, end

            end
 
            %------------------------------------------------------
            % TextLabel above
            if ~isempty(ccn.Text) 

 
                pos        = get(h,'position');
                pos([1 2]) = pos_n;
                pos(3)     = ceil(pos(3));
                pos(4)     = hht;

               [Msg,htext] = tag_frame(par,'add',nr,'uicontrol' , ...
                           'position'    ,  pos          , ...
                           'style'       , 'text'        , ...
                           'string'      , ccn.Text      , ...
                           'fontname'    , ct.FontName   , ...
                           'fontunits'   , ct.FontUnits  , ...
                           'fontsize'    , ct.FontSize   , ...
                           'fontweight'  , ct.FontWeight , ... 
                           'cdata'       , []            , ...
                           'foregroundcolor'     , fcp   , ...
                           'backgroundcolor'     , bcp   , ...
                           'horizontalalignment' , 'center'         , ...
                           'callback'  , ''     );

              if ~isempty(Msg) , return, end

            end

            %------------------------------------------------------
            % Search for following ControlButtons, cc(:,1) == NaN

            ind = ( iic+1 : nc );
 
            jj = find( cumprod( double( isnan( hp(ind) ) ) ) );

            if ~isempty(jj)
      
               ind = ind(jj);

               for jj = ind

                 ccb = cc{jj,2};

                 if isnan(ccb.Color)
                    ccb.Color = bcp;
                 end

                 fc1 = get( h , 'foregroundcolor' );
                 if isequal( fc1 , ccb.Color )
                    fc1 = 1 - fc1;
                 end

                 CB = '';

                 if ~isempty(ccb.CBFcn)
       
                   CB = [ ccb.CBFcn '('  HP ','''  field0{iig} ''',''' ...
                                               field01{jj} ''',1);' ];

                 end

                 if ~ischar(ccb.Option)
                    ccb.Option = '';
                 end

                 [Msg,hb] = tab_list(h,'button', ...
                            'style'       , ccb.Type      , ...
                            'string'      , ccb.String    , ...
                            'cdata'       , []            , ...
                            'value'       , ccb.Value(1)  , ...
                            'min'         , ccb.Value(2)  , ...
                            'max'         , ccb.Value(3)  , ...
                            'sliderstep'  , ccb.Value([4 5])/diff(ccb.Value([2 3])) , ...
                            'userdata'    , ccb.UserData      , ...
                            'tooltipstring'       , ccb.Option , ...
                            'foregroundcolor'     , fc1       , ...
                            'backgroundcolor'     , ccb.Color , ...
                            'horizontalalignment' , 'center'  , ...
                            'callback'            , CB        , ...
                            'Interruptible'       , 'off'     , ...
                            'BusyAction'          , 'cancel'      );

                  if ~isempty(Msg) , return, end

                  H1 = setfield( H1 , field01{jj} , hb );

               end
               % jj

            end
            %  ~isempty(jj)

            pos = get(h,'position');

            pos_n(2) = pos_n(2) - ceil( pos(4) + hht );


          %*******************************************************************
          otherwise
          % uicontrol


            CB = '';

            if ~isempty(ccn.CBFcn)
      
               CB = [ ccn.CBFcn '('  HP ','''  field0{iig}  ''',''' ...
                                           field01{iic} ''',1);' ];

            end


            %------------------------------------------------------
            % Position

            cw = real(ccn.Width);  % CharacterWidth
            ch = imag(ccn.Width);  % CharacterHight

 
            ww0 = ceil( cw * c.PixelFit(1,1) + c.PixelFit(2,1) + ...
                        c.LineOffset(1) + c.SliderWidth  );
 
            hh0 =  ch * ( c.PixelFit(1,2) + c.LineOffset(2) ) + ...
                   c.PixelFit(2,2);

            hh0 = ceil( ( ch >  0 ) * hh0 + ...
                        ( ch == 0 ) * c.SingleHeight  );

            %------------------------------------------------------
            % PixelPosition if NEGATIVE !!!

            ww0 = ww0 - ( ww0 + cw ) * ( cw < 0 );
            hh0 = hh0 - ( hh0 + ch ) * ( ch < 0 );

            ww0 = ww0 + 2 * c.BorderWidth;
            hh0 = hh0 + 2 * c.BorderWidth;

            %------------------------------------------------------
            % SeparatorPosition
            if strcmp(typ,'separator');

               typ = 'frame';

               ww0 = ~( cw == 0 ) * ww0  + ...
                      ( cw == 0 ) * cnf.Separator.SingleHeight;

               hh0 = ~( ch == 0 ) * hh0  + ...
                      ( ch == 0 ) * cnf.Separator.SingleHeight;

               if ~isempty(ccn.Style)

                  if strcmp( lower(ccn.Style) , 'list' )

                    ww0 = ww0 + ~( cw == 0 ) * 2 * cnf.Control.BorderWidth;
                    hh0 = hh0 + ~( ch == 0 ) * 2 * cnf.Control.BorderWidth;

                  end

               end


               fc0       = fcp;   % ForeGroundColor
               ccn.Color = bcp;   % BackGroundColor


            %------------------------------------------------------
            elseif ~strcmp(typ,'listbox')

               if ~isempty(ccn.Style)

                  if strcmp( lower(ccn.Style) , 'list' )

                    ww0 = ww0 + 2 * cnf.Control.BorderWidth;

                  end

               end

            end
            %------------------------------------------------------

            % HorizontalAlignment
            ht = 'left';
            if   strcmp(typ,'pushbutton')  | ...
               ( strcmp(typ,'edit') & ( ch == 0 ) )
              ht = 'center';
            end


            %------------------------------------------------------
            % Frame around if ListBox !!!

            bw = cnf.Control.BorderWidth * strcmp( typ , 'listbox' );


            hht  = ceil( ct.SingleHeight  * (~isempty(ccn.Text)) );

            %------------------------------------------------------
            % TextLabel above
            if ~isempty(ccn.Text)

                pos = [ pos_n+[bw 0]  ww0 hht ];

               [Msg,htext] = tag_frame(par,'add',nr,'uicontrol' , ...
                           'position'    ,  pos          , ...
                           'style'       , 'text'        , ...
                           'string'      , ccn.Text      , ...
                           'fontname'    , ct.FontName   , ...
                           'fontunits'   , ct.FontUnits  , ...
                           'fontsize'    , ct.FontSize   , ...
                           'fontweight'  , ct.FontWeight , ... 
                           'cdata'       , []            , ...
                           'foregroundcolor'     , fcp   , ...
                           'backgroundcolor'     , bcp   , ...
                           'horizontalalignment' , 'center'         , ...
                           'callback'            , ''   , ...
                           'tag' , [ field0{iig} '.' field01{iic}  '_Text' ]   );

               if ~isempty(Msg) , return, end

               pos_n(2) = pos_n(2) - pos(4);

            end
 

            %------------------------------------------------------

            if ~ischar(ccn.Option)
               ccn.Option = '';
            end


            %------------------------------------------------------

            pos = [ pos_n  ww0  hh0 ];

            %------------------------------------------------------
            % Move vertical if Text single in Vertical

            if strcmp( typ , 'text' )

              dh = ( 2*cnf.Control.BorderWidth + cnf.Control.SingleHeight - ...
                     c.SingleHeight ) / 2 * strcmp( typ , 'text' ) ;
 
              dh = ( 1 - 1/2*text_only ) * dh;

              pos(2) = pos(2) - dh * ( 1 + ( nv == 0 ) * ( nn == 1 ) );

              pos_n(2) = pos_n(2) + dh;

            end

            %------------------------------------------------------
            % Frame around if ListBox !!!

            if strcmp( typ , 'listbox' )

               %------------------------------------------------------
               % Search for following Buttons, cc(:,1) == NaN

               dh = 0;

               ind = ( iic+1 : nc );
 
                jj = find( cumprod( double( isnan( hp(ind) ) ) ) );
 
               if ~isempty(jj)

                 ind = ind(jj);

                 for jj = ind

                     ccb = cc{jj,2};

                    typ1 = lower(ccb.Type);

                   %--------------------------------------------
                   switch typ1
                      case 'listbox'
                        c1 = cnf.List;
                      case 'edit'
                        c1 = cnf.Edit;
                      case 'text'
                        c1 = cnf.Text;
                      otherwise
                        c1 = cnf.Control;
                   end

                   %--------------------------------------------
                   % Check for use Height
                   if ~isempty(ccb.Style)

                      n = size(ccb.Style,2);

                      if strcmp( ccb.Style(n) , upper(ccb.Style(n)) );

                          ccb.Style = cat( 2 , lower(ccb.Style(1:n-1)) , ccb.Style(n) );

                          c1 = setstyle( c1 , ccb.Style , cnf );

                      end

                   end

                    ch1 = imag(ccb.Width);  % CharacterHight

                     hb = ch1 * ( c1.PixelFit(1,2) + c1.LineOffset(2) ) + c1.PixelFit(2,2);

                     hb = ceil( ( ch1 >  0 ) * hb + ...
                                ( ch1 == 0 ) * c1.SingleHeight );

                     hb = hb + ( c1.SliderHeight - hb ) * strcmp(typ1,'slider');

                     hb = hb - ( hb + ch1 ) * ( ch1 < 0 );  % PixelPosition if Negative

                     hb = hb + 2 * c1.BorderWidth;

                     dh = dh + hb;

                  end
                  % jj

                end
                % ~isempty(jj)

                posf        = pos;
                posf([3 4]) = posf([3 4]) + [ 0  dh ] + 2*bw;
 
                pos_n       = pos_n + [ 1  -1 ] * 2*bw;

                pos([1 2])  = pos([1 2]) + [ 1  -1 ] * bw;

                [Msg,hf] = tag_frame(par,'add',nr,'uicontrol' , ...
                            'position'    , posf         , ...
                            'style'       , 'frame'      , ...
                            'string'      ,  ''          , ...
                            'cdata'       , []           , ...
                            'userdata'    , []           , ...
                            'tooltipstring'       , ''   , ...
                            'foregroundcolor'     , fc0       , ...
                            'backgroundcolor'     , ccn.Color , ...
                            'callback'            , ''        , ...
                            'Interruptible'       , 'off'     , ...
                            'BusyAction'          , 'cancel'  , ...
                            'tag'  , [ field0{iig} '.' field01{iic}  ListBoxExt ]  );

               if ~isempty(Msg) , return, end

               H1 = setfield( H1 , [ field01{iic}  ListBoxExt ] , hf );
                
            end
            % listbox-frame

            %------------------------------------------------------

            [Msg,h] = tag_frame(par,'add',nr,'uicontrol' , ...
                            'position'    , pos          , ...
                            'style'       , typ          , ...
                            'string'      , ccn.String   , ...
                            'fontname'    , c.FontName   , ...
                            'fontunits'   , c.FontUnits  , ...
                            'fontsize'    , c.FontSize   , ...
                            'fontweight'  , c.FontWeight , ...
                            'cdata'       , []           , ...
                            'value'       , ccn.Value(1) , ...
                            'min'         , ccn.Value(2) , ...
                            'max'         , ccn.Value(3) , ...
                            'sliderstep'  , ccn.Value([4 5])/diff(ccn.Value([2 3])) , ...
                            'userdata'    , ccn.UserData      , ...
                            'tooltipstring'       , ccn.Option , ...
                            'foregroundcolor'     , fc0       , ...
                            'backgroundcolor'     , ccn.Color , ...
                            'horizontalalignment' , ht        , ...
                            'callback'            , CB        , ...
                            'Interruptible'       , 'off'     , ...
                            'BusyAction'          , 'cancel'  , ...
                            'tag'                 , [ field0{iig} '.' field01{iic} ]  );


             if ~isempty(Msg) , return, end


             H1 = setfield( H1 , field01{iic} , h );


             pos_n(2) = pos_n(2) - pos(4);


             %------------------------------------------------------
             % Search for following Buttons, cc(:,1) == NaN

             ind = ( iic+1 : nc );
 
             jj = find( cumprod( double( isnan( hp(ind) ) ) ) );
 
             if ~isempty(jj)
      
               ind = ind(jj);

               for jj = ind

                   ccb = cc{jj,2};

                 typ = lower(ccb.Type);

                 if isnan(ccb.Color)
                    ccb.Color = bcp;
                 end

                 fc1 = fcp;
                 if isequal( fc1 , ccb.Color )
                    fc1 = 1 - fc1;
                 end

                 %-------------------------------------
                 switch typ
                    case 'listbox'
                      c1 = cnf.List;
                    case 'edit'
                      c1 = cnf.Edit;
                    case 'text'
                      c1 = cnf.Text;
                    otherwise
                      c1 = cnf.Control;
                 end

                 %-------------------------------------
                 % Check for use Height
                 if ~isempty(ccb.Style)

                    n = size(ccb.Style,2);

                    if strcmp( ccb.Style(n) , upper(ccb.Style(n)) );

                        ccb.Style = cat( 2 , lower(ccb.Style(1:n-1)) , ccb.Style(n) );

                        c1 = setstyle( c1 , ccb.Style , cnf );

                    end

                 end

                  pos(2) = pos_n(2);


                  cw = real(ccb.Width);  % CharacterWidth
                  ch = imag(ccb.Width);  % CharacterHight

                  pos(4) = ch * ( c1.PixelFit(1,2) + c1.LineOffset(2) ) + c1.PixelFit(2,2);

                  pos(4) = ceil( ( ch >  0 ) * pos(4) + ...
                                 ( ch == 0 ) * c1.SingleHeight  );

                  pos(4) = pos(4) + ( c1.SliderHeight - pos(4) ) * ...
                                       strcmp(typ,'slider');

                  pos(4) = pos(4) - ( pos(4) + ch ) * ( ch < 0 );  % PixelPosition

                  pos(4) = pos(4)  + 2 * c1.BorderWidth;

                  pos_n(2) = pos_n(2) - pos(4);
  

                 CB = '';

                 if ~isempty(ccb.CBFcn)

                   CB = [ ccb.CBFcn '('  HP ','''  field0{iig} ''',''' ...
                                               field01{jj} ''',1);' ];

                 end


                 if ~ischar(ccb.Option)
                    ccb.Option = '';
                 end


                 [Msg,hb] = tag_frame(par,'add',nr,'uicontrol' , ...
                            'position'    , pos           , ...
                            'style'       , typ           , ...
                            'string'      , ccb.String    , ...
                            'fontname'    , c1.FontName   , ...
                            'fontunits'   , c1.FontUnits  , ...
                            'fontsize'    , c1.FontSize   , ...
                            'fontweight'  , c1.FontWeight , ...
                            'cdata'       , []            , ...
                            'value'       , ccb.Value(1)  , ...
                            'min'         , ccb.Value(2)  , ...
                            'max'         , ccb.Value(3)  , ...
                            'sliderstep'  , ccb.Value([4 5])/diff(ccb.Value([2 3])) , ...
                            'userdata'    , ccb.UserData      , ...
                            'tooltipstring'       , ccb.Option , ...
                            'foregroundcolor'     , fc1       , ...
                            'backgroundcolor'     , ccb.Color , ...
                            'horizontalalignment' , ht        , ...
                            'callback'            , CB        , ...
                            'Interruptible'       , 'off'     , ...
                            'BusyAction'          , 'cancel'  , ...
                            'tag'                 , [ field0{iig} '.' field01{jj} ] );

                 if ~isempty(Msg) , return, end

                 H1 = setfield( H1 , field01{jj} , hb );

               end                 
               % jj

             end
             % ~isempty(jj)

             pos(3) = pos(3) + 2*bw;  % Used for wmax           
      
        end
        % typ


        drawnow

        %*******************************************************************

        wm =  pos(3) + c.HSpace + ( hd - c.HSpace ) * ...
                                  ( hd < c.HSpace ) * ...
                         ~strcmp(typ,'tag_frame') * ( iih == nh );

        vm = pos_n(2) - c.VSpace;
  
        wmax_v = max( wmax_v , wm );
        vmin_n = min( vmin_n , vm );

        %*******************************************************************

      end
      % iin

      pos_v(2) = vmin_n;            % Move vertical before next V-Element
    
    end    
    % iiv, Go Vertical

    pos_h(1) = pos_h(1) + wmax_v;   % Move horizontal before next H-Element

    vmin_h = min( vmin_h , pos_v(2) );
      
  end
  % iih, Go Horizontal 

  wmax_g = max( wmax_g , pos_h(1) );

  pos0(2) = vmin_h;    % Move Vertical before next Group


%**********************************************************
elseif isnumeric(cfg1)

  %---------------------------------------------------------
  % Separator

  ok = ( prod(size(cfg1)) == 1 );
  if ok
     ok = isnan(cfg1);
  end

  if ok

     c = cnf.Separator;

     pos  = [ c.HSpace  pos0(2)  -c.HSpace  c.SingleHeight ];

     [Msg,H1] = tag_frame(par,'add',nr,'uicontrol', ...
               'position'        ,  pos          , ...
               'style'           , 'frame'       , ...
               'backgroundcolor' ,  bcp          , ...
               'foregroundcolor' ,  fcp          , ...
               'tag'             ,  'Separator'        );

      if ~isempty(Msg) , return, end


      pos0(2) = pos0(2) - pos(4) - vd;

  end
        
end
%**********************************************************

  HC = setfield( HC , field0{iig} , H1 );

end
% iig, Group

%------------------------------------------------------------ 
% Return Position

pos_out = [ wmax_g  -pos0(2) ];

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function c = setstyle(c,style,cnf);

% Style  use Width
% stylE  use Hight
% StylE  use Width and Hight
% STYLE  use all  
% STYLe  use all, except SingleHeight  

if isempty(style)
   return
end

n = size(style,2);

UseWidth = strcmp( style(1) , upper(style(1)) );

UseHight = strcmp( style(n) , upper(style(n)) );

UseAll   = strcmp( style(2:n-1) , upper(style(2:n-1)) );
   
try

  cf = fieldnames(cnf);
  jj = find( strcmp( lower(style) , lower(cf) ) );

  if ~isempty(jj)

      c1 = getfield( cnf , cf{jj(1)} );

      if UseAll

         h = c.SingleHeight;

         c = c1;

         if ~UseHight

             c.SingleHeight = h;

         end

      else

         if UseWidth

             c.SliderWidth   = c1.SliderWidth;
             c.PixelFit(:,1)  = c1.PixelFit(:,1);
             c.LineOffset(1)  = c1.LineOffset(1);
             c.HSpace         = c1.HSpace + 2*(c1.BorderWidth-c.BorderWidth);

         end

         if UseHight

             c.SliderHeight   = c1.SliderHeight;
             c.PixelFit(:,2)  = c1.PixelFit(:,2);
             c.LineOffset(2)  = c1.LineOffset(2);
             c.SingleHeight   = c1.SingleHeight;
             c.VSpace         = c1.VSpace + 2*(c1.BorderWidth-c.BorderWidth);

         end

      end

  end

end 

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function cnf = make_cnf(fig)

% MAKE_CNF   Font-, PixelConfiguration for UIControls


%************************************************

[fs,scr_si,ppi,is_win,fn] = bestfont;


FontName   = 'helvetica';
FontUnits  = 'points';
FontSize   =  fs;
FontWeight = 'normal';

        fs = fs - is_win;

SliderWidth  = 0;
SliderHeight = 10 + is_win;
BorderWidth  = 3;

 TagFontSize    = fs + 2;
 TagFontWeight  = 'bold';
 TagFrameOffset = 3;

ListFontSize    = FontSize;
ListFontName    = fn;

ListSliderWidth = 15 + ( 3 + 3*(ppi>96) ) * is_win;
ListLineOffset  = 2 - [ 0  2 ] * is_win;

TextFontSize    = fs + 2;
TextFontWeight  = 'bold';
TextBorderWidth = 0;

SeparatorHeight  = 3;
SeparatorOffset  = 5;


z2  = zeros(1,2);
z22 = zeros(2,2);

%************************************************

%-----------------------------------------------
Tag = struct(  ...
               ...
'FontName'     , {    FontName     } , ...
'FontUnits'    , {    FontUnits    } , ...
'FontSize'     , { TagFontSize     } , ...
'FontWeight'   , { TagFontWeight   } , ...
'BorderWidth'  , {  BorderWidth    } , ...
'LineOffset'   , { z2              } , ...
'SliderWidth'  , {  SliderWidth    } , ...
'SliderHeight' , {  SliderHeight   } , ...
'PixelFit'     , { z22             } , ...
'SingleHeight' , { 0               } , ...
'HSpace'       , { TagFrameOffset  } , ...
'VSpace'       , { TagFrameOffset  }        );

%-----------------------------------------------
List = struct( ...
               ...
'FontName'     , { ListFontName    } , ...
'FontUnits'    , {     FontUnits   } , ...
'FontSize'     , { ListFontSize    } , ...
'FontWeight'   , {     FontWeight  } , ...
'BorderWidth'  , {     BorderWidth } , ...
'LineOffset'   , { ListLineOffset  } , ...
'SliderWidth'  , { ListSliderWidth } , ...
'SliderHeight' , { ListSliderWidth } , ...
'PixelFit'     , { z22             } , ...
'SingleHeight' , { 0               } , ...
'HSpace'       , { 0               } , ...
'VSpace'       , { 0               }        );

%-----------------------------------------------
Text = struct( ...
               ...
'FontName'     , {     FontName    } , ...
'FontUnits'    , {     FontUnits   } , ...
'FontSize'     , { TextFontSize    } , ...
'FontWeight'   , { TextFontWeight  } , ...
'BorderWidth'  , { TextBorderWidth } , ...
'LineOffset'   , { z2              } , ...
'SliderWidth'  , {     SliderWidth } , ...
'SliderHeight' , { 0               } , ...
'PixelFit'     , { z22             } , ...
'SingleHeight' , { 0               } , ...
'HSpace'       , { 0               } , ...
'VSpace'       , { 0               }        );

%-----------------------------------------------
Edit = struct( ...
               ...
'FontName'     , { FontName    } , ...
'FontUnits'    , { FontUnits   } , ...
'FontSize'     , { FontSize    } , ...
'FontWeight'   , { FontWeight  } , ...
'BorderWidth'  , { BorderWidth } , ...
'LineOffset'   , { z2          } , ...
'SliderWidth'  , { SliderWidth } , ...
'SliderHeight' , { 0           } , ...
'PixelFit'     , { z22         } , ...
'SingleHeight' , { 0           } , ...
'HSpace'       , { 0           } , ...
'VSpace'       , { 0           }        );

%-----------------------------------------------
Control = struct( ...
                  ...
'FontName'     , { FontName     } , ...
'FontUnits'    , { FontUnits    } , ...
'FontSize'     , { FontSize     } , ...
'FontWeight'   , { FontWeight   } , ...
'BorderWidth'  , { BorderWidth  } , ...
'LineOffset'   , { z2           } , ...
'SliderWidth'  , { SliderWidth  } , ...
'SliderHeight' , { SliderHeight } , ...
'PixelFit'     , { z22          } , ...
'SingleHeight' , { 0            } , ...
'HSpace'       , { 0            } , ...
'VSpace'       , { 0            }        );

%-----------------------------------------------
Button = struct( ...
                 ...
'FontName'     , {     FontName    } , ...
'FontUnits'    , {     FontUnits   } , ...
'FontSize'     , { TextFontSize    } , ...
'FontWeight'   , { TextFontWeight  } , ...
'BorderWidth'  , { BorderWidth  } , ...
'LineOffset'   , { z2           } , ...
'SliderWidth'  , { SliderWidth  } , ...
'SliderHeight' , { 0            } , ...
'PixelFit'     , { z22          } , ...
'SingleHeight' , { 0            } , ...
'HSpace'       , { 0            } , ...
'VSpace'       , { 0            }        );

%-----------------------------------------------
Separator = struct( ...
                    ...
'FontName'     , { FontName     } , ...
'FontUnits'    , { FontUnits    } , ...
'FontSize'     , { FontSize     } , ...
'FontWeight'   , { FontWeight   } , ...
'BorderWidth'  , { 0            } , ...
'LineOffset'   , { z2           } , ...
'SliderWidth'  , { SliderWidth  } , ...
'SliderHeight' , { 0            } , ...
'PixelFit'     , { z22          } , ...
'SingleHeight' , { SeparatorHeight } , ...
'HSpace'       , { SeparatorOffset } , ...
'VSpace'       , { 0               }        );


%************************************************

cnf = struct( 'Tag'       , { Tag       } , ...
              'List'      , { List      } , ...
              'Text'      , { Text      } , ...
              'Edit'      , { Edit      } , ...
              'Control'   , { Control   } , ...
              'Button'    , { Button    } , ...
              'Separator' , { Separator }       );
              

%***************************************************
% Right Positioning of UIControls

  % DummyAxes

   axe = axes( 'parent'   , fig          , ...
               'units'    , 'normalized' , ...
               'position' , [ 0 0 1 1 ]  , ...
               'visible'  , 'off'              );
  % DummyText

  ht = text('parent'     , axe        , ...
            'units'      ,'pixels'    , ...
            'position'   , [ 1 1 0 ]  , ...
            'string'     , ''         , ...
            'interpreter', 'none'     , ...
            'visible'    , 'off'               );


%------------------------------------------------

for ff = { 'Tag'  'List'  'Text'  'Edit'  'Control'  'Button' }

    c = getfield( cnf , ff{1} );

    set( ht , 'fontname'   , c.FontName   , ...
              'fontunits'  , c.FontUnits  , ...
              'fontsize'   , c.FontSize   , ...
              'fontweight' , c.FontWeight        );

    ext = zeros(2,4);

    for ii = [ 1  2 ]

      set( ht , 'string' , char( double('H') * ones(ii,ii) ) );

      ext(ii,:) = get( ht , 'extent' );

    end

    m =   ext(2,[3 4]) - ext(1,[3 4]);
    n = ( ext(1,[3 4]) - m );

    c.PixelFit = cat( 1 , m , n );

    cnf = setfield( cnf , ff{1} , c );

end

delete(ht);
delete(axe);

%***************************************************

c = cnf.Control.PixelFit;

hd   = ceil( 1.2 * ( 1 * c(1,1) + c(2,1) ) );  % Horizontal Distance between
vd   = ceil( 2/3 * ( 1 * c(1,2) + c(2,2) ) );  % Vertical Distance between



cnf.List.PixelFit(2,2) = 0;  % List: no vertical Offset


for ff = { 'Separator'  'List'  'Text'  'Edit'  'Control'   'Button' }

    c = getfield( cnf , ff{1} );

    switch ff{1}

     %------------------------------------------
     case 'List'

       c.HSpace = 2*hd;
       c.VSpace =   vd;

       c.SingleHeight = c.PixelFit(1,2) + c.LineOffset(2);

     %------------------------------------------
     case { 'Control'  'Text'  'Edit'  'Button' }

       c.HSpace = hd + hd * strcmp( ff{1} , 'Button' );
       c.VSpace = vd;

       f = 1 + 0.2 * ( strcmp( ff{1} , 'Control' ) * (~is_win) + ...
                       strcmp( ff{1} , 'Button'  )                 ) ;

       c.SingleHeight = ceil( f * ( 1 * c.PixelFit(1,2) + c.PixelFit(2,2) ) );

       c.SingleHeight = c.SingleHeight - c.PixelFit(2,2) * ...
                                strcmp( ff{1} , 'Edit' ) * is_win;

     %------------------------------------------
     case 'Tag'

       f = 1 + 0.2;

       c.SingleHeight = ceil( f * ( 1 * c.PixelFit(1,2) + c.PixelFit(2,2) ) );


     %------------------------------------------
     case 'Separator'

       c.VSpace =   vd;


   end

    cnf = setfield( cnf , ff{1} , c );

end

%***************************************************

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
%    A positive complex Value for DIM, to removes Blanks only from Start,
%    A negative complex Value for DIM, to removes Blanks only from End.
%       
%  CHAR specifies BlankCharacters to remove
%       default:  [ 160  32  13  10  9  0 ];  % [ NBSP Space CR LF TAB ZERO ]
%

  
msg = '';
 nl = char(10);

Nin = nargin;

if Nin < 1
  error('Not enough Input Arguments.')
else
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
    dim = dim(:);
    if ~all( ( abs(dim) == 1 ) |  ( abs(dim) == 2 ) )
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
             'Values for Input DIM must define 1. or 2. Dimension.' ];
    end
  end 
end

if Nin < 3
  cc = [ 160  32  13  10  9  0 ];  % [ NBSP  Space CR LF TAB ZERO ]
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
  str = double(str);
  return
end


  blank  = 0*str;

  for ii = cc
    blank = ( blank | ( str == ii ) );
  end

  si = size(str);

  for ii = 1 : size(dim,1)

    d = dim(ii);

    s = sign(imag(d));  % Remove from wich Side:  1  0  -1 
 
    d = abs(d);

    jj = find( sum(blank,3-d) == si(3-d) );  % Columns with full Blanks

    if ~isempty(jj) 

         p  = [ 3-d  d ];
        str = permute(str,p);

         jj = jj(:)';
         nb = size(jj,2);

        %--------------------------------------------
        % Blank at Begin

        ind = ( 1 : nb );
        jj1 = find( ( ( jj == ind ) & ( s >= 0 ) ) );

        %--------------------------------------------
        % Blank at End

        ind = ind + si(d) - nb;
        jj2 = find( ( ( jj == ind ) & ( s <= 0 ) ) );

        %--------------------------------------------

        str(:,jj([jj1 jj2])) = [];

        str = permute(str,p);

    end
    
  end

  str = char(str);

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [fs,ssi,ppi,is_win,fn] = bestfont

% BESTFONT  Returns an optimal FontSize for this computer
%
% [ FontSize , ScreenSize , ScreenPixelsPerInch ] = BESTFONT
%
% [ ... , IsWin , FixedWidthFontName ] = BESTFONT
%
% returns true for PCWIN-System and the FixedWidthFontname
%

is_win = strcmp( upper(computer) , 'PCWIN' );

uni = get(0,'units');       
      set(0,'units','pixels')
ssi = get(0,'ScreenSize');  
      set(0,'units',uni);
          
ppi = get(0,'ScreenPixelsPerInch');

is_tall = -1 + ( ssi(4) >=  480 ) + ...
               ( ssi(4) >=  600 ) + ...
             1*( ssi(4) >= 1050 );

fs =  8 + 2 * is_tall - 1 * is_win;


if is_win
   fn = 'courier';
else
   fn = { 'arrial'  get(0,'fixedwidthfontname') };
   fn = fn{ 1 + ( ssi(4) >= 1050 ) } ;
end


%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function wygiwys(fig)

% WYGIWYS  WhatYouGetIsWhatYouSee
%
% Switch the PaperPosition to the same View like on the Screen
% The Figure will positioned in the Center of the Paper
%
% WYGIWYS( FigureHandle )
%
%


if nargin < 1
  fig = get(0,'currentfigure');
end

if isempty(fig)
 return
end

ok = ( isnumeric(fig)  &  ( prod(size(fig)) == 1 ) );
if ok
   ok = ishandle(fig);
   if ok
      ok = strcmp( get(fig,'type') , 'figure' );
   end
end

if ~ok
   error('Input must be a FigureHandle.');
end


figuni = get(fig,'units');
papuni = get(fig,'paperunits');

set(fig,     'units' , 'pixels' , ...
        'paperunits' , 'inches'       );

figpos = get(fig,'position');
pappos = get(fig,'paperposition');
pap_si = get(fig,'papersize');

ppi    = get(0,'screenpixelsperinch');

pappos([3 4]) = figpos([3 4]) / ppi;
pappos([1 2]) = (pap_si-pappos([3 4])) / 2;

set(fig,'paperposition',pappos);

set(fig,     'units' , figuni   , ...
        'paperunits' , papuni         );

