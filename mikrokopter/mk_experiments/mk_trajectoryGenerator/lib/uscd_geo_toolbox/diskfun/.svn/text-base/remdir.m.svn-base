function [Msg,ok] = remdir(pfad,force);

% REMDIR  Removes Directory 
%
% [Msg,Ok] = REMDIR( Pfad , force )
%
% Removes recursivly Directory with name Pfad.
%
%  force == 0   Question, if Directory not empty (default)
%           1    forced removing
%
% Ok = -1  removing  Cancelled
% Ok =  0  Directory not removed
% Ok =  1  Directory removed
%

Msg = '';
ok  =  0;
 
nl = char(10);

if nargin < 1
  Msg = 'Input Pfad is missing.';
  return
end

if ~( ischar(pfad)  &  ( prod(size(pfad)) == size(pfad,2) )  )
  Msg = 'Input Pfad must be a String.';
  return
end


if ~( exist(pfad,'dir') == 7 )
  Msg = 'Input Pfad must be a DirectoryName.';
  return
end


if nargin < 2
  force = 0;
end


force = isequal(force,1);


%------------------------------------

fs = filesep;

is_fs = strcmp(pfad(end),fs);
 

if strcmp( getenv('HOME') , pfad( 1 : (end-is_fs) ) );
   Msg = cat( 2 , 'Can''t remove HOME-Directory: ' , pfad );
   return
end


pfad = cat( 2 , pfad , fs(1:(end*(~is_fs))) );


%------------------------------------
% Get DirectoryContents

d = dir(pfad);

if isempty(d)

   c = cell(0,1);
   d = zeros(0,1);

else

   c = { d.name };
   d = cat( 1 , d.isdir );
     
   bad = find( strcmp( c , '.' ) | strcmp( c , '..' ) );

   c(bad) = [];
   d(bad) = [];

end

%------------------------------------
% Empty Directory

if isempty(c)

  Msg = rm_dir(pfad);
  
  return

end


%------------------------------------
% Question

if ~force

   txt = [ 'Directory: ' pfad nl nl ...
           ' is not empty! Remove anyway?' nl nl ];

   ok = questdlg(txt,'Stop','Remove','Cancel','Cancel');

   if ~strcmp(ok,'Remove')
      ok = -1;
      return
   end

   if any(d)

     txt = [ 'Directory: ' pfad nl nl ...
             ' contains SubDirectories! Remove anyway?' nl nl ];

     ok = questdlg(txt,'Stop','Remove','Cancel','Cancel');

     if ~strcmp(ok,'Remove')
        ok = -1;
        return
     end

   end

end

%------------------------------------
% Remove Directory recursivly

ii = find(d);

for jj = ii(:)'

    pf = cat( 2 , pfad , c{jj} );

    %-----------------------------
    % Check for symbolic Links

    ok = ~isunix;

    if ~ok
        ok = ~islink(pf);
    end

    if ok

      fprintf(pf);
      fprintf(nl);

      msg = remdir( pf , 1 );

      Msg = [ Msg nl(1:(end*(~isempty(Msg)))) msg ];

    else

      d(jj) = 0;

    end

end

%------------------------------------
% Files

ii = find(~d);

for jj = ii(:)'

   file = cat(2,pfad,c{jj});

   msg = '';

   try
     delete(file)
   catch 
     msg = lasterr;
   end

   if ~isempty(msg)
  
       Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
               'Error call DELETE( ' file ' ).' nl msg nl ];

   elseif ~isempty(dir(file))

       Msg = [ Msg  nl(1:(end*(~isempty(Msg)))) ...
               'Can''t delete File:'  file ];

   end

end
     

%------------------------------------
% Get DirectoryContents again

d = dir(pfad);

if isempty(d)

   c = cell(0,1);
   d = zeros(0,1);

else

   c = { d.name };
   d = cat( 1 , d.isdir );
     
   bad = find( strcmp( c , '.' ) | strcmp( c , '..' ) );

   c(bad) = [];
   d(bad) = [];

end

%------------------------------------
% Empty Directory

if isempty(c)

  msg = rm_dir(pfad);
  
else

  msg = [ 'Directory: '  pfad  '  is not empty.' ];

end

ok = ~( exist(pfad,'dir') == 7 );

Msg = [ Msg nl(1:(end*(~isempty(Msg)))) msg ];



%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function Msg = rm_dir(pfad);

nl = char(10);

Msg = '';

%-------------------------------------------------------
% for Windows !!!

pfad = pfad( 1 : ( end - strcmp(pfad(end),filesep) ) );

%-------------------------------------------------------


command = cat(2,'remdir ',pfad);
 
if isunix

    [status,Msg] = unix(command);

elseif strcmp(computer,'PCWIN') 

    [status,Msg] =  dos(command);

else

    command = [ '! ' command ];

    eval(command,'status = -1; Msg = lasterr;')

end

if isempty(Msg)  &  ~isempty(dir(pfad))

  Msg = [  'Can''t remove Directory:'  nl pfad ...
           nl(1:(end*(~isempty(Msg))))   Msg    ];

end  
 

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,typ,link,name] = islink(name);

% ISLINK  Checks if argument is a symbolic Link
%
% [ IsLink , Type , Value , Target ] = ISLINK( Name )
%
%  IsLink:  1 if Name is a symbolic link by UNIX: test -L,
%           0 if Name is not a symboliclink or on non-UNIX systems
%          [] if Name is empty
%
%  Type     2 if Name or the target of the link is a file,
%           7 if Name or the taget of the link is a directory
%           0 if Name is not a File or Directory
%    
%  Value:  Value of symbolic Link (Derefer) by UNIX: readlink -m
%
%  Target: Target of symbolic Link by UNIX: ls -ld; or origin File or Directory
%
%
% uses UNIX commands: test, readlink, ls
%

Nout = nargout;

ok     = NaN;
typ    = 0;
link   = '';

if ischar(name) & isempty(name)
   ok = []; 
   return
elseif ~( ischar(link) &  ( size(typ,2) == prod(size(typ)) ) )
   error('Input must be a String.');
end

%------------------------------------------------------
% Check if Directory or File exist

typ = 2 * ( exist(name,'file') == 2 ) + ...
      7 * ( exist(name,'dir')  == 7 );

ok = ~( typ == 0 );

if ok
   link = name;
end

ok = ( ok & isunix );

if ~ok
    return
end


%------------------------------------------------------
% Remove FileSeparator from End in case of Directory !!!

name = name( 1 : end-( strcmp( link(end) , filesep ) & ( typ == 7 ) ) );

%------------------------------------------------------
% Check for Link, use UNIX: test

[s,w] = unix(sprintf('test -L "%s"',name));

ok = ( s == 0 );

if ~ok | ( Nout < 3 )
    return
end

%------------------------------------------------------
% Derefer Link, use UNIX: readlink

[s,link] = unix(sprintf('readlink -n -m "%s"',name));

if ~( s == 0 )
    link = '';
end

if ( Nout < 4 )
    return
end

%------------------------------------------------------
% Get Target of Link, use UNIX: ls

pre = '->';

[s,l] = unix(sprintf('ls -ld "%s"',name));

ii = findstr( l , pre );

if isempty(ii)
   return
end

ii = max(ii) + size(pre,2);

name = rmblank(l(ii:end),2);

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ButtonName=questdlg(Question,Title,Btn1,Btn2,Btn3,Default)
%QUESTDLG Question dialog box.
%  ButtonName=QUESTDLG(Question) creates a modal dialog box that 
%  automatically wraps the cell array or string (vector or matrix) 
%  Question to fit an appropriately sized window.  The name of the 
%  button that is pressed is returned in ButtonName.  The Title of 
%  the figure may be specified by adding a second string argument.  
%  Question will be interpreted as a normal string.  
%
%  The default set of buttons names for QUESTDLG are 'Yes','No' and 
%  'Cancel'.  The default answer for the above calling syntax is 'Yes'.  
%  This can be changed by adding a third argument which specifies the 
%  default Button.  i.e. ButtonName=questdlg(Question,Title,'No').
%
%  Up to 3 custom button names may be specified by entering
%  the button string name(s) as additional arguments to the function 
%  call.  If custom ButtonName's are entered, the default ButtonName
%  must be specified by adding an extra argument DEFAULT, i.e.
%
%    ButtonName=questdlg(Question,Title,Btn1,Btn2,DEFAULT);
%
%  where DEFAULT=Btn1.  This makes Btn1 the default answer.
%
%  To use TeX interpretation for the Question string, a data
%  structure must be used for the last argument, i.e.
%
%    ButtonName=questdlg(Question,Title,Btn1,Btn2,OPTIONS);
%
%  The OPTIONS structure must include the fields Default and Interpreter.  
%  Interpreter may be 'none' or 'tex' and Default is the default button
%  name to be used.
%
%  A sample application of this function is:
%
%    ButtonName=questdlg('What is your wish?', ...
%                        'Genie Question', ...
%                        'Food','Clothing','Money','Money');
%
%  
%    switch ButtonName,
%       case 'Food', 
%        disp('Food is delivered');
%      case 'Clothing',
%        disp('The Emperor''s  new clothes have arrived.')
%      case 'Money',
%        disp('A ton of money falls out the sky.');
%    end % switch
%
%  See also TEXTWRAP, INPUTDLG.

%  Author: L. Dean
%  Copyright (c) 1984-98 by The MathWorks, Inc.
%  $Revision: 5.37 $

if nargin<1,error('Too few arguments for QUESTDLG');end

Interpreter='none';
if ~iscell(Question),Question=cellstr(Question);end

if strcmp(Question{1},'#FigKeyPressFcn'),
  QuestFig=get(0,'CurrentFigure');
  AsciiVal= abs(get(QuestFig,'CurrentCharacter'));
  if ~isempty(AsciiVal),
    if AsciiVal==32 | AsciiVal==13,
      set(QuestFig,'UserData',1);
      uiresume(QuestFig);
    end %if AsciiVal
  end %if ~isempty
  return
end
%%%%%%%%%%%%%%%%%%%%%
%%% General Info. %%%
%%%%%%%%%%%%%%%%%%%%%
Black      =[0       0        0      ]/255;
LightGray  =[192     192      192    ]/255;
LightGray2 =[160     160      164    ]/255;
MediumGray =[128     128      128    ]/255;
White      =[255     255      255    ]/255;

%%%%%%%%%%%%%%%%%%%%
%%% Nargin Check %%%
%%%%%%%%%%%%%%%%%%%%
if nargout>1,error('Wrong number of output arguments for QUESTDLG');end
if nargin==1,Title=' ';end
if nargin<=2, Default='Yes';end
if nargin==3, Default=Btn1;end
if nargin<=3, Btn1='Yes'; Btn2='No'; Btn3='Cancel';NumButtons=3;end
if nargin==4, Default=Btn2;Btn2=[];Btn3=[];NumButtons=1;end
if nargin==5, Default=Btn3;Btn3=[];NumButtons=2;end
if nargin==6, NumButtons=3;end
if nargin>6, error('Too many input arguments');NumButtons=3;end

if isstruct(Default),
  Interpreter=Default.Interpreter;
  Default=Default.Default;
end


%%%%%%%%%%%%%%%%%%%%%%%
%%% Create QuestFig %%%
%%%%%%%%%%%%%%%%%%%%%%%
FigPos=get(0,'DefaultFigurePosition');
FigWidth=75;FigHeight=45;
FigPos(3:4)=[FigWidth FigHeight];
QuestFig=dialog(                                               ...
               'Visible'         ,'off'                      , ...
               'Name'            ,Title                      , ...
               'Pointer'         ,'arrow'                    , ...
               'Units'           ,'points'                   , ...
               'Position'        ,FigPos                     , ...
               'KeyPressFcn'     ,'questdlg #FigKeyPressFcn;', ...
               'UserData'        ,0                          , ...
               'IntegerHandle'   ,'off'                      , ...
               'WindowStyle'     ,'normal'                   , ... 
               'HandleVisibility','callback'                 , ...
               'Tag'             ,Title                        ...
               );

%%%%%%%%%%%%%%%%%%%%%
%%% Set Positions %%%
%%%%%%%%%%%%%%%%%%%%%
DefOffset=3;

IconWidth=32;
IconHeight=32;
IconXOffset=DefOffset;
IconYOffset=FigHeight-DefOffset-IconHeight;
IconCMap=[Black;get(QuestFig,'Color')];

DefBtnWidth=40;
BtnHeight=20;
BtnYOffset=DefOffset;
BtnFontSize= bestfont + 2;

BtnWidth=DefBtnWidth;

ExtControl=uicontrol(QuestFig   , ...
                     'Style'    ,'pushbutton', ...
                     'String'   ,' '         , ...
                     'FontUnits','points'   , ...                     
                     'FontSize' ,BtnFontSize , ...
                     'FontWeight' , 'bold'   ...
                     );
                     
for lp=1:NumButtons,
  eval(['ExtBtnString=Btn' num2str(lp) ';']);
  set(ExtControl,'String',ExtBtnString);
  BtnExtent=get(ExtControl,'Extent');
  BtnWidth=max(BtnWidth,BtnExtent(3)+8);
end % lp
delete(ExtControl);

MsgTxtXOffset=IconXOffset+IconWidth;

FigWidth=max(FigWidth,MsgTxtXOffset+NumButtons*(BtnWidth+2*DefOffset));
FigPos(3)=FigWidth;
set(QuestFig,'Position',FigPos);

BtnXOffset=zeros(NumButtons,1);

if NumButtons==1,
  BtnXOffset=(FigWidth-BtnWidth)/2;
elseif NumButtons==2,
  BtnXOffset=[MsgTxtXOffset
              FigWidth-DefOffset-BtnWidth];
elseif NumButtons==3,
  BtnXOffset=[MsgTxtXOffset
              0
              FigWidth-DefOffset-BtnWidth];
  BtnXOffset(2)=(BtnXOffset(1)+BtnXOffset(3))/2;
end

MsgTxtYOffset=DefOffset+BtnYOffset+BtnHeight;
MsgTxtWidth=FigWidth-DefOffset-MsgTxtXOffset-IconWidth;
MsgTxtHeight=FigHeight-DefOffset-MsgTxtYOffset;
MsgTxtForeClr=Black;
MsgTxtBackClr=get(QuestFig,'Color');

CBString='uiresume(gcf)';
for lp=1:NumButtons,
  eval(['ButtonString=Btn',num2str(lp),';']);
  ButtonTag=['Btn' num2str(lp)];
  
  BtnHandle(lp)=uicontrol(QuestFig            , ...
                         'Style'              ,'pushbutton', ...
                         'Units'              ,'points'    , ...
                         'Position'           ,[ BtnXOffset(lp) BtnYOffset  ...
                                                 BtnWidth       BtnHeight   ...
                                               ]           , ...
                         'CallBack'           ,CBString    , ...
                         'String'             ,ButtonString, ...
                         'HorizontalAlignment','center'    , ...
                         'FontUnits'          ,'points'    , ...
                         'FontSize'           ,BtnFontSize , ...
                         'FontWeight'         , 'bold'  , ...
                         'Tag'                ,ButtonTag     ...
                         );
                                   
end

MsgHandle=uicontrol(QuestFig            , ...
                   'Style'              ,'text'         , ...
                   'Units'              ,'points'       , ...
                   'Position'           ,[MsgTxtXOffset      ...
                                          MsgTxtYOffset      ...
                                          0.95*MsgTxtWidth   ...
                                          MsgTxtHeight       ...
                                         ]              , ...
                   'String'             ,{' '}          , ...
                   'Tag'                ,'Question'     , ...
                   'HorizontalAlignment','left'         , ...    
                   'FontUnits'          ,'points'       , ...
                   'FontWeight'         ,'bold'         , ...
                   'FontSize'           ,BtnFontSize    , ...
                   'FontWeight'         , 'bold'        , ...
                   'BackgroundColor'    ,MsgTxtBackClr  , ...
                   'ForegroundColor'    ,MsgTxtForeClr    ...
                   );

[WrapString,NewMsgTxtPos]=textwrap(MsgHandle,Question,75);

NumLines=size(WrapString,1);

% The +2 is to add some slop for the border of the control.
MsgTxtWidth=max(MsgTxtWidth,NewMsgTxtPos(3)+2);
MsgTxtHeight=NewMsgTxtPos(4)+2;

MsgTxtXOffset=IconXOffset+IconWidth+DefOffset;
FigWidth=max(NumButtons*(BtnWidth+DefOffset)+DefOffset, ...
             MsgTxtXOffset+MsgTxtWidth+DefOffset);

        
% Center Vertically around icon  
if IconHeight>MsgTxtHeight,
  IconYOffset=BtnYOffset+BtnHeight+DefOffset;
  MsgTxtYOffset=IconYOffset+(IconHeight-MsgTxtHeight)/2;
  FigHeight=IconYOffset+IconHeight+DefOffset;    
% center around text    
else,
  MsgTxtYOffset=BtnYOffset+BtnHeight+DefOffset;
  IconYOffset=MsgTxtYOffset+(MsgTxtHeight-IconHeight)/2;
  FigHeight=MsgTxtYOffset+MsgTxtHeight+DefOffset;    
end    

if NumButtons==1,
  BtnXOffset=(FigWidth-BtnWidth)/2;
elseif NumButtons==2,
  BtnXOffset=[(FigWidth-DefOffset)/2-BtnWidth
              (FigWidth+DefOffset)/2      
              ];
          
elseif NumButtons==3,
  BtnXOffset(2)=(FigWidth-BtnWidth)/2;
  BtnXOffset=[BtnXOffset(2)-DefOffset-BtnWidth
              BtnXOffset(2)
              BtnXOffset(2)+BtnWidth+DefOffset
             ];              
end

ScreenUnits=get(0,'Units');
set(0,'Units','points');
ScreenSize=get(0,'ScreenSize');
set(0,'Units',ScreenUnits);

FigPos(1)=(ScreenSize(3)-FigWidth)/2;
FigPos(2)=(ScreenSize(4)-FigHeight)/2;
FigPos(3:4)=[FigWidth FigHeight];

set(QuestFig ,'Position',FigPos);

BtnPos=get(BtnHandle,{'Position'});BtnPos=cat(1,BtnPos{:});
BtnPos(:,1)=BtnXOffset;
BtnPos=num2cell(BtnPos,2);  
set(BtnHandle,{'Position'},BtnPos);  

delete(MsgHandle);
AxesHandle=axes('Parent',QuestFig,'Position',[0 0 1 1],'Visible','off');

MsgHandle=text( ...
    'Parent'              ,AxesHandle                      , ...
    'Units'               ,'points'                        , ...
    'FontUnits'          ,'points'       , ...
    'FontSize'           ,BtnFontSize    , ...
    'FontWeight'         , 'bold'        , ...
    'HorizontalAlignment' ,'left'                          , ...
    'VerticalAlignment'   ,'bottom'                        , ...
    'HandleVisibility'    ,'callback'                      , ...
    'Position'            ,[MsgTxtXOffset MsgTxtYOffset 0], ...
    'String'              ,WrapString                      , ...
    'Interpreter'         ,Interpreter                     , ...
    'Tag'                 ,'Question'                        ...
    );

IconAxes=axes(                                      ...
             'Units'       ,'points'              , ...
             'Parent'      ,QuestFig              , ...  
             'Position'    ,[IconXOffset IconYOffset  ...
                             IconWidth IconHeight], ...
             'NextPlot'    ,'replace'             , ...
             'Tag'         ,'IconAxes'              ...
             );         
 
set(QuestFig ,'NextPlot','add');

IconData= ...
[2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2;
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2;
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2;
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2;
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 2 2 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 0 0 0 2 2 2 2 0 0 0 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 0 0 0 2 2 2 2 0 0 0 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2; 
 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];

Img=image('CData',IconData,'Parent',IconAxes);
set(QuestFig, 'Colormap', IconCMap);
set(IconAxes, ...
   'Visible','off'           , ...
   'YDir'   ,'reverse'       , ...
   'XLim'   ,get(Img,'XData'), ...
   'YLim'   ,get(Img,'YData')  ...
   );
set(findobj(QuestFig),'HandleVisibility','callback');
set(QuestFig ,'WindowStyle','modal','Visible','on');
drawnow;

uiwait(QuestFig);

TempHide=get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on');

if any(get(0,'Children')==QuestFig),
  if get(QuestFig,'UserData'),
    ButtonName=Default;
  else,
    ButtonName=get(get(QuestFig,'CurrentObject'),'String');
  end
  delete(QuestFig);
else
  ButtonName=Default;
end

set(0,'ShowHiddenHandles',TempHide);


%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [fs,ssi,ppi,is_win] = bestfont

% BESTFONT  Returns an optimal FontSize for this computer
%
% [ FontSize , ScreenSize , ScreenPixelsPerInch ] = BESTFONT
%

is_win = strcmp( upper(computer) , 'PCWIN' );

uni = get(0,'units');       
      set(0,'units','pixels')
ssi = get(0,'ScreenSize');  
      set(0,'units',uni);
          
ppi = get(0,'ScreenPixelsPerInch');

is_tall = -1 + ( ssi(4) >=  480 ) + ...
               ( ssi(4) >=  600 ) + ...
             0*( ssi(4) >= 1024 );

fs =  8 + 2 * is_tall - 1 * is_win;
