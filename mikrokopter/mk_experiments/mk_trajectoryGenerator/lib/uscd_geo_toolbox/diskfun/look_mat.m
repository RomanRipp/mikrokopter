function [msg,varargout] = look_mat(varargin);

% LOOK_MAT  Display Informations of a MAT-File
%
% [Msg,CellText,VarStruct,FileStruct ] = LOOK_MAT( FILENAME )
%
%  Display Informations of MAT-File, specified by FILENAME.
%
% CellText: Contents of MatFile in a CellStringArray
%
%  VarStruct: name size bytes class
%
% FileStruct: name date bytes isdir creation
%
%  Msg contains the ErrorMessages.
%
%-----------------------------------------------------------------------
% DisplayMode
%
%  LOOK_MAT( ... , MODE ) , where Mode is a scalar, specifies the Display:
%
%  MODE == 0   displays the TEXT on the Terminal,
%  MODE ~= 0   displays the TEXT in a ListBox of a Figure.
%
%  In case of ( MODE ~= 0 ), the second Output is the Handle 
%   of the LOOK_MAT-Figure: [ Msg , FIG , ... ] = LOOK_MAT( ... )
%
%  Use the Value of FIG for MODE, to reset the LOOK_MAT-Figure with
%    the Informations about a new MAT-File.
%
%-----------------------------------------------------------------------
% FileMode
%
%  LOOK_MAT( ... , OUTFILENAME )  saves the TEXT into an ASCII-File,
%                                 specified by OUTFILENAME.
%
%  By default, a new File is created, an existing File with the same Name
%    will be overwritten. If the OUTFILENAME begins with an '@',
%    the TEXT will be append on an existing File.
%
%  The Inputs MODE and OUTFILENAME are optional 
%   and could be used together.
%
%-----------------------------------------------------------------------
%
%  See also: WHOSFILE
%

Nin  = nargin;
Nout = nargout;

NoOut = ( Nout == 0 );

Nout = Nout - 1;

msg = '';

varargout = cell(1,Nout);

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
msg0 = sprintf('%s: ',upper(fcn));



nl   = char(10);

%*************************************************************
% Check for CallBack

if Nin > 1
   fcn = varargin{1};
   if chkstr(fcn,1)
      if ( double(fcn(1)) == double('#') ) & ( size(fcn,2) > 1 )
         try
            feval(fcn(2:end),varargin{2:end});
         catch
            warning(sprintf('%s Error Callback: %s\n%s',msg0,fcn,lasterr));
         end
         return
      end
   end
end

%*************************************************************
% Check Inputs

[msg,file,mode,outfile] = checkin(varargin{:});

if ~isempty(msg)
    msg = sprintf('%s%s',msg0,msg);
    if NoOut
       fprintf(1,char(7));
       fprintf( 1 , sprintf('%s%s%s%s',nl,msg,nl,nl) );
       clear msg
    end
    return
end

if NoOut  &  isempty(outfile)  &  isempty(mode)
   mode = -1;    % Default: ListBoxWindow
end

%*************************************************************
% Get Info

is_mode = ~isempty(mode);
is_file = ~isempty(outfile);

[msg,txt,vlst,flst,id] = look( file );

if ~isempty(msg)
    msg = sprintf('%s%s',msg0,msg);
    if NoOut
       fprintf(1,char(7));
       fprintf( 1 , sprintf('%s%s%s%s',nl,msg,nl,nl) );
       clear msg
    end
    return
end

%*************************************************************
% Write File

if is_file
   msg = write(outfile,txt);
   if ~isempty(msg)
       msg = sprintf('%s%s',msg0,msg);
   end
end

%*************************************************************
% Display

fig = [];
if is_mode
   if mode == 0
   % Terminal
     fprintf(1,char(10));
     fprintf(1,'%s\n',txt{:});
   else
   % Figure with ListBox
     fig = list(file,txt,id,mode);
   end
end

%*************************************************************
% OutPut

if NoOut
   if ~isempty(msg)
       fprintf(1,char(7));
       fprintf( 1 , sprintf('%s%s%s%s',nl,msg,nl,nl) );
   end
   clear msg
end

if Nout <= 0
   return
end

if  ~isempty(fig)
   v = { fig txt vlst flst };
else
   v = { txt vlst flst };
end

n = min( Nout , prod(size(v)) );

varargout(1:n) = v(1:n);

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [msg,file,mode,outfile] = checkin(varargin)

% CHECKIN  First Check of Inputs

Nin = nargin;

msg     = '';
file    = '';
mode    = [];
outfile = '';

%-------------------------------------------------------------
% Check FileName

if Nin < 1
   file = '';
else
   file = varargin{1};
end

if isempty(file)
   [f,p] = uigetfile('*.mat','Select a MAT-File ...');
   if isequal(f,0)
      return
   end
   file = fullfile(p,f);
end
       
if ~chkstr(file,0)
    msg = 'Input FileName must be a String.';
    return
end

for e = { ''  '.mat' };
    f  = cat( 2 , file , e{1} );
    ok = ( exist(f) == 2 );
    if ok
       break
    end
end

if ok
   file = f;
else
   msg = sprintf('File "%s" doesn''t exist.',file);
   return
end

%-------------------------------------------------------------
% Check Options

for ii = 2 : Nin
    if chkstr(varargin{ii},1)
       outfile = varargin{ii};
    elseif isnumeric(varargin{ii}) & ( prod(size(varargin{ii})) == 1 )    
          mode = varargin{ii};
    end
end


%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,str,v,d,id] = look(file);

% LOOK (WHOSFILE)  returns Contents of MAT-File as CellString
%
% [Msg,CellText,VarStruct,FileStruct] = WHOSFILE( MatFileName )
%
%  VarStruct: name size bytes class
% FileStruct: name date bytes isdir
%
% see also: WHOS, LOOK_MAT
%

Msg  = '';

 str = cell(0,1);
  id = cell(0,1);

   v = [];
   d = [];

  nl = char(10);

%---------------------------------------------------------
% Check File

if ~( ischar(file) & ~isempty(file) & ...
      ( prod(size(file)) == size(file,2) ) );
   Msg = 'Input File must be a String.';
   return
end

ok = ( exist(file,'file') == 2 );
if ~ok
    f  = cat(2,file,'.mat');
    ok = ( exist(f,'file') == 2 );
    if ok
       file = f;
    end
end

if ~ok
    Msg = ['File '  file ' does not exist.'];
    return
end

f = which(file); 
if ~isempty(f)
    file = f;
end

%---------------------------------------------------------
% DIR

d = dir(file);

if isempty(d)
   Msg = [ 'Can''t read File '  file '.' ];
elseif ~( prod(size(d)) == 1 )
   Msg = [ 'Multiple read of File '  file '.' ];
elseif isequal(d.isdir,1)
   Msg = [ 'File '  file ' is a Directory.' ];
elseif isequal(d.bytes,0)
   Msg = [ 'File '  file ' is empty.' ];
end

if ~isempty(Msg)
    return
end

d.name = file;

%---------------------------------------------------------
% WHOS

try
  v = whos('-file',file);
catch 
  Msg = lasterr;
end

if ~isempty(Msg) | isempty(v)
  Msg = [ 'File '  file ' must be a MAT-File.' ...
            nl(1:(end*(~isempty(Msg))))  Msg ];
  return
end

%---------------------------------------------------------

bl = char(32);

nn = prod(size(v));
   
n = str2mat( 'Name'  , bl , v.name  , bl , 'Summary' );  % Name
c = str2mat( 'Class' , bl , v.class , bl , bl        ); % Class

s = cell(nn+4,1);   % Size
b = cell(nn+4,1);   % Bytes

s{   1} = cat( 2 , bl(1,ones(1,8)) , 'Size' );
s{   2} = bl;
s{nn+3} = bl;

b{   1} = cat( 2 , bl(1,ones(1,4)) , 'Bytes' );
b{   2} = bl;
b{nn+3} = bl;

ss = 0;
bb = 0;

for ii = 1 : nn

    s{ii+2} = sprintf('%9.0f',v(ii).size(1));

    for jj = 2 : size(v(ii).size,2);

       s{ii+2} = cat( 2 , s{ii+2} , ...
                     sprintf(' x %.0f',v(ii).size(jj)) );
   
    end

    b{ii+2} = sprintf('%9.0f',v(ii).bytes);

    bb = bb + v(ii).bytes;
    ss = ss + prod(v(ii).size);

    s{nn+4} = sprintf('%13.0f',ss);
    b{nn+4} = sprintf('%9.0f',bb);

end

n = cat( 2 ,      n  , bl(ones(nn+4,1),ones(1,2)) );
s = cat( 2 , char(s) , bl(ones(nn+4,1),ones(1,2)) );
b = cat( 2 , char(b) , bl(ones(nn+4,1),ones(1,3)) );


str = cat(2,bl(ones(nn+4,1),1),n,s,b,c);

str(   2,:) = '-';
str(nn+3,:) = '-';

str = cellstr(str);


id = cell(nn+4,1);
id(:) = {''};

id(2+(1:nn)) = {v.name}';

%--------------------------------------------------
% FileInfo first

siz = d.bytes;

scl = 10 * [ 0  1  2 ];
scl = 2 .^ scl;

ii = sum( siz >= scl );
ii = max(ii,1);

pre = { ''  'k'  'M' };
fmt = sprintf('%s.%.0ff','%',3*(ii>1));
fmt = sprintf(' %s %sBytes',fmt,pre{ii});

siz = sprintf(fmt,siz/scl(ii));

dat = info(file);  % Read CreationDate

txt = cat( 1 , {sprintf('File: %s',d.name)} , ...
               {sprintf('Size: %s',siz)} , ...
               {sprintf('Date: %s',dat)}       );

txt = txt( 1 : end-isempty(dat) );

str = cat( 1 , {''} , txt , {char('-'*ones(1,size(char(txt),2)))} , ...
               {''} , str , {''} );

ii = cat( 2 , ones(1,size(txt,1)+3) , (1:size(id,1)) , 1 );

id = id(ii);

d.creation = dat;

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function dat = info(file);

% Read CreationDate from MAT-File

dat = '';

fid = fopen(file,'r');

if fid == -1
   return
end

b = fread(fid,[1 124],'char');

fclose(fid);

bd = ~( (  32 < b  &   b <= 126 ) | ...
        ( 160 < b  &   b <= 255 )        );

if all(bd) | ( size(b,2) < 124 )
   return
end

bd = find(bd);

b(bd) = 32;

b = char(b);

str = 'Created on:';

ii = findstr(b,str);

if ~( prod(size(ii)) == 1 )
      return
end

dat = rmblank( b( ii+size(str,2) : end ) , 2 );

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function msg = write(file,txt)

% WRITE  Write Text to OutFile
 

msg = '';

 if file(1) == '@'
    file = file(2:end);
    mode = 'a';
 else
    mode = 'wt';
 end

 fid  = fopen(file,mode);

  if fid ~= -1
     fprintf( fid , '%s\n' , txt{:} );
     fclose(fid);
  else
     msg = sprintf('Error open OutFile "%s" for writing.',file);
  end

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  fig = list(file,txt,id,fig)

% LIST  Display Info in ListBoxFigure

Nin = nargin;

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end

tag = upper(fcn);

new = ( Nin < 4 );
if ~new
    new = ~chkhndl(fig,'figure',tag);
end

if new
   fig = newfig(tag,fcn);
end

set( fig , 'name' , file );

hl = findobj(fig,'type','uicontrol','tag','LIST');

set( hl , 'listboxtop' , 1   , ...
               'value' , 1   , ...
              'string' , txt , ...
            'userdata' , id        );

hc = findobj(fig,'type','uicontextmenu','tag','CONTEXT');

set( hc , 'visible' , 'off' , 'userdata' , '' );

ha = findobj(fig,'type','uimenu','tag','ASSIGN');

set( ha , 'userdata' , file      , ...
          'label'    , '&Assign' , ...
          'visible'  , 'off'           );

drawnow

refresh(fig);

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function fig = newfig(tag,fcn);

% NEWFIG   Create New ListBoxFigure


  %--------------------------------------
  scr_uni = get(0,'units');      set(0,'units','pixels')
  scr_si  = get(0,'ScreenSize'); set(0,'units',scr_uni);
  
  ppi    = get(0,'ScreenPixelsPerInch');
        
  %--------------------------------------

   is_tall = -1 + ( scr_si(4) >=  480 ) + ...
                  ( scr_si(4) >=  600 ) + ...
                  ( scr_si(4) >= 1024 );

   is_win = strcmp( upper(computer) , 'PCWIN' );

   fontsize =  8 + 2 * is_tall - 2 * is_win;

 
   if is_win
      fontname = 'courier';
   else
      fontname = { 'arrial'  get(0,'fixedwidthfontname') };
      fontname = fontname{ 1 + ( scr_si(4) >= 1050 ) } ;
   end

  %--------------------------------------
  
  fs = floor( 25/18 * ppi/100 * fontsize );  % Points --> Pixels

  ww = 80;                   % Character
  hh = 50;
                          
  fig00 = ceil([ 0.55*ww  1.2*hh ] * fs);                         
  fig11 = floor([ 1/2  2/3 ].*scr_si(3:4));
  
  figpos = NaN*ones(1,4);
  
  figpos(3:4)= fig00 + ( fig11 - fig00 ) .* ( fig11 < fig00 );
  
  voffs = max( 60 , min( ceil(1/6*scr_si(4)) , 80 ) );

  figpos(1) = 50;
  figpos(2) = scr_si(4)-voffs-figpos(4);
  

 fig  = figure('position'   , figpos , ...
               'numbertitle', 'off'  , ...
             'integerhandle', 'off'  , ...
               'menubar'    , 'none' , ...
               'toolbar'    , 'none' , ...
               'name'       , tag    , ...
               'tag'        , tag    , ...
               'createfcn'  , ''     , ...
          'handlevisibility','callback' );

    ParentCB = 'get(gcbo,''parent'')';

       NewCB = sprintf('%s(''#new'',%s);',fcn,ParentCB);
      DownCB = sprintf('%s(''#down'',gcbo);',fcn);
    SelectCB = sprintf('%s(''#select'',gcbo);',fcn);
    AssignCB = sprintf('%s(''#assign'',gcbo);',fcn);
   RefreshCB = sprintf('refresh(%s);',ParentCB);
     CloseCB = sprintf('delete(%s);',ParentCB);

  hx = uicontextmenu( 'parent'  , fig , ...
                      'tag'     , 'CONTEXT' , ...
                     'callback' , SelectCB , ...
                      'visible' , 'off' );

  ha = uimenu('parent' , hx           , ...
              'label'  , '&Assign'    , ...
            'callback' ,  AssignCB    , ...
                 'tag' , 'ASSIGN'     , ...
             'visible' , 'off'       );

  hl = uicontrol( 'parent'    , fig        , ...
                  'style'     , 'listbox'  , ...
            'backgroundcolor' , [1 1 1]    , ...
            'foregroundcolor' , [0 0 0]    , ...
                  'units'     ,'normalized', ...
                  'position'  , [0 0 1 1] , ...
                  'min'       , 0         , ...
                  'max'       , 1         , ...
                  'fontunits' , 'points'  , ...
                  'fontsize'  , fontsize  , ...
                  'fontname'  , fontname  , ...
                  'string'    , ''        , ...
                  'tag'       , 'LIST'    , ...
               'callback'     , DownCB    , ...
              'uicontextmenu' , hx        , ... 
        'horizontalalignment' , 'left'           );

   hn = uimenu( 'parent'      , fig   , ...
                'label'       , '&New' , ... 
                'tag'         , 'NEW' , ...
                'callback'    , NewCB       );

   hr = uimenu( 'parent'      , fig      , ...
                'label'       , '&Refresh' , ...
                'tag'         , 'REFRESH' , ...
                'callback'    , RefreshCB     );
 
   hh = uimenu( 'parent'      , fig     , ...
                'label'       , '&Help' , ...
                'tag'         , 'HELP' , ...
                'callback'    , ''           );

        uimenu( 'parent' , hh , ...
                'label'  , '&Click right to assign selected Variable' , ...
                'tag'    , 'Help01' , ...
              'callback' , ''           );

   hb = uimenu( 'parent'      , fig     , ...
                'label'       , '   '   , ...
                'tag'         , 'BLANK' , ...
                'callback'    , ''           );

   hc = uimenu( 'parent'      , fig     , ...
                'label'       , '&Close' , ...
                'tag'         , 'CLOSE' , ...
                'callback'    , CloseCB       );

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function new(fig)

% CallBack for New-Menu

   [f,p] = uigetfile('*.mat','Select a MAT-File ...');
   if isequal(f,0)
      return
   end

   file = fullfile(p,f);

   msg = look_mat(file,fig);

   if ~isempty(msg)
       msg = cat( 2 , 'Error using LOOK_MAT' , char(10) , msg );
       warndlg(msg,'Error','warn');
       return
   end

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function down(hl)

fig = get(hl,'parent');

sel = get(fig,'selectiontype');

ok = strcmp(sel,'normal');

if ok

   id = get(hl,'userdata');
   sz = size(id);

   vl = get(hl,'value');

   ok = ( chkcstr(id,1) & isequal(sz,size(get(hl,'string'))) );
   if ok
      ok = ( vl <= prod(sz) );
      if ok
         ok = ~isempty(id{vl});
      end
   end

end

hc = findobj( fig , 'type' , 'uicontextmenu' , 'tag' , 'CONTEXT' );
ha = get( hc , 'children' );

if ~ok
    set( hl , 'uicontextmenu' , [] );
    set( hc , 'visible' , 'off' , 'userdata' , '' );
    set( ha , 'visible' , 'off' , 'label' , '&Assign' ); 
    return
end

set( hc , 'userdata' , id{vl} );

set( ha , 'label' , sprintf('&Assign   %s  ',id{vl}) , ...
        'visible' , 'on' );

set( hl , 'uicontextmenu' , hc );

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function select(hc)

% SELECT  CallBack for Context-Menu

   ha = get( hc , 'children' );

   file = get( ha , 'userdata' );
   name = get( hc , 'userdata' );

   if strcmp(get(ha,'visible'),'off') | ...
      isempty(file) | isempty(name)
      set(hc,'visible','off');
   end

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function assign(ha)

% ASSIGN  CallBack for Assign-Menu

   hc = get( ha , 'parent' );

   file = get( ha , 'userdata' );
   name = get( hc , 'userdata' );

   if strcmp(get(ha,'visible'),'off') | ...
      isempty(file) | isempty(name)
      return
   end

   msg = '';

   if ~chkstr(name,1)
       msg = 'Invalid UserData of UIMenu.';
   else
       try
          v = load(file,name);
       catch
          msg = sprintf('Error call LOAD( %s , %s ).\n%s',file,name,lasterr);
       end 
   end

   if isempty(msg)
      if ~isequal(fieldnames(v),{name})
          msg = sprintf('Didn''t found Variable "%s" in File "%s"',name,file);
      else
          try
             v = getfield(v,name);
             assignin('base',name,v);
          catch
             msg = lasterr;
          end
      end
   end
  
   if ~isempty(msg)
       msg = sprintf('Error Assign Variable %s\n%s',name,msg);
       warndlg(msg,'Error','warn');
       return
   end


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
%    A positive complex Part of DIM, to removes Blanks only from Start,
%    A negative complex Part of DIM, to removes Blanks only from End.
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
    if ~all( ( real(dim) == 1 ) |  ( real(dim) == 2 ) )
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
 
    d = real(d);

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

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
