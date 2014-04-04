function [msg,varargout] = look_hdf(varargin);

% LOOK_HDF  Reads Basic Informations from a HDF multifile scientific Dataset-File
%
% [ Msg , DIM , VAR , GlobalATT , TXT ] = LOOK_HDF( FILENAME )
%
%  Reads Informations from HDF-File, specified by FILENAME.
%
%  VAR : Information for Variables,  [ Nvar by 13 ] CellArray
%
%  VAR =   1. Name
%          2. Rank
%          3. Type
%          4. Ndim
%          5. Size
%          6. Natt 
%          7. Range       : valid_range = [ Min Max ]
%          8. Fill        : _FillValue
%          9. DIM         : { Name  Length  Type  Scale  ATT    }
%         10. ATT         : { Name  Count   Type  Value  String }
%         11. Compressed  : 0 | 1
%         12. Chunked     : 0 | 1
%         13. ChunkLenght
%
%  DIM : Information for Dimensions, [ Ndim by  5 ] CellArray
%
%  DIM = { Name  Length  Type  Scale  ATT   }
%
%  ATT : Information for Attributes, [ Natt by  5 ] CellArray
%
%  ATT = { Name  Count   Type  Value  String }
%
%
%  TXT is an CellArray of Strings,
%        containing the Informations as ASCII-Text.
%
%  Msg contains the ErrorMessages.
%
%-----------------------------------------------------------------------
% DisplayMode
%
%  LOOK_HDF( ... , MODE ) , where Mode is a scalar, specifies the Display:
%
%  MODE == 0   displays the TEXT on the Terminal,
%  MODE ~= 0   displays the TEXT in a ListBox of a Figure.
%
%  In case of ( MODE ~= 0 ), the second Output is the Handle 
%   of the LOOK_HDF-Figure: [ Msg , FIG , ... ] = LOOK_HDF( ... )
%
%  Use the Value of FIG for MODE, to reset the LOOK_HDF-Figure with
%    the Informations about a new HDF-File.
%
%-----------------------------------------------------------------------
% FileMode
%
%  LOOK_HDF( ... , OUTFILENAME )  saves the TEXT into an ASCII-File,
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
%  See also: READ_HDF, ASSIGN_HDF, HDFSD, LOOK_CDF
%
%-----------------------------------------------------------------------
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

is_text  = ( is_mode | is_file | ( Nout > 3 ) );

[msg,dim,var,att] = look( file , max(Nout,3*is_text) );

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
% Get Text

if is_text
   txt = gettext(file,dim,var,att);
else
   txt = cell(0,1);
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
     fprintf(1,strhcat(txt,char(10)));
   else
   % Figure with ListBox
     fig = list(file,dim,var,txt,mode);
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
   v = { fig dim var att txt };
else
   v = { dim var att txt };
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
   [f,p] = uigetfile('*.nc','Select a HDF-File ...');
   if isequal(f,0)
      return
   end
   file = fullfile(p,f);
end
       
if ~chkstr(file,0)
    msg = 'Input FileName must be a String.';
    return
end

if ~( exist(file,'file') == 2 );
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

function [msg,dim,var,att] = look(file,mode);

% LOOK  Get Information about VAR, ATT

msg = '';
dim = cell(0,5);
var = cell(0,10);
att = cell(0,1);

%-------------------------------------------------------------
% Open HDF file

fid = hdfsd('start',file,'read');

if fid == -1
  ff = which(file);
  fid = hdfsd('start',ff,'read');
  if fid == -1
     msg = sprintf('Error open File "%s" as HDF-File.',file);
     return
  end
end

if mode == 0
   status = hdfsd('end',fid);
   return
end

%-------------------------------------------------------------
% HDF-Inquire

[nvar,natt,status] = hdfsd('fileinfo',fid);

if status == -1
   msg = 'FILEINFO: status = -1 '; 
   return
end


%--------------------------------------------------------------
% Global Attributes

att = getatt(fid,natt);

%--------------------------------------------------------------
% get Information about Variables
%

% var = { 1.Name 2.Rank 3.Type 4.Ndim 5.[Size] 6.Natt 7.Fill 8.[Range] 9.{dim} 10.{att} ...
%         11.Compressed 12. Chunked  13.ChunkLenght }

var = cell(nvar,10);

var(:,7) = {NaN};         % Fill
var(:,8) = {[NaN NaN]};   % Range

var(:,[11 12]) = {0};      % Compressed Chunked

for ii = 1 : nvar

    vid = hdfsd('select',fid,ii-1);

    % Name     Rank      Size      Type      Natt 

    [var{ii,1},var{ii,2},var{ii,5},var{ii,3},var{ii,6},status] = hdfsd('getinfo',vid);

    var{ii, 4} = size(var{ii,5},2);  % NDim

    var{ii, 9} = getdim(vid,var{ii,4});
    var{ii,10} = getatt(vid,var{ii,6});

    %---------------------------------------------------
    % Compressed / Chunked

    try
       [cl,ch,cm,status] = hdfsd('getchunkinfo',vid);
    catch
       status = -1;
    end

    if ~( status == -1 )
        var(ii,[11 12 13]) = { cm ch cl };
    end

    %---------------------------------------------------
    % FillValue

     try
        [f,status] = hdfsd('getfillvalue',vid);
     catch
        status = -1;
     end

     if ~( status == -1 )
          var{ii,8} = f;
     end

    %---------------------------------------------------
    % Range

     try
        [r1,r0,status] = hdfsd('getrange',vid);
     catch
        status = -1;
     end

     if ~( status == -1 )
          var{ii,7} = [ r0  r1 ];
     end

    %---------------------------------------------------

    status = hdfsd('endaccess',vid);

end

status = hdfsd('end',fid);

%*************************************************************
% Cat Dimensions

dim = cat(1,var{:,9});

if size(dim,1) == 1
   return
end

[h,si] = sort(dim(:,1));

dim = dim(si,:);

nd = size(dim,1);
ok = ones(nd,1);

for ii = 2 : size(dim,1)
    ok(ii) = ~isequal( dim(ii-1,[1 2 3 5]) , dim(ii,[1 2 3 5]) );
end

if ~all(ok)
    ok = find(ok);
    dim = dim(ok,:);
end

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function dim = getdim(vid,n)

dim = cell(n,5);  % { Name Length Type Scale {att} }

dim(:,4) = {NaN};

for ii = 1 : n

    did = hdfsd('getdimid',vid,ii-1);

    [dim{ii,1},dim{ii,2},dim{ii,3},natt,status] = hdfsd('diminfo',did);

     try
        [s,status] = hdfsd('getdimscale',did);
     catch
        status = -1;
     end

     if ~( status == -1 );
          dim{ii,4} = s;
     end

     dim{ii,5} = getatt(did,natt);

end

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function att = getatt(vid,n)

att = cell(n,5);  % { Name Count Type Value String }

for ii = 1 : n

    [att{ii,1},att{ii,3},att{ii,2},status] = hdfsd('attrinfo',vid,ii-1);

    [att{ii,4},status] = hdfsd('readattr',vid,ii-1);

     att{ii,5} = att2str(att{ii,4},att{ii,3});

end
 
%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function txt = gettext(file,dim,var,att);

% GETTEXT   Create Text

ndim = size(dim,1);
nvar = size(var,1);

tab  = char(32*ones(1,4));

ntbl  = 8;    % Blanks for Type
nmax  = 240;  % Max Length of AttributeString
nsplt = 60;   % Split AttributeString

%----------------------------------------------------------------
% FileName

txt = { ... 
  ''
 sprintf( 'File: %s' , file )
 sprintf( '------%s' , char('-'+0*file) )
  '' };

%----------------------------------------------------------------
% Global Attributes

if ~isempty(att);

    txt = cat( 1 , txt , ...
              { ' Global Attributes:' 
                ' ------------------' 
                ''   } );

    s0 = size( char(att(:,1)) , 2 );
    nt = sprintf( '\n %s%s' , tab , char(32*ones(1,s0+5+ntbl)) );

    for jj = 1 : size(att,1)

        s2 = size(att{jj,1},2);

        tb = sprintf( '%s%s' , tab , char(32*ones(1,s0-s2)) );

        str = att{jj,5};

        str = splitstr(str,nmax,nsplt);

        str = strrep( str , char(10) , ...
                    [ nt char(32*ones(1,strcmp(str(1),'"'))) ] );

        s3 = size(att{jj,3},2);
        typ = sprintf('%s%s' , att{jj,3} , char(32*ones(1,ntbl-s3)) );

        str = sprintf( '%s%s = #%s  %s' , tb , att{jj,1} , typ , str );

        [m,str] = char2cell(str);  
   
        txt = cat( 1 , txt , str );

    end

    txt = cat( 1 , txt , ...
              { '' 
                '' } );

end


%----------------------------------------------------------------
% Dimensions

%  DIM = { 1.Name 2.Length 3.Type 4.Scale 5.{ATT} }

txt = cat( 1 , txt , ...
{ sprintf(' %2.0f Dimensions:',ndim)
  ' --------------' 
  ''  } );

   s0 = cat(1,dim{:,5});
   if ~isempty(s0)
       s0 = size( char(s0(:,1)) , 2 );
   else
       s0 = 0;
   end
   s0 = max(s0,6);
   nt = sprintf( '\n %s%s%s' , tab , tab , char(32*ones(1,s0+5+ntbl)) );

   s1 = size(char(dim(:,1)),2);

sep = tab(1:end-1);

for ii = 1 : ndim

    s2 = size(dim{ii,1},2);

    tb = char(32*ones(1,s1-s2));

   dim_txt = sprintf('%2.0f)%s%s%s = %.0f',ii,sep,dim{ii,1},tb,dim{ii,2});

   txt = cat( 1 , txt , { dim_txt ; '' } );

   att = dim{ii,5};  % { Name Count Type Value String } 

   if ~strcmp(dim{ii,3},'none')
       att = cat( 1 , { 'Type' 0 '_' [] dim{ii,3} } , att );
   end

   for jj = 1 : size(att,1)

      s2 = size(att{jj,1},2);

      tb = sprintf( '%s%s' , tab , char(32*ones(1,s0-s2)) );

      str = att{jj,5};

      str = splitstr(str,nmax,nsplt);

      str = strrep( str , char(10) , ...
                    [ nt char(32*ones(1,strcmp(str(1),'"'))) ] );

      s3 = size(att{jj,3},2);
      typ = sprintf('%s%s' , att{jj,3} , char(32*ones(1,ntbl-s3)) );

      str = sprintf( '%s%s%s = #%s  %s' , tab , tb , att{ii}{jj,1} , typ , str );

      [m,str] = char2cell(str);

      txt = cat( 1 , txt , str );

  end

  if ~isempty(att)
      txt = cat( 1 , txt , {''} );
  end

end

%----------------------------------------------------------------
% Variables

% var = { 1.Name 2.Rank 3.Type 4.Ndim 5.[Size] 6.Natt 7.Fill 8.[Range] 9.{dim} 10.{att} ...
%         11.Compressed 12. Chunked  13.ChunkLenght }

txt = cat( 1 , txt , ...
{ ''
  ''
  sprintf(' %2.0f Variables:',nvar)
  ' -------------'
  '' } );


   s0 = cat(1,var{:,10});
   if ~isempty(s0)
       s0 = size( char(s0(:,1)) , 2 );
   else
       s0 = 0;
   end
   s0 = max(s0,6);

   nt = sprintf( '\n %s%s%s' , tab , tab , char(32*ones(1,s0+5+ntbl)) );

sep = tab(1:end-1);

inf = { 'compressed' 'chunked' };

for ii = 1 : nvar

   var_txt = sprintf('%2.0f)%s%s',ii,sep,var{ii,1});

   if ~isempty(var{ii,9})
      % Dimensions
      var_txt = sprintf( '%s( %s )' , var_txt , strhcat(var{ii,9}(:,1),' , ') );
    end
 
   txt = cat( 1 , txt , { var_txt ; '' } );

   att = var{ii,10};  % { Name Count Type Value String } 


   if var{ii,11} | var{ii,12}
      if var{ii,12}
         att = cat( 1 , { 'ChunkLength' 0 '_' [] sprintf('%.0f ',var{ii,13}) } , att );
      end
      att = cat( 1 , { 'Info' 0 '_' [] ...
                       strhcat(inf(find(cat(2,var{ii,[11 12]}))),', ') } , att );  
   end

   att = cat( 1 , { 'Size' 0 '_' [] sprintf('%.0f ',var{ii,5}) } , att );

   att = cat( 1 , { 'Rank' 0 '_' [] int2str(var{ii,2}) } , att );

   att = cat( 1 , { 'Type' 0 '_' [] var{ii,3} } , att );

   for jj = 1 : size(att,1)

      s2 = size(att{jj,1},2);

      tb = sprintf( '%s%s' , tab , char(32*ones(1,s0-s2)) );

      str = att{jj,5};

      str = splitstr(str,nmax,nsplt);

      str = strrep( str , char(10) , ...
                    [ nt char(32*ones(1,strcmp(str(1),'"'))) ] );

      s3 = size(att{jj,3},2);
      typ = sprintf('%s%s' , att{jj,3} , char(32*ones(1,ntbl-s3)) );

      str = sprintf( '%s%s%s = #%s  %s' , tab , tb , att{jj,1} , typ , str );

      [m,str] = char2cell(str);

      txt = cat( 1 , txt , str );

  end

  if ~isempty(att)
     txt = cat( 1 , txt , {''} );
  end

  txt = cat( 1 , txt , {''} );

end

 txt = strrep(txt,char(0),'\0');


%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = splitstr(str,nmax,nsplt);

if isempty(str)
   return
end

n = size(str,2);

is_char = all( str([1 n]) == '"' );

if all( n <= [ nsplt+2*is_char nmax+2*is_char ] )
   return
end

if n > nmax+2*is_char

   suf = ' ...';
   ns  = size(suf,2);

   ind = cat( 2 , ( 1 : nmax+is_char ) ,  n*ones(1,is_char) );

   str = str(ind);

   ind = ( 1 : ns ) + nmax-ns+is_char;

   str(ind) = suf;

end

n = size(str,2);
i0 = cat( 2 , 1 , find(str==10)+1 );
lg = diff(cat(2,i0,n+1));

if ~any( lg > nsplt )
    return
end

i1 = zeros(1,0);
for cc = ',;.'
    i1 = cat(2,i1,findstr(str,sprintf('%s ',cc)));
end

if isempty(i1)
   return
end

i1 = sort(i1);

str = insert(str,i1+1,char(10));


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
     fprintf( fid , strhcat(txt,char(10)) );
     fclose(fid);
  else
     msg = sprintf('Error open OutFile "%s" for writing.',file);
  end

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  fig = list(file,dim,var,txt,fig)

% LIST  Display Info in ListBoxFigure

Nin = nargin;

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end

tag = upper(fcn);

new = ( Nin < 5 );
if ~new
    new = ~chkhndl(fig,'figure',tag);
end

if new
   fig = newfig(tag,fcn);
end

set( fig , 'name' , file );

hl = findobj(fig,'type','uicontrol','tag','LIST');

set( hl , 'listboxtop' , 1 , ...
               'value' , 1 , ...
              'string' , txt     );

ha = findobj(fig,'type','uimenu','tag','ASSIGN');

delete(get(ha,'children'));

set( ha , 'userdata' , file );

nn = 0;

for ii = 1 : size(var,1)

    h0 = uimenu( 'parent'   , ha , ...
                 'label'    , sprintf(' %s',var{ii,1}) , ...
                 'callback' , '' , ...
                 'userdata' , var{ii,1}  );

    if var{ii,4} == 1
       si = [ var{ii,5}  1  ];
    else
       si = var{ii,5};
    end

    si = si(var{ii,4}:-1:1);

    lab = sprintf( '%.0f x ' , si  );
    lab = sprintf('%s   %s',lab(1:end-3),var{ii,3});

    h1 = uimenu( 'parent'   , h0  , ...
                 'label'    , lab , ...
                 'userdata' , si  , ...
                 'callback' , sprintf('%s(''#assign'',gcbo);',fcn) );

    nn = nn + prod(si);

end

    h0 = uimenu( 'parent'   , ha         , ...
                 'label'    , ' All ... ' , ...
                'separator' , 'on' , ...
                 'callback' , ''   , ...
                 'userdata' , var(:,1)   );

    lab = sprintf( 'max %.1f MBytes' , nn*8 / 2^20  );

    h1 = uimenu( 'parent'   , h0  , ...
                 'label'    , lab , ...
                 'userdata' , si  , ...
                 'callback' , sprintf('%s(''#assign'',gcbo);',fcn) );

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
               'menubar'    , 'none' , ...
               'toolbar'    , 'none' , ...
               'name'       , tag    , ...
               'tag'        , tag    , ...
               'createfcn'  , ''     , ...
          'handlevisibility','callback' );


  hl = uicontrol( 'parent'    , fig        , ...
                  'style'     , 'listbox'  , ...
            'backgroundcolor' , [1 1 1]    , ...
            'foregroundcolor' , [0 0 0]    , ...
                  'units'     ,'normalized', ...
                  'position'  , [0 0 1 1] , ...
                  'min'       , 0         , ...
                  'max'       , 2         , ...
                  'fontunits' , 'points'  , ...
                  'fontsize'  , fontsize  , ...
                  'fontname'  , fontname  , ...
                  'string'    , ''        , ...
                  'tag'       , 'LIST'    , ...
        'horizontalalignment' , 'left'           );


    ParentCB = 'get(gcbo,''parent'')';

       NewCB = sprintf('%s(''#new'',%s);',fcn,ParentCB);
   RefreshCB = sprintf('refresh(%s);',ParentCB);
     CloseCB = sprintf('delete(%s);',ParentCB);
   
   hn = uimenu( 'parent'      , fig   , ...
                'accelerator' , 'N'   , ...
                'label'       , 'New' , ... 
                'tag'         , 'NEW' , ...
                'callback'    , NewCB       );

   ha = uimenu( 'parent'      , fig      , ...
                'accelerator' , 'A'      , ...
                'label'       , 'Assign' , ...
                'tag'         , 'ASSIGN' , ...
                'callback'    ,  ''            );
 
   hr = uimenu( 'parent'      , fig      , ...
                'accelerator' , 'R'      , ...
                'label'       , 'Refresh' , ...
                'tag'         , 'REFRESH' , ...
                'callback'    , RefreshCB     );
 
   hc = uimenu( 'parent'      , fig     , ...
                'accelerator' , 'C'     , ...
                'label'       , 'Close' , ...
                'tag'         , 'CLOSE' , ...
                'callback'    , CloseCB       );

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function new(fig)

% CallBack for New-Menu

   [f,p] = uigetfile('*.nc','Select a HDF-File ...');
   if isequal(f,0)
      return
   end

   file = fullfile(p,f);

   [msg,dim,var,att,txt] = look_hdf(file);

   if ~isempty(msg)
       msg = cat( 2 , 'Error using LOOK_HDF' , char(10) , msg );
       warndlg(msg,'Error','warn');
       return
   end

   list(file,dim,var,txt,fig);

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function assign(h1)

% ASSIGN  CallBack for Assign-Menu

   h0 = get( h1 , 'parent' );
   ha = get( h0 , 'parent' );

   file = get( ha , 'userdata' );
   var  = get( h0 , 'userdata' );

   [ok,var] = chkcstr(var,0);

   if ~ok
       msg = 'Invalid UserData of UIMenu.';
   else
       try
          [msg,d,v] = read_hdf(file,'var',var);
          if isempty(msg)
             for ii = 1 : size(v,1)
                 assignin('base',strrep(v{ii,1},' ','_'),v{ii,14});
             end
          end
       catch
          msg = lasterr;
       end 
   end
  
   if ~isempty(msg)
       msg = cat( 2 , 'Error Assign Variable' , char(10) , msg );
       warndlg(msg,'Error','warn');
       return
   end

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = type_cmp(typ,varargin);

ok = 0;
if nargin < 2
   return
end

[vok,cmp] = chkcstr(varargin(:));

if ~vok
    return
end

for ii = 1 : size(cmp,1)
    if size(typ,2) >= size(cmp{ii},2)
       ok = ~isempty(findstr(typ,cmp{ii}));
       if ok
          break
       end
    end
end

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = att2str(val,typ)

%  STRING = ATT2STR( Value , Type )
%
%  Converts AttributeValues to String
%

str = '';

if strcmp(typ,'none');
   typ = class(val);
end

%-----------------------------------------------------------
if type_cmp(typ,'char')

  if isequal(val,0)
    val = '\0';
  end

  val = checkchar(val);

  str = sprintf('"%s"',val);

%-----------------------------------------------------------
elseif type_cmp(typ,'int','long','short','byt')

  str = int2str(double(val));

%-----------------------------------------------------------
else

  str = num2str(double(val));

end


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = checkchar(str,mode)

% CHECKCHAR  Checks for valid Characters

if nargin < 2
   mode = 1;
end

  str = double(str);

  bad = find( ~( ( str ==  9 ) |  ...
                 ( str == 10 ) |  ...
                 ( str == 13 ) |  ...
                 (  28 <= str  &   str <= 126 ) | ...
                 ( 160 <= str  &   str <= 255 )        ) );

str(bad) = 32;


if mode
  str = rmblank(str,2);
else
  str = char(str);
end

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [x,ind,ix] = insert(x,ind,y,dim)

% INSERT Inserts Values into Array
%
% [ Z , IY , IX ] = INSERT( X , Index , Y , DIM )
%
% inserts Y into X at the Index in the Dimension DIM.
%
% defaults: DIM = 1, for Vectors X the first non-singleton Dimension
%
% Index >  0  ==>  Y will inserted behind X(Index)
% Index <  0  ==>  Y will inserted before X(Index)
% Index == 0  ==>  Y will inserted before X(1)
%
% Make sure, that Index match the Size of X in DIM.
% 
% Y can be Empty, a Scalar, Vector or Matrice.
%
%  Make sure, that the Size of Y match the Size of X and 
%   the Length of Index in DIM.
%
% If Y is empty or not given, X(Index) will inserted (duplicated).
%
%--------------------------------------------------------------
% The Output-IndexVectors IY and IX refers to the inserted Y 
%   and the original X.
%
% IY = [ N by 2 ];  with ( N <= length(Index) )  
%
%   The 1. Column of IY is the Index in Z of the Inserted  Y
%   The 2. Column of IY is the Index in Y of the Inserted  Y
% 
% IX = [ size(X,DIM) by 1 ] is the Index in Z of the Original X
%
%--------------------------------------------------------------
% for Vectors X and Vectors or Scalars Y
%
%  Z( IY(:,1) ) == Y( IY(:,2) ) 
%
%  Z(IX)        == X
% 
%--------------------------------------------------------------
% for Matrice X and Y
%
%  Z( : , ... , IY(:,1) , ... , : ) == Y( : , ... , IY(:,2) , ... , : ) 
%
%  Z( : , ... , IX , ... , : )      == X
%


Nin = nargin;

Nout = nargout;

ix = [];

if Nin < 2
   error('Not enough Input Arguments.');
end

if Nin < 3
   y = [];
end

if Nin < 4
   dim = [];
end

%--------------------------------------------------------
% Check Class of X and Y

if ~strcmp(class(x),class(y))
    error('X and Y must be of same Class.');
end

%--------------------------------------------------------
% Check Size of X and Index

sx = size(x);

if isempty(dim)
   [ms,dim] = max(sx);
   if ~( ms == prod( sx + ( sx == 0 ) ) )  % Not a Vector
      dim = 1;
   end
else
   sx = cat( 2 , sx , ones(1,dim-size(sx,2)) );
end

if all( sx == 0 )
   sx      = ones( 1 , dim+(dim==1) );
   sx(dim) = 0;
   x       = zeros(sx);
end

ind = floor(ind(:));
 ni =  size(ind,1);

ind = abs(ind) - 1 * ( ind < 0 );  % Negative Inserts before !!!

if ~all( ( 0 <= ind ) & ( ind <= sx(dim) ) );
    error(sprintf('Index exceeds Matrix Dimension %.0f.',dim));
end

%--------------------------------------------------------
% Check Size of Y

sy = size(y);

not_y =  all( sy == 0 );   % ZERO-Size
one_y =  all( sy == 1 );   % Scalar

vec_y = ( max(sy) == prod( sy + ( sy == 0 ) ) ); % Scalar | Vector
vec_y = ( vec_y & ~one_y );                      % Vector

mat_y = ~( not_y | one_y );                      % Vector | Matrice

%--------------------------------------------------------
% Vector ==> reshape into Dimension DIM

if vec_y 

   ss = cat( 2 , sy , ones(1,size(sx,2)-size(sy,2)) );
   ss(dim) = sx(dim);

   if ~isequal(ss,sx)

       [ny,dy] = max(sy);

       if ~( dy == dim )

           p      = ( 1 : size(sy,2) );
           p(dim) = dy;
           p(dy ) = dim;

                y = permute( y , p );

               sy = size(y);

       end

    end

end

%--------------------------------------------------------
% Check Size of Y with Size of X and Length of Index at DIM

if mat_y

    sy = cat( 2 , sy , ones(1,dim-size(sy,2)) );

    s      = sx;
    s(dim) = ni;

    if ~isequal( s , sy )
       error(sprintf('Size of Y must match Size of X and Length of Index in Dimension %.0f.',dim));
    end

end


%--------------------------------------------------------
% Permute and Reshape X and Y to 1. Dimension at DIM

ns = size(sx,2);

perm = cat( 2 , dim , ( 1 : dim-1 ) , ( dim+1 : ns ) );

x = reshape( permute(x,perm) , sx(perm(1)) , prod(sx(perm(2:ns))) );

if mat_y
   y = reshape( permute(y,perm) , sy(perm(1)) , prod(sy(perm(2:ns))) );
end

%--------------------------------------------------------
% Remove Multiple Indize

[ind,si] = sort(ind,1);

    bad  = find( diff(ind,[],1) == 0 );

ind(bad) = [];

ni = size(ind,1);

if mat_y
     y        = y(si,:);
     y(bad,:) = [];
    si(bad)   = [];
else
    si        = ones(ni,1);
    if not_y
       si = NaN * si;
    end
end

%--------------------------------------------------------
% Build IndexVector to Insert

ii = ones( sx(dim)+ni , 1 );

ind = ind + ( 1 : ni )';

ii(ind) = 0;

if Nout > 2
   ix = find(ii);
end

if ~isempty(ii)

  ii    = cumsum(ii,1);

  ii(1) = ii(1) + ( ii(1) == 0 );   % Check for ZERO-Index

end

% Insert X (duplicate)
if ~isempty(x)
   x = x(ii,:);
end

% Insert Y
if ~not_y
   x(ind,:) = y;
end

if Nout > 1
   ind = cat( 2 , ind , si );
end

%--------------------------------------------------------
% Reshape and Permute back

sx(dim) = sx(dim) + ni;

x = reshape( x , sx(perm) );

perm = cat( 2 , ( 1 : dim-1 ) + 1 , 1 , ( dim+1 : ns ) );

x = permute( x , perm );

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
 n = [];
end
if Nin < 2
 del = char((32*ones(1,3)));
end


if isempty(str)
 str = '';
 return
end


if ischar(str)
  str = cellstr(str);
end

str = str(:);

if isempty(n)
   n = size(str,1) + 1;
end

str( : , 2 ) = { del };

str(n:n:size(str,1),2) = { nl };

str(    size(str,1),2) = { '' };


str = permute( str , [ 2  1 ] );

str = cat(2,str{:});


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [Msg,bb,ref] = char2cell(bb,mnl,ref01);

% CHAR2CELL Converts CharArray to CellStringArray
%
% [Msg,CellString] = CHAR2CELL( String )
%
% Converts the CharacterArray String into a CellStringArry.
%
%---------------------------------------------------------
% CHAR2CELL( String , MarkerNewLine )  
%
% Replace the String's defined by MarkerNewLine with NewLine,
%   default: EMPTY
%
%---------------------------------------------------------
% [Msg,CellString,Reference] = ...
%    CHAR2CELL( String , MarkerNewLine , ReferenceCell)
%  
% Returns References in String, using GET_REF, 
%   ReferenceCell = { S1 S2 C1 C2 }, type: >> help get_ref
%


 Msg = '';
 
 Msg0 = 'CHAR2CELL: ';

 ref = cell(0,2);

 Nin  = nargin;
 Nout = nargout;
 
 nl = char(10);

  if isempty(bb)
     bb = cell(0,1);
     return
  end

  % Marker for NewLine
  if Nin < 2
    mnl = '';
  end

%*****************************************************
% Check Inputs

  if iscellstr(bb)
     bb = cellstr(char(bb));
     bb = strhcat(bb,'',1);
  end

  if isempty(bb)
     bb = cell(0,1);
     return
  end

  ok = ( isnumeric(bb) | ischar(bb) );
  if ok

    if ischar(bb)
       bb = double(bb);
    end

    ok = all( ( mod(bb,1) == 0 )  & ( bb >= 0 ) & isfinite(bb)  );

  end

  if ~ok

      Msg = [ Msg0 ...
              'Input String must be a CharacterArray or ASCII-Codes.'];
 
      bb = cell(0,1);

      return

  end

  %---------------------------------------------------
  % Check MarkerNewLine

  if ~( isempty(mnl) |  ...
        ( ischar(mnl) &  ( prod(size(mnl)) == size(mnl,2) ) ) )

    Msg = [ Msg0  'Input MarkerNewLine must be a String.'];

    return

  end

  %---------------------------------------------------
  % Check Reference

  if ( Nin == 3 )  &  ( Nout == 3 )

     [Msg,ref] = get_ref('',ref01{:});

    if ~isempty(Msg)
       Msg = [ Msg0  'Invalid Input for Reference.' nl Msg ];
       return
    end

  end

%*****************************************************

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
    Msg = [Msg0 'Invalid Characters in String.' ];
    return
  end


%*****************************************************
 

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

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [ok,msg] = chkhndl(h,typ,tag);

% CHKHNDL(H,Type,Tag)  Checks, if H is a Handle of specified Type and Tag
%


Nin = nargin;

ok  = 0;
msg = [];


if Nin == 0
   return
end

ok = ( isnumeric(h) & ( prod(size(h)) == 1 ) );

if ~ok
   return
end

ok = ishandle(h);
if ~ok | ( Nin < 2 )
   return
end

%-------------------------------------------------------------------------
% Check with Type

is_typ = chkstr(typ,0);

if is_typ

  ok = isempty(typ);
  if ~ok
      ok = strcmp( get(h,'type') , typ  );
  end

  if ~ok | ( Nin < 3 )
     return
  end

elseif ( Nin == 3 )

  msg = 'Input Type must be a String.';
  ok  = 0;
  return

else
 
  tag = typ;
  
end


%-------------------------------------------------------------------------
% Check Tag

[ok,tag] = chkcstr(tag,0);

if ~ok
    msg = 'Input Tag must be a CharArray or CellStringArray.';
    return
end   


%-------------------------------------------------------------------------
% Check with Tag

 t = get(h,'tag');

nt = size(t,2);

ok = strcmp(t,tag);


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
   str = cellstr(str);
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

   
