function [msg,fl] = mkcont(pfad,opt,name,mp,cp)

% MKCONT  Creates "Contents.m"- File of Directory
%
% MKCONT( Directory , [Option] , [Prefix] )
%
% Option = '-r' | '-rl'  
%
%          '-r'   works recursivly
%          '-rl'  works recursivly, follow Links
%
% First Line of "Contents.m":
%
% %# Contents of Directory "[Prefix]/DirName/SubDir/..."
%
% If a File "Comments.txt" exist in the Directory,
%  it's text is written at the begin of the "Contents.m"
%
%  see also:  DIRCONT
%

Nin  = nargin;
Nout = nargout;


pre = '# Contents of "';
suf = '"';
prf = '%   ';
sep = '  - ';

ntf = 'No help comments found';

file = 'Contents.m';
alt  = 'Content.m';
comm = 'Comments.txt';

fs = filesep;
ps = pathsep;

if Nin
   mp = matlabpath;
   cp = seperate(mp,ps);
end

fl = '';  % First Line of Comments.m

%--------------------------------------------

if Nin < 1
   pfad = '';
end

if Nin < 2
   opt = '';
end

lnk = 0;
rec = ~isempty(opt);
if rec
   rec = ( isequal(opt,'-r') | isequal(opt,'-rl') );
   if rec
      lnk = isequal(opt(end),'l');
   end
end

if Nin < 3
   name = '';
end

%--------------------------------------------
% Try to get absolute Pfad, Ziel

fs = filesep;

if isempty(pfad)
   pfad = cd;
end

if ~( exist(pfad,'dir') == 7 )
    error('Directory doesn''t exist.');
end

% Absolute Path !!!

p0 = cd;

try
   cd(pfad);
   pfad = cd;
end

cd(p0);
   
%-------------------------------------------------------------

pfad = cat( 2 , pfad , fs(1:(end*(~strcmp(pfad(end),fs)))) );

pfd  = pfad(1:(end-1));

pname = pfd;

ii = find( double(pname) == double(fs) );
if ~isempty(ii)
   pname = pname(max(ii+1):end);
end

name = cat( 2 , name , fs(1:(end*(~isempty(name)))) , pname );

%-------------------------------------------------------------

d = dir(pfad);

if isempty(d)
   fprintf(1,'%s   can''t read Directory\n',pfad)
   return
end

d = d(find(cat(1,d.isdir)));

d = { d.name };

d = d(find(~( strcmp(d,'.') | strcmp(d,'..') )));

%-------------------------------------------------------------
% Recurse

if rec & ~isempty(d)

   d = sort(d(:));

   n = size(d,1);

   c  =  cell(n,1);
   ok = zeros(n,1);
   nl = char(10);

   for ii = 1 : n
 
       p = cat(2,pfad,d{ii});

       ok(ii) = ( lnk | ~isunix );
       if ~ok(ii)
           ok(ii) = ~islink(p);
       end

       if ok(ii)
          [msg,c{ii}] = mkcont( p , opt , name , mp , cp );
          ok(ii) = ~( isempty(c{ii}) | strcmp('private',d{ii}) );
          if ok(ii)
              d{ii} = cat( 2 , d{ii} , fs );
          end             
       end

   end

   if ~any(ok)
       c = '';
   else

       if ~all(ok)
            n = sum(ok);
           ok = find(ok);
            c = c(ok);
            d = d(ok);
       end

       d = char(d);
       d(find(d==' ')) = '@';
       d = cellstr(d);
       d = strrep(d,'@',' ');

       o = ones(n,1);

       c = cat( 2 ,  prf(o,:) , char(d) , sep(o,:) , char(c) );

       c = cat(1,cellstr(c),{prf});

       c = sprintf('%s\n',c{:});

   end

else

   c = '';

end

%-------------------------------------------------------------
% Look for M-Files

fprintf(1,'%s',pfad);

d = dir(cat(2,pfad,'*.m'));

if ~isempty(d)
    d = d(find(~cat(1,d.isdir)));
end

if ~isempty(d)
    d = {d.name};
    [p,n,e] = fileparts(alt);
    jj = strmatch(n,d);
    if ~isempty(jj)
        d(jj) = [];
    end
else
    d = cell(1,0);    
end

cm = dir(cat(2,pfad,comm));

if prod(size(cm)) == 1
   if ~cm.isdir
       [p,n,e] = fileparts(cm.name);
        d = cat( 2 , {[n e]} , d );
   end
end

if isempty(d) & isempty(c)

   msg = 'no M-Files or Comments found';

else

   % Check if PFAD in MatlabPath !!!
   %  to work correct using HELP

   ok = ( any(strcmp(cp,pfd)) | any( pfd == '@' ) );
   ok = ( ok | ~isempty(findstr([fs 'private'],pfd)) );
   if ~ok 
       matlabpath(cat(2,mp,ps,pfd));
   end

   [msg,fl] = wrt_cont(pfad,file,alt,comm,c,d,name,pre,suf,prf,sep,ntf);

   if ~ok
       matlabpath(mp);
   end

end

fprintf(1,'   %s\n',msg);

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,fl] = wrt_cont(pfad,file,alt,comm,c,d,name,pre,suf,prf,sep,ntf)

msg = '';

fl  = '';  % First Line of Comments.txt

str = '';  % String of FIleHelp

%********************************************************
if ~isempty(d)
%********************************************************

ic = strcmp(d{1},comm);

if ic
   cm = d{1};
    d = d(2:end);
end

n  = prod(size(d));

d = sort(d);

str = cell(n,1);

fc = 'function';

nc = size( fc , 2 );

nl = char(10);

ok = zeros(n,1);

for ii = 1 : n

    f = lower(d{ii});

    ok(ii) = ( ~isempty(f) & ~strcmp(f,lower(file)) );

    if ok(ii)
        
        f = f( 1 : end-2 );

        nf = size(f,2);

        try
          txt = help( cat(2,pfad,f) );  % Works only if PFAD in MatlabPath !!!
        catch
          txt = '';
        end

        ok(ii) =  ~isempty(txt);

        if ok(ii)

            %-----------------------------------------------
            % Extract first HelpLine

            txt = cat( 2 , double( rmblank(txt,2) ) , 10 );

            txt( find( txt == 13 ) ) = [];

            txt = txt( 1 : min( find( txt == 10 ) ) - 1 );

            txt = char(txt);

            if ~isempty(txt)
                txt = txt( 1 : (end-strcmp(txt(end),'.')) );
            end

            %-----------------------------------------------
            % Remove FunctionName in Help
         
            nt = size(txt,2);

            if ~isempty(txt) & ( nt >= nf )
                if strcmp( lower(txt(1:nf)) , f )  &  ...
                  ~strcmp( txt(min(nt,nf+1)) , '(' )
                   txt = rmblank( txt( nf+1 : end ) );
                end
            end

            %-----------------------------------------------
            % Start with Letter or Number

            if ~isempty(txt)
                jj = ~( isletter(txt) | ...
                        ( ( '0' <= txt ) & ( txt <= '9' ) ) );
                if jj(1)
                   if all(jj)
                      txt = '';
                   else
                      jj = sum(cumprod(double(jj),2),2) + 1;
                      txt = txt(jj:end);
                   end
                end
             end
        end
        % ~isempty(txt)

        ok(ii) =  ~isempty(txt);

        if ok(ii)

            %-----------------------------------------------
            % Check for FUNCTION-Statement in Help

            gd = ( size(txt,2) < nc );
            if ~gd
                gd = ~strcmp( txt(1:nc) , fc );
            end

            if gd
               txt(1) = upper(txt(1));
            end

        else

            txt = ntf;

        end

        %----------------------------------------------------- 

        str{ii} =  txt;

          d{ii} = d{ii}( 1 : end-2 );

    end
    % ok

end
% ii

if ~any(ok)

    str = '';

else

    if ~all(ok)
         n = sum(ok);
        ok = find(ok);
       str = str(ok);
         d =   d(ok);
    end

    d = char(d);
    d(find(d==' ')) = '@';
    d = cellstr(d);
    d = strrep(d,'@',' ');

    o = ones(n,1);

   str = cat( 2 ,  prf(o,:) , char(d) , sep(o,:) , char(str) );

   str = cat(1,cellstr(str),{prf});

   str = sprintf('%s\n',str{:});

   msg = sprintf('%.0f files',n);

end

if  ~isempty(c)
    str = cat( 2 , c , str );
end


%********************************************************
% Check for Comments.txt

if ic
    [m,cm] = loadfile(cat(2,pfad,cm),'char',20*2^10); % Max. 20 kB
    if isempty(m)
       cm(find(cm==char(13))) = [];
       cm = cat(2,nl,cm,nl);
       cm = strrep(cm,[nl '%'],nl);
       cm = rmblank(cm,2);
    else
       cm = '';
    end 
    ok = ~isempty(cm);
    if ok
       jj = ( cm == nl );
       if ~any(jj)
          fl = cm;
          cm = '';
       else
          jj = min(find(jj));
          fl = cm( 1 : jj-1 );
          cm = rmblank(cm((jj+1):end) ,2 );
       end
       ok = ~isempty(cm);
    end
    if ok
       cm = strrep(cm,nl,[nl '% ']);
       cm = sprintf('%% %s\n%s\n',cm,prf);
       if ~isempty(str)
           sp = 0;
           jj = ( str == nl );
           if sum(jj) > 2
              sp = max(diff(find(jj)));
           else
              jj = ( cm == nl );
              if sum(jj) > 2
                 sp = max(diff(find(jj)));
              end
           end
           if sp == 0
              sp = '%';
           else
              sp = ones(1,sp);
              sp = char('-'*sp);
              sp(1) = '%';
           end
           cm = sprintf('%s%s\n%s\n',cm,sp,prf);
       end
        str = cat( 2 , cm , str );
    end
end

if isempty(str) & isempty(fl)
   msg = 'no Help found';
   return
end

%********************************************************
elseif ~isempty(c)
%********************************************************

   str = c;

%********************************************************
end
%********************************************************

%********************************************************
% Check for Contents.m

chk = chkfile(pfad,file,pre);
if ~chk
    m  = sprintf('%s allready exist',file);
    cc = -1;
    while ~chk & ( cc < 9 )          
       cc = cc + 1;      
       ff  = sprintf('%s%.0f%s',file(1:end-2),cc,file(end-1:end));
       chk = chkfile(pfad,ff,pre);
       if chk
          file = ff;
       end
    end
    if ~chk
        msg = m;
        return
    else
        fprintf(1,'   %s, write %s',m,file);
    end
end


%********************************************************
% Write File

fid = fopen( cat(2,pfad,file) , 'wt' );

if isempty(fl)
   pre = sprintf('%%%s',pre);
else
   pre = sprintf('%% %s %s',fl,pre);
end

if fid == -1
   msg = sprintf('can''t open%s',file);
else
   fprintf(fid,'%s%s%s\n%s\n%s',pre,name,suf,prf,str);
   fclose(fid);
end


%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = chkfile(pfad,file,pre)

ok = 0;

file = cat(2,pfad,file);

if ( exist(file,'dir')  == 7 )
   return
end

ok = ~( exist(file,'file') == 2 );
if ok
   return
end


fid = fopen(file,'r');

if ( fid == -1 )
   return
end


bb =  fgetl(fid);
ok = ~isempty(findstr(pre,bb));

fclose(fid);


%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,bb] = loadfile(file,varargin);

% LOADFILE  Load binary Data from File, using FREAD
%
% [ Msg , V ] = LOADFILE( FileName , MaxSize , Precision , Extension )
%
%   MaxSize     MaximumFileSize [Bytes],   default: 2097152 == 2MBytes
%   Precision   Precision for using FREAD, default: 'char'
%   Extension   optional Extension to expand FileName,id neccessary
%               '.ext'
%
%  For more Informations to Precision  type: >> help fread
%
%  In case of Precision 'char', a CharacterString will returned,
%   if all Bytes are valid Characters for conversion, if some
%   Characters are invalid, Msg is not empty.
%

msg = '';
bb  = [];

msg0 = 'LOADFILE: ';

%------------------------------------------------
% Defaults

m   = 2*2^20; % 2 MBytes;  % MaximumFileSize [Bytes]
p   = 'char';
ext = '';

Nin = nargin;

%************************************************
% Check File

if Nin == 0
   msg = [msg0 'Input File is missing.'];
   return
end

if isempty(file)
   return
end

if ~chkstr(file,1)
   msg = [msg0 'Input File must be a String.'];
   return
end

%------------------------------------------------
% Get MaxSize and Precision

for ii = 1 : Nin-1

    v = varargin{ii};

    if chkstr(v,1) & ~isempty(v)

       if strcmp(v(1),'.')
          ext = v;
       else
          p = v;
       end
 
    elseif ( isnumeric(v)  &  ( prod(size(v)) == 1 ) )

       if ( isfinite(v)  &  ( v >= 0 ) )
          m = v;
       end

    end

end

%************************************************
% Open and Read File

  fid = fopen(file,'r');

  if fid == -1 
     if ~isempty(ext)
        ie = 1 + strcmp( file(end) , ext(1) );
        ff = cat( 2 , file , ext(ie:end) );
        fid = fopen(ff,'r');
     end
     if fid == -1   
        msg = [ msg0  'Can''t open File.' ];
        return
     end
     file = ff;
  end

 %----------------------------------------------
 % Check Size of File

  d = dir(file);

  if isempty(d)
   d = dir( which(file) );  % which(file) gives the full Name
                            %  [ PathName FileName ]
  end

  if d.bytes > m
    msg = [ msg0 'File too large, Limit = '  ...
            sprintf('%.0f Bytes',m) '.' ];
    fclose(fid);
    return
  end


 %----------------------------------------------

  try
     bb = fread(fid,p);
  catch
     msg = [ msg0 'Error call FREAD.' char(10) lasterr ];
  end

  fclose(fid);

 %----------------------------------------------
 % Precision: 'char'  ==>  Transform to String

  if isempty(msg) & isequal( p , 'char' )
    
    % Check Characters

    if any( ( 126 < bb ) & ( bb < 160 ) );

       % Old DOS !!!
       old = [ 132   148   129   142   153   154  225 ];
       new = [ 'ä'   'ö'   'ü'   'Ä'   'Ö'   'Ü'  'ß' ];

       n  = size(old,2);
       ok = zeros(1,n);
 
       for ii = 1 : n
           if ( ( ii < n ) | any(ok(1:(n-1*(ii>1)))) )
              jj = find(bb==old(ii));
              ok(ii) = ~isempty(jj);
              if ok(ii)
                 bb(jj) = double(new(ii));
              end
           end
       end 
          
    end

    bb(find(bb==12)) = [];   % FormFeed

    if all( ( bb ==  9 ) |  ...
            ( bb == 10 ) |  ...
            ( bb == 13 ) |  ...
            (  28 <= bb  &   bb <= 126 ) | ...
            ( 160 <= bb  &   bb <= 255 )        );

       bb = char(bb(:)');

    else

       msg = [ msg0 'Invalid Characters in File.' ];

    end


  end

%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

