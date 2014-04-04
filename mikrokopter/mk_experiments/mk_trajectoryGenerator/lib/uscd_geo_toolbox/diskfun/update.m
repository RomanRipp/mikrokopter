function c = update(pfad,varargin)

% UPDATE   Update Files in Directory  by Suffix
%
% C = UPDATE( Path , STR1 , STR2 , ... , Options )
%
% Update Files in Path with the Syntax: *SUF#
%   with Files returns by WHICH(FileName)
% 
%  STR# = .ext  Search for "*.ext"
%  STR# = name  Search for "name.m" or "name.m" , WildCard '*' is allowed.
%  STR# = name.ext  Search for "name.ext", WildCard '*' is allowed.
%
% C = UPDATE( Path , NAME1 , NAME2 , ... , '-e' )
%
% Updates only Files with the exact Name: NAME#
%
%----------------------------------------------------------------------------
%
% !!! Please Handle with Care !!!
%
% !!! It's strictly recommended to run a TestMode with the Option '-t' before
%     to verify the Actions and Updates
%     and to  create a Backup with the Option '-b' !!!!!!
% 
% The Inputs for FileSuffixes SUF# must be Strings.
%
% Default: Path = pwd
%          SUF  = '.m';
%
% C returns a CellStringArray with 3 Columns:
%
%  { File UpdateFile BackupFile }
%
% The Actions to create BackupDirectories and BackupFiles and copy Files
%  are displayed in the MatlabCommandWindow.
%
% Options:
%
%  '-t': !!! Test Mode, display Files to backup/update but don't copy
%
%  '-b': !!! Create Backups of updated Files
%                             The User is requested to Enter a DirectoryName for Backup
%
%  '-bNAME': Create Backups of updated Files in Directory NAME
%
%  '-n': updates only with newer Files
%
%  '-e': updates only Files matching exact Name
%
%  '-r': recurse into SubDirectories
%
%  '-R': recurse into SubDirectories which are Links
%
%  '-c': recurse into Class-Directories, starts with '@'
%
%  '-p': recurse into Private-Directories 'private/'
%
% The Options '-c' and '-p' should be used together with '-r' or '-R'
%
% The M-File "Contents.m" will never updated.
%
%
% see also: WHICH, WHICHFILE, DIR, MKDIR, COPYFILE, RECPATH, BACKUP
%

Nout = nargout;

c = cell(0,3);    % { Org New Bcp }

fs = filesep;

%***************************************************************
% Check Inputs

%-----------------------------------------------------
% Pfad

if nargin == 0
   pfad = [];
end

if  isempty(pfad)
    pfad = cd;
elseif ~chkstr(pfad);
    error('Path must be a String.');
elseif ~( exist(pfad,'dir') == 7 )
    error(sprintf('Path doesn''t exist: %s',pfad));
end

%---------------------------------------------
% Get full PathName

p0 = pwd;

try
   cd(pfad);
   pfad = cd;
   cd(p0);
end
 
%-----------------------------------------------------
% Options

rec = NaN;

opt = { '-t'  '-n'  '-p'  '-c' };

def = zeros(size(opt));

bcp = '';

pref = '*';

if isempty(varargin)
   v = {};
else
   [ok,v] = chkcstr(varargin);
   if ~ok
       error('Inputs must be Strings');
   end
   ii  = ( strcmp(v,'-r') | strcmp(v,'-R') );
   if any(ii)
      rec = v(max(find(ii)));
      rec = rec{end}(2);
      rec = ( rec == upper(rec) ); % True for follow Links
      rec = ( rec | ~isunix );     %  !!!
      v   = v(find(~ii));
   end
   if ~isempty(v)
       for jj = 1 : prod(size(opt))
           ii  = strcmp(v,opt{jj});
           def(jj) = any(ii);
           if def(jj)
              v   = v(find(~ii));
           end
           if isempty(v)
              break
           end
       end
   end
   if ~isempty(v)
       ii = strmatch('-b',v);
       if ~isempty(ii)
           bcp = v{max(ii)};
           v(ii) = [];
       end
   end      
   if ~isempty(v)
       ii = strmatch('-e',v);
       if ~isempty(ii)
           pref  = '';
           v(ii) = [];
       end
   end      
end

if isempty(v)
   v = {'.m'};
end

nv = prod(size(v));

tst = def(1);
new = def(2);
prv = def(3);
cls = def(4);

if tst
   fprintf(1,'!!!!!!  TestMode  !!!!!!\n');
end

%-----------------------------------------------------
% BackupDirectory

if ~isempty(bcp)
    if size(bcp,2) == 2
       bcp = '';
    else
       bcp = rmblank(bcp(3:end),2);
    end
    if isempty(bcp)
       fprintf(1,'%s',char(7));
       bcp = input('Enter DirectoryName for Backup: ','s');
       if isempty(bcp)
          fprintf(1,'Aborted\n');
          if Nout == 0
             clear c
          end
          return
       end
       bcp = rmblank(bcp,2);
    end
end

if isempty(bcp)
   if ~tst
       fprintf(1,'%s',char(7));
       wrn = 'Warning: NO Backup! Continue anyway? Yes/{No}: ';
       ok = input(wrn,'s');
       if ~strcmp(lower(ok),'yes')
          fprintf(1,'Aborted\n');
          if Nout == 0
             clear c
          end
          return
       end
   end
elseif ~( exist(bcp,'dir') == 7 )
   error(sprintf('BackupDirectory doesn''t exist: %s',bcp));
else
   % Get full PathName
   p0 = pwd;
   try
      cd(bcp);
      bcp = cd;
      cd(p0);
   end
   if ~( bcp(end) == fs );
         bcp = cat( 2 , bcp , fs );
   end
   [p,n] = fileparts(pfad(1:(end-1)));
   if isequal(p,bcp(1:(end-1)))
      error(sprintf('BackupDirectory and Path are the same: %s',pfad));
   end
end

%***************************************************************

%---------------------------------------------
% Get Files

fs = filesep;

d = getfiles(pfad,v,new,rec,prv,cls,'',fs,pref);

if isempty(d)
   if Nout == 0
      clear c
   end
   return
end

%---------------------------------------------
% Copy Files

c = cpfiles(d,bcp,c,tst);

if tst & ~isempty(c)
   fprintf(1,'!!!!!!  TestMode  !!!!!!\n');
end

if Nout == 0
   clear c
end

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function lst = getfiles(pfad,v,new,rec,prv,cls,pre,fs,pref)

lst = cell(1,3);

if ~( pfad(end) == fs );
    pfad = cat( 2 , pfad , fs );
end

[r,p] = fileparts(pfad(1:(end-1)));

pre   = cat(2,pre,fs,p);

lst{3} = pre;

%---------------------------------------------

c = cell(0,3);   % { OrgFile  NewFile Name }

for ii = 1 : prod(size(v))

    d = dir(cat(2,pfad,pref,v{ii}));

    if ~isempty(d)
        ok = ~cat(1,d.isdir); 
        if any(ok)
           nk = sum(ok);
           ok = find(ok);
           fn = {d(ok).name};
           % Check for ClassName
           [r,p]  = fileparts(pfad(1:(end-1)));
           if ~isempty(p)
               if ( p(1) == '@' ) & ( size(p,2) > 1 )
                  p = cat(2,p(2:end),fs);
               else
                  p = '';
               end
           end
           for ff = fn
               m = ff{1};
               n = which([p m]);
              ok = ~( isempty(n) | ...
                      isequal(n,'built-in') | ...
                      isequal(n,'variable') | ...
                      isequal(lower(m),'contents.m') );
               m = cat(2,pfad,m);
               if ok
                  % Check for equal FileNames but not equal FullNames
                  [p1,nn,e] = fileparts(n); nn = [ nn  e ];
                  [p1,nm,e] = fileparts(m); nm = [ nm  e ];
                  if isunix
                     ok = ( ~isequal(n,m) & isequal(nn,nm) );
                  else
                     ok = ( ~isequal(lower(n),lower(m)) & isequal(lower(nn),lower(nm)) );
                  end
               end
               if ok & new %, if strcmp(nm,'toposurf.m'),keyboard,end
                  d1 = dir(n);
                  d2 = dir(m);
                  ok = ( ( prod(size(d1)) == 1 ) & ( prod(size(d2)) == 1 ) );
                  if ok
                     ok = ~( d1.isdir | d2.isdir );
                     if ok
                        ok = ( datenum(d1.date) > datenum(d2.date) );
                     end
                  end
               end
               if ok
                  c = cat( 1 , c , { m  n  nm } );
               end
           end
        end
    end
 
end

lst{1} = c;

%---------------------------------------------
% Recurse

c = cell(0,3);

if ~isnan(rec)

     d = dir(pfad);

     if ~isempty(d)

         isd = cat(2,d.isdir);
         dnm = { d.name };

         d   =  char(dnm);
         l   = lower(dnm);

         isd = ( isd & ~( strcmp(dnm,'.') | strcmp(dnm,'..') | ...
                          ( strcmp(l,'private') & ~prv )     | ...
                          ( ( d(:,1)' == '@' )  & ~cls )            ) );

         if any(isd)

            isd = find(isd);

            dnm = dnm(isd);

            for d = dnm(:)'

                pf = cat( 2 , pfad , d{1} );

                ok = rec;

                if ~ok
                    ok = ~islink(pf);
                end

                if ok
                    n = getfiles(pf,v,new,rec,prv,cls,pre,fs,pref);
                    c = cat( 1 , c , n );
                end

            end

         end

    end

end   

lst{2} = c;

if isempty(lst{1}) & isempty(lst{2})
   lst = cell(0,3);
end

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function c = cpfiles(d,bcp,c,tst);

for ii = 1 : size(d,1);

    kk = isempty(bcp);

    msg = '';

    if ~kk

        dst = fullfile(bcp,d{ii,3});

        kk = ( exist(dst,'dir') == 7 );

        if kk

           fprintf(1,'Use BackupDirectory:  %s\n',dst);

        else

           fprintf(1,'Create BackupDirectory:  %s',dst);

           kk = tst;
           if ~kk
              [kk,msg] = mkdir(bcp,d{ii,3});
           end

           if ~kk
               fprintf(1,'  error');
           end
           if ~isempty(msg)
               fprintf(1,'\n%s',msg);
           end
           fprintf(1,'\n');

        end

    end

    if kk 

       if ~isempty(d{ii,1})

           m = d{ii,1}(:,1);  % Org
           n = d{ii,1}(:,2);  % New
           f = d{ii,1}(:,3);  % Name

           nn = size(d{ii,1},1);
           ok = ones(nn,1);

           for jj = 1 : nn

               if isempty(bcp)
                  f{jj} = '';
               else
                  f{jj} = fullfile(dst,f{jj});
                  fprintf(1,'Backup  %s  to  %s',m{jj},f{jj});
                  if ~tst
                      [ok(jj),msg] = copyfile(m{jj},f{jj});
                  end
                  if ~ok(jj)
                      fprintf(1,'  error');
                  end
                  if ~isempty(msg)
                      fprintf(1,'\n%s',msg);
                  end
                  fprintf(1,'\n');
               end

               if ok(jj)

                  fprintf(1,'Update  %s  by  %s',m{jj},n{jj});
                  if ~tst
                      [ok(jj),msg] = copyfile(n{jj},m{jj});
                  end
                  if ~ok(jj)
                      fprintf(1,'  error');
                  end
                  if ~isempty(msg)
                      fprintf(1,'\n%s',msg);
                  end
                  fprintf(1,'\n');
               end

           end

           if any(ok)
              cc = [ m  n  f ];
              if ~all(ok)
                  ok = find(ok);
                  cc = cc(ok,:);
              end
              c = cat( 1 , c , cc );
           end

       end

       if ~isempty(d{ii,2})
           c = cpfiles(d{ii,2},bcp,c,tst);
       end

    end

end

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

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

