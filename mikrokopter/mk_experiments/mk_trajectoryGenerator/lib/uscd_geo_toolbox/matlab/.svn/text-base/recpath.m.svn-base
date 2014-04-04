function  [msg,pfade,mpath] = recpath(varargin)

% RECPATH  Adds/Removes Directory recursivly to Matlab's-SearchPath
%
% [Msg,Dirs,MatlabPath] = RECPATH( Directory , Options... )
%
% Adds or removes the Directory and recursivly the SubDirectories
%  to Matlab's-SearchPath. A missing or empty Input for Directory 
%   use the current working directory (PWD).
%
% Msg    contains a ErrorMessage, empty if all ok.
% Dirs   is the CellstringArray of the added Directories
%
% valid Options (multiple allowed):
%
%  '-begin'     '-b'  appends  the directories  (default)
%  '-end'       '-e'  prepends the directories
%
%  '-remove'    '-r'  removes  the directories
%
%  '-nofollow'  '-n'  don't follow DirectoryLinks
%     
%  '-display'   '-d'  display's added/removed Directories 
%                      in MatlabCommandWindow
%
%  '|name'            don't add SubDirectories which Names ends with "name"
%                      or matching "name" if it contains the wildcard "*"
%
%
% Note: "private"-Directories and Class-Directories starting with "@"
%        will never added to Matlab's-SearchPath
%
% see also: ADDPATH, REMPATH, PATH, MATLABPATH, STRWCMP
%


msg0    = 'RECPATH: ';
pfade   = {};
nl      = char(10);
tb      = char(32*ones(1,size(msg0,2)));  % TabSpace
fsp     = filesep;
psp     = pathsep;
comp    = computer;
mpath   = matlabpath;

wc      = '*';         % WildCard for DirNames to exclude

is_win = strcmp( upper(comp) , 'PCWIN' );

Nin  = nargin;
Nout = nargout;

pfad = '';

if Nin > 0
   v  = varargin{1};
   ok = isempty(v);
   if ~ok
       ok = chkstr(v,1);
       if ok
          ok = ~( v(1) == '-' );
       end
   end
   if ok
      pfad = v;
      varargin{1} = 'none';
   end
end

msg = cell(0,1);

%-------------------------------------
% Check Inputs

ok = ( isempty(pfad) | chkstr(pfad,0) );

if ~ok
    msg  = cat(1,msg,{'Pfad-Input must be a String'});
else

   if isempty(pfad)
      pfad = cd;
   else
      if is_win
         pfad = strrep(pfad,'/',fsp);  % Take care before UNIX-FileSeperator
      end
      ok = ( ( exist(pfad,'dir') == 7 ) | ( ~is_win  & strcmp(pfad,fsp) ) );
      if ok
         orgpfad = cd;
         try
            cd(pfad);
         catch
            ok = 0;
         end
         if ok
            pfad = cd;
            cd(orgpfad)
         end
      end
      if ~ok
          msg  = cat(1,msg,{['Invalid Directory: '  pfad ]});
      end
   end

   if ok & isempty(pfad)
      if is_win
         pfad = [ 'C:'  fsp ];
      else
         pfad = fsp;
      end
   end

end

if ~( chkcstr(varargin) | ( Nin <= 1 ) )
      msg = cat(1,msg,{'Option-Inputs must be strings.'});
end

%-------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    msg = sprintf('Invalid Inputs.\n%s',msg);
    if Nout == 0
       error(msg);
    end
    return
end

msg = '';

%-------------------------------------

pfad = cat( 2 , pfad , fsp(1:double(~strcmp(pfad(end),fsp))) );

remsep = ~( strcmp(pfad,fsp) | strcmp(comp([1 2]),'MA') );

remsep = double(remsep);

  pfad = pfad( 1 : (end-remsep) );


%-------------------------------------
% Check Options

opt  = 'b';  % Begin
mode = 0;
folw = 1;
excl = {};

if Nin > 0

   v = char(varargin);

   if ~isempty(v)

       %-------------------------------------------------------

       jj = ( v(:,1) == '-' );

       if any(jj)

          jj = find(jj);
           w = v(jj,min(2,size(v,2)));    % 2. Letter

          jj = ( ( w == 'b' ) | ( w == 'e' ) | ( w == 'r' ) );
          if any(jj)
             jj = max(find(jj));
             opt = w(jj);
          end

          mode = any( w == 'd' );

          folw = ~( isunix & any( w == 'n' ) );

       end

       %-------------------------------------------------------
       % DirNames to Exclude

       jj = ( v(:,1) == '|' );

       if any(jj)
          jj = find(jj);
           w = char(varargin(jj));
          if size(w,2) >= 2
              w = w(:,2:end);
             jj = ~( sum( w == ' ' , 2 ) == size(w,2) );
             if any(jj)
                if ~all(jj)
                    jj = find(jj);
                     w = w(jj,:);
                end
                w = cellstr(w);
                if is_win
                   w = strrep(w,'/',fsp);  % Take care before UNIX-FileSeperator
                end
                excl = w(:)';
                for ii = 1 : size(excl,2)
                    excl{ii} = cat( 2 , excl{ii} , fsp(1:double(~any(strcmp(excl{ii}(end),{fsp wc})))) );
                end
             end
          end
       end

   end

end

%-------------------------------------
%  GetPath

if mode
   fprintf(1,'\n%sGet DirectoryStructure of:\n\n%s%s\n',msg0,tb,pfad);
end

pfade = getpath0([pfad fsp(1:remsep)],remsep,fsp,folw,excl,wc);


%-------------------------------------
%  SetPath

if mode

   str = { 'Append'  'Prepend'  'Remove' };

   str = str{ 1 + ( opt == 'b' ) + 2 * ( opt == 'r' ) };

   txt = { 'to'  'from' };

   txt = txt{ 1 + ( opt == 'r' ) };

   fprintf( 1 , '\n%s%s Directories %s Matlab''s SearchPath:\n\n%s%s\n' , ...
                 msg0 , str , txt , tb , strhcat(pfade,'',1,[nl tb]) );

end

mpath = setpath(pfade,psp,opt,is_win);

if mode
   fprintf(1,'\n');
end

if Nout == 0
   clear msg
end

%********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  mpath = setpath(pfade,psp,Option,is_win);

% SETPATH  Set MatlabPath
%
%  MatlabPath = SETPATH( pfade , PathSep , Option , IsWin );
%

mpath = matlabpath;  % Current MatlabPath

%-------------------------------------------------
% Append and Prepend PathSep

mpath = cat(  2   , psp(1:(end*double(~strcmp(mpath( 1 ),psp)))) , ...
            mpath , psp(1:(end*double(~strcmp(mpath(end),psp)))) );

%-------------------------------------------------
% Check for existing Path's in current MatlabPath

ind = zeros(0,2);  % [ Start  Lenght ] of existing Path's

mp = mpath;
pf = pfade;

if is_win
   mp = lower(mp);  % !!!
   pf = lower(pf);  % !!!
end

for p = pf(:)'

    jj = findstr( mp , cat( 2 , psp , p{1} , psp ) );

    if ~isempty(jj)
       for kk = jj(:)'
           ind = cat( 1 , ind , [ kk  size(p{1},2)+1 ] );
       end
    end

end

%-------------------------------------------------
% Remove Path's from current MatlabPath

if ~isempty(ind)

   ii = grp2ind(ind(:,1),ind(:,2));

   mpath(ii) = [];

end

%-------------------------------------------------
% Append or Prepend Pfade

pfade = strhcat( pfade , psp , size(pfade,1)+1 );

switch Option

  case 'b'

     mpath = cat( 2 , psp , pfade , mpath );
    
  case 'e'

     mpath = cat( 2 , mpath , pfade , psp );

end

%-------------------------------------------------
% Remove PathSep from begin and End

mpath = mpath( 2 : ( end-1) );

%-------------------------------------------------
% Set new MatlabPath

matlabpath( mpath );


%********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function pfade = getpath0(pfad,remsep,fsp,folw,excl,wc);

% function pfade = getpath0(pfad,RemoveSeparator,FileSeperator,Follow);

pfade   = {};

%-------------------------------------
%  Recurse

d = dir(pfad);

if isempty(d)
   return
end

is_dir = cat(1,d.isdir);

if ~any(is_dir)
    return
end

is_dir = find(is_dir);

 pfade    = cell(size(is_dir(:),1),1);
 pfade(:) = { {} };

for ii = 1 : size(pfade,1)

    jj = is_dir(ii);

    if ~any(strcmp(d(jj).name,{ '.'  '..'  'private' }))  & ...
       ~strcmp(d(jj).name(1),'@')

        p = cat( 2 , pfad , d(jj).name , fsp(1:double(~strcmp(d(jj).name(end),fsp))) );

        ok = folw;
        if ~ok
            ok = ~islink(p);
        end

        if ok & ~isempty(excl)
           sp = size(p,2);
           for e = excl
               if any( e{1} == wc )
                  ok = ~strwcmp(p,e{1});
               else
                  se = size(e{1},2);
                  ok = ( se <= sp );
                  if ok
                     ok = ~strcmp( p((sp-se+1):end) , e{1} );
                  end
               end
               if ~ok
                   break
               end
           end
        end

        if ok
           pfade{ii}= getpath0( p , remsep , fsp , folw , excl , wc );
        end

    end

end


pfade = cat( 1 , {pfad(1:(end-remsep))} , cat(1,pfade{:}) );


%********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


if ischar(str)
  str = cellstr(str);
end

str = str(:);

str( : , 2 ) = { del };

str(n:n:size(str,1),2) = {nl};

str(    size(str,1),2) = { '' };


str = permute( str , [ 2  1 ] );

str = cat(2,str{:});


%********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,cmp,wc,nn,cc] = strwcmp(str,cmp,wc,nn,cc)

% STRWCMP   Compare Strings, including WildCards
%
% OK = STRWCMP( STR , CMP , WildCard )
%
% STR: CharacterArray or CellStringArray to compare with Comp
%
% CMP: CharacterArray or CellStringArray with strings to compare with STR
%         strings can contains WildCards
%
% WildCard: specify WildCard to use, default: '*'
%
% OK : logical Array with same size of STR, contains 1 if 
%        any strings of CMP match string of STR
%
% Special Output:  [ok,cmp,wc,nn,cc] = STRWCMP( STR , CMP , [WildCard] )
%
%  to use in follwing statements: ok = STRWCMP( STR , cmp , wc , nn , cc );
%
%  which makes it a bit faster.
%
% see also: STRCMP, FINDSTR
%

Nin  = nargin;
Nout = nargout;

Msg = '';
 nl = char(10);

%***************************************************************
% Check Inputs

%---------------------------------------------------------------
% String

if Nin < 1
   str = cell(0,0);
end

if ~( iscell(str) & isempty(str) )
   [ok,str] = chkcstr(str,0);
   if ~ok
      Msg=[ Msg  nl(1:(end*(~isempty(Msg)))) ...
           'First Input must be a CharacterArray or CellStringArray.' ];
   end
end

%---------------------------------------------------------------
% CompareString

if Nin < 2
   cmp = cell(0,0);
end

if ~( iscell(cmp) & isempty(cmp) )
   [ok,cmp] = chkcstr(cmp,0);
   if ~ok
      Msg=[ Msg  nl(1:(end*(~isempty(Msg)))) ...
            'Second Input must be a CharacterArray or CellStringArray.' ];
   end
end

%---------------------------------------------------------------
% WildCard

if Nin < 3
   wc = '*';
elseif ~( ischar(wc) & ( prod(size(wc)) == 1 ) )
   Msg=[ Msg  nl(1:(end*(~isempty(Msg)))) ...
         'WildCard must be a Single Character.' ];
end
  
%---------------------------------------------------------------

if ~isempty(Msg)
   error(Msg)
end

%***************************************************************

si = size(str);

if ( isempty(str) | isempty(cmp) ) & ( Nout <= 1 ) 
   ok = zeros(si);
   return
end

cmp = cmp(:);

if ~isempty(cmp)
    if any(strcmp(cmp,wc)) & ( Nout <= 1 )  % Wildcard only
       ok = ones(si);
       return
    end
end

%***************************************************************
% Analyze CompareStrings

nc = size(cmp,1);

ok = ( Nin == 5 );

if ok
   ok = ( isequal(size(cc),[nc 1]) & iscell(cc) & ...
          isequal(size(nn),[nc 2]) & isnumeric(nn)    );
   if ok
      try
         ok = cat( 1 , cc{:} ); 
         ok = ( size(ok,2) == 3 );
      catch
         ok = 0;
      end
   end
end


%--------------------------------------------------
if ~ok
%--------------------------------------------------

  cc    = cell(nc,1);
  cc(:) = { zeros(0,3) };  % { [ Start End N  ] }

  nn = zeros(nc,2);        %   [ Ncmp  sum(N) ] 

  for ii = 1 : nc

    if ~isempty(cmp{ii})

       iwc = ( double(cmp{ii}) == double(wc) );

       if ~any( iwc )

          nn(ii,:) = size(cmp{ii},2);
          cc{ii}   = [ 1  nn(ii,:) ];

       else

         %--------------------------------------------
         % Remove Duplicate WildCards

         iwc = find( iwc );
         if ~( prod(size(iwc)) == 1 )
            jj = find( diff(iwc) == 1 );
            if ~isempty(jj)
               cmp{ii}(iwc(jj+1)) = [];
            end
         end

         %--------------------------------------------
         % Get Start End
     
         iwc = ( double(cmp{ii}) == double(wc) );
  
          n  = size(iwc,2);

          if ( n == 1 ) & ( iwc == 1 ) & ( Nout <= 1 ) % Wildcard only
             ok = ones(si);
             return
          end

          i0 = ~iwc(1);
          i1 = ~iwc(n);

         iwc = cat( 2 , ones(1,i0) , iwc , ones(1,i1) );

         iwc = find( iwc );

         iwc = iwc(:);

           n = size(iwc,1) - 1;

         cc{ii} = zeros(n,3);

         if n > 0      

            cc{ii}(:,[1 2]) = cat( 2 , iwc(1:n)+1 , iwc((1:n)+1)-1 ) - i0;

            cc{ii}(:,3) = cc{ii}(:,2) - cc{ii}(:,1) + 1;
 
         end

         nn(ii,:) = cat( 2 , size(cmp{ii},2) , sum(cc{ii}(:,3),1) );

       end

    end

  end

%--------------------------------------------------
end
%--------------------------------------------------

if ( Nout > 1 )

  if ( isempty(str) | isempty(cmp) )
     ok = zeros(si);
     return
  end

  if any(strcmp(cmp,wc))  % Wildcard only
     ok = ones(si);
     return
  end

end

%***************************************************************
% Compare

ok = zeros(si);

for ii = 1 : prod(si)
 
    s2 = size(str{ii},2);

    for jj = 1 : nc

        ok(ii) = ( ( s2 == 0 ) & ( nn(jj,1) == 0 ) );
       
        if ok(ii)
           break
        end
       
        if ( s2 >= nn(jj,2) ) & ~( nn(jj,1) == 0 )

           ok(ii) = compare(str{ii},cmp{jj},cc{jj});

           if ok(ii)
              break
           end
            
        end

    end

end

%********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = compare( str , cmp , cc )

sc = size(cmp,2);

ok = 1;

for ii = 1 : size(cc,1) 

    s2 = size(str,2);

    ok = ( ok &  ( s2 >= cc(ii,3) ) );

    if ~ok
       break
    end

    ic  = ( cc(ii,1) : cc(ii,2) );

    i01 = ( cc(ii,[1 2]) == [ 1  sc ] );

    i0  = ( 1 : cc(ii,3) );
    i1  = s2 - cc(ii,3) + i0;

    i2  = cat( 2 , findstr( str , cmp(ic) ) , 0 );
    i2 = i2(1);
    i2 = ( i2 + cc(ii,3) - 1 ) * ~( i2 == 0 );

    i2  = s2       * (  i01(1) &  i01(2) & strcmp(str    ,cmp    ) ) + ...
          cc(ii,3) * (  i01(1) & ~i01(2) & strcmp(str(i0),cmp(ic)) ) + ...
          s2       * ( ~i01(1) &  i01(2) & strcmp(str(i1),cmp(ic)) ) + ...
          i2       * ( ~i01(1) & ~i01(2) );

    ok = ( ok & ~( i2 == 0 ) );

    if ~ok
       break
    end

    str = str( i2+1 : s2 );

end

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

