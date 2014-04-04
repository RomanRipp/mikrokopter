function out = bgrep(varargin)

% BGREP  Search in binary Files for specific string
% 
% bgrep(String,Pfad,'-r','-l','-d','-n','-K#','.ext1','.ext2', ... )
% bgrep(String,File,'-r','-l','-d','-n','-K#','.ext1','.ext2', ... )
%
%  '-r'  ==  recurse
%  '-l'  ==  follow DirectoryLinks
%  '-d'  ==  no Display
%  '-n'  ==  no IntroDisplay
%  '-B#' ==  Max. Number of Bytes to read from begin of File
%            neg. Number of Bytes to read from end   of File
%  '-K#' ==  Max. Size of Files in KBytes 
%  '-M#' ==  Max. Size of Files in MBytes 
%
%  defaults:   no recurse 
%              no follow DirectoryLinks
%              Display
%              IntroDisplay
%             '.ext' == '.m'
%             '-B#'  == '-BInf'   % Read full File
%             '-K#'  == '-K1024'  % 1-MBytes MaxSize
%             '-M#'  == '-M1'     % 1-MBytes MaxSize
%
%  note:  give '.'  to grep all Extensions
%    
%  out = { FileName  Date  Bytes  ByteNumber  String }
%


out = cell(0,5);

%****************************************************************************
% Get and Check Inputs

try
   [msg,str,Pfad,ext,recurse,follink,nodispl,nointro,bt,kb] = checkin(varargin{:});
catch
   msg = sprintf('Error call BGREP:CHECKIN\n%s',lasterr);
end
 
if ~isempty(msg)
   error(msg);
end

%****************************************************************************

fs = filesep;
nl = char(10);
wc = '*';

 
if ~nointro & ~nodispl

  txt = 'BGREP: ';

  if ~recurse
    txt = [ txt 'no ' ];
  end

    txt = [ txt 'recurse, extension: ' ];
  
  if isempty(ext{1})
    txt =  [ txt '''*'', ' ];
  else
    txt = [ txt  ''''  strhcat(ext,''', ''')  ''', ' ];
  end   

  if isfinite(bt)
     txt = [ txt sprintf(' read %.0f Bytes, ',bt) ];
  end

    txt = [ txt ' MaxSize: ' sprintf('%g ',kb) 'KBytes' ];

    fprintf(1,[ nl  txt nl nl ]);

end

if ~nodispl
   fprintf(1,'%s',Pfad);
end

d = dir(Pfad);

if isempty(d)
   if ~nodispl
       fprintf(1,' ... can''t read File or Directory\n');
   end
   return
end

%---------------------------------------------------------
% Check for File

IsFile = ( ( size(d,1) == 1 )  &  ( d(1).isdir == 0 ) );

if IsFile

   if isempty(fileparts(d(1).name))
      d(1).name = fullfile(fileparts(Pfad),d(1).name);
   end

   ext = { char(ones(1,0)) };
   wc  = ext{1};

else

    Pfad = cat( 2 , Pfad , fs(1:(end*(~strcmp(Pfad(end),fs)))) );
  
    is_dir = find(cat(1,d.isdir));

    dd = { d(is_dir).name };

end

if ~nodispl
    fprintf(1,nl)
end

ns = size(str,2);

bytes = kb * 1024;

if ~isfinite(bt)
    bt = Inf;
end

for ee = ext(:)'

  if ~IsFile
      d = dir(cat(2,Pfad,wc,ee{1}));
  end

  is_file = find(~cat(1,d.isdir));

  for ii = is_file(:)'

     if d(ii).bytes <= bytes

        name = d(ii).name;

        file = [ Pfad((1:(end*(~IsFile))))  name ];

        fid = fopen(file,'r');

        if fid == -1

          if ~nodispl
              fprintf(1,'   Error read %s\n',file);
          end

          bb = [];

        else

          if bt < 0
             fseek(fid,max(-d(ii).bytes,bt),'eof');
          end

          bb = fread(fid,[1 abs(bt)],'char');

          fclose(fid);

        end


        if ~isempty(bb)

          bb = char(bb); % num2char( bb , 32 );

          nn = cat(2,0,find(bb==32),size(bb,2)+1);

            if ns <= size(bb,2)

              ff = findstr(bb,str);
              for jj = ff(:)'

                 i0 = sum( nn <= jj+0    );
                 i1 = sum( nn <= jj+ns-1 ) + 1;

                 i0 = nn(i0) + 1;
                 i1 = nn(i1) - 1;

                 i0 = max( i0 , jj-64);
                 i1 = min( i1 , jj+ns-1+32 );

                 ll = bb(i0:i1);

                 i0 = i0 + max( d(ii).bytes + bt , 0 ) * ( bt < 0 );

                out = cat(1,out,{ file  d(ii).date  d(ii).bytes jj ll}); %%% i0  ll });

                if ~nodispl
                   fprintf(1,'%s: %8.0f: %s\n',file,i0,ll);
                end

              end
              % jj
            end
            % ns ok
 
        end 
        % fid == -1
     end
     % bytes ok
  end
  % ii
end
% ee

if ~recurse
    return
end

opt = { '-l' '-r' '-n' sprintf('-B%.0f',bt) sprintf('-K%.5f',kb) '-d' };

opt = opt( 2-follink : end-(~nodispl) );

is_dir = find(cat(1,d.isdir));

for d = dd(:)'

    if ~any(strcmp(d{1},{ '.'  '..' }))

      pf = cat( 2 , Pfad , d{1} );

      ok = follink;

      if ~ok
         ok = ~islink(pf);
      end

      if ok
         out1 = bgrep(str,pf,opt{:},ext{:});
          out = cat(1,out,out1);
      end

    end

end

   

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,str,Pfad,ext,recurse,follink,nodispl,nointro,bt,kb] = ...
            checkin(str,Pfad,varargin);


% CHECKIN  Check of Inputs 


msg = {};
nl  = char(10);

fs  = filesep;

Nin = nargin;

%*********************************************************
% Defaults

recurse = 0;
follink = 0;
nodispl = 0;
nointro = 0;

bt      = inf;

kb      = 1024;

ext  = cell(1,0);

ext0 = '.m';

%*********************************************************
% Get and Check Inputs

%---------------------------------------------------------
% Check Pfad and String

if Nin < 2
   Pfad = '';
end

if Nin < 1
   str = '';
   msg = 'Input String is missing.';
   return
end

if isempty(Pfad)
   Pfad = cd;
end

if ~chkstr(Pfad,0)
    msg = cat( 1 , msg , {'Pfad must be a String.'} );
elseif isempty(Pfad)  |  any(strcmp(Pfad,{ '.'  '*' }))
    Pfad = cd;
end

if ~chkstr(str,0)
    msg = cat( 1 , msg , {'String must be a NonEmpty String.'} );
end



%---------------------------------------------------------
% Get Options from VarArg


if Nin <  3

  ext = { ext0 };

else

  for ii = 1 : Nin-2

      v = varargin{ii};

      ok = chkstr(v,1);
      if ok
         s  = size(v,2);
         ok = ( ( ( v(1) == '-' ) & ( s >= 2 ) ) | ( v(1) == '.' ) ); 
      end

      if ok

         %--------------------------------------------------------

          recurse = ( recurse | strcmp(v,'-r') );
          follink = ( follink | strcmp(v,'-l') );
          nodispl = ( nodispl | strcmp(v,'-d') );
          nointro = ( nointro | strcmp(v,'-n') );

          %--------------------------------------------------------

          if strcmp( v(1) , '.' );
             ext = cat(2,ext,{v});
          end

          %--------------------------------------------------------
          if ( v(1) == '-' ) & ( s > 2 )

             w = v(3:end);
             v = v(2);

             if any( v == 'BKM' )

                m = '';
                eval(['w = ' w ';'],'m = lasterr;');

                if isempty(m)
                   ok = ( isnumeric(w) & ( prod(size(w)) == 1 ) );
                   if ok
                      ok = ( ( w > 0 ) | ( v == 'B' ) );
                   end
                   if ~ok
                       m = 'Value must be single positive finite.';
                   end
                end

                if isempty(m)
                   switch v
                       case 'B'
                          bt = ceil(w);
                       case 'K'
                          kb = w;
                       case 'M'
                          kb = w * 1024;
                   end
                else
                   m = sprintf('Invalid Value of %.0f. Input for Option "%s".\n%s',ii+2,v,m);
                   msg = cat(1,msg,{m});
                end

             end

          end

      else

          msg = cat(1,msg,{sprintf('%.0f Input must be a nonempty String.',ii+2)});

      end 
      % if char

  end
  % ii

  %---------------------------------------------------------------------------
  if isempty(ext)

     ext = { ext0 };

  elseif any( strcmp( ext , '.' ) )

     ext = { char(ones(1,0)) };

  end

end
% if Nin

%---------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
else
    msg = '';
end

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function c = num2char(c,r)

% NUM2CHAR   Converts numeric Array's to CharacterArrays
%
%  C = NUM2CHAR( Array , Replace );
%
%  default: Replace == ' ';  % SpaceCharacter
%
%  The values of Array will transformed into Interval [ 0 .. 255 ]
%  Bad Characters will replaced.
%
%  Good Characters are:  9, 10, 13, [ 28 .. 126 ], [ 160 .. 255 ]
%

 
Nin = nargin;
Msg = '';
 nl = char(10);

if Nin == 0
   c = char(zeros(0,0));
   return
end

%---------------------------------------------------------------

if ~isnumeric(c) | ischar(c)
   Msg = 'Input must be a Numeric- or Character-Array.';
end

%---------------------------------------------------------------

if Nin < 2
   r = 32;
elseif isempty(r)
   r = [];
else
   if ~( ( isnumeric(r) | ischar(r) ) & ( prod(size(r)) == 1 ) )
      Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
              'Value for REPLACE must be empty or a single Numeric or Character.' ];
   end
end

%---------------------------------------------------------------

if ~isempty(Msg)
   error(Msg)
end

%---------------------------------------------------------------

if isempty(c)
   c = char( zeros(size(c)) );
   return
end


  if ~isa(c,'double')
     c = double(c);
  end

c = floor( c - 256 * floor( c / 256 ) );


if ~isempty(r)

   if ~isa(r,'double')
      r = double(r);
   end

   r = floor( r - 256 * floor( r / 256 ) );

end


c( find( ~( ( c ==  9 ) |  ...
            ( c == 10 ) |  ...
            ( c == 13 ) |  ...
            (  28 <= c  &   c <= 126 ) | ...
            ( 160 <= c  &   c <= 255 )        ) ) ) = r;


c = char(c);

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
