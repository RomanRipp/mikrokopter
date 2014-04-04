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
