function [Msg,str,v,d] = whosfile(file);

% WHOSFILE  returns Contents of MAT-File as CellString
%
% [Msg,CellText,VarStruct,FileStruct] = WHOSFILE( MatFileName )
%
%  VarStruct: name size bytes class
% FileStruct: name date bytes isdir creation
%
% see also: WHOS, LOOK_MAT
%

Nout = nargout;

Msg  = '';

 str = cell(0,1);
   v = [];
   d = [];

  nl = char(10);

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
Msg0 = sprintf('%s: ',upper(fcn));

%---------------------------------------------------------
% Check File

if ~( ischar(file) & ~isempty(file) & ...
      ( prod(size(file)) == size(file,2) ) );
   Msg = [ Msg0 'Input File must be a String.' ];
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
    Msg = [ Msg0 'File '  file ' does not exist.' ];
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
   Msg = [ Msg0 'Can''t read File '  file '.' ];
elseif ~( prod(size(d)) == 1 )
   Msg = [ Msg0 'Multiple read of File '  file '.' ];
elseif isequal(d.isdir,1)
   Msg = [ Msg0 'File '  file ' is a Directory.' ];
elseif isequal(d.bytes,0)
   Msg = [ Msg0 'File '  file ' is empty.' ];
end

if ~isempty(Msg)
    return
end

d.name = file;

d.creation = '';

%---------------------------------------------------------
% WHOS

try
  v = whos('-file',file);
catch 
  Msg = lasterr;
end

if ~isempty(Msg) | isempty(v)
  Msg = [ Msg0 'File '  file ' must be a MAT-File.' ...
            nl(1:(end*(~isempty(Msg))))  Msg ];
  return
end

%---------------------------------------------------------

bl = char(32);

N = prod(size(v));
   
n = str2mat( 'Name'  , bl , v.name  , bl , 'Summary' );  % Name
c = str2mat( 'Class' , bl , v.class , bl , bl        ); % Class

s = cell(N+4,1);   % Size
b = cell(N+4,1);   % Bytes

s{  1} = cat( 2 , bl(1,ones(1,8)) , 'Size' );
s{  2} = bl;
s{N+3} = bl;

b{  1} = cat( 2 , bl(1,ones(1,4)) , 'Bytes' );
b{  2} = bl;
b{N+3} = bl;

ss = 0;
bb = 0;

for ii = 1 : N

    s{ii+2} = sprintf('%9.0f',v(ii).size(1));

    for jj = 2 : size(v(ii).size,2);

       s{ii+2} = cat( 2 , s{ii+2} , ...
                     sprintf(' x %.0f',v(ii).size(jj)) );
   
    end

    b{ii+2} = sprintf('%9.0f',v(ii).bytes);

    bb = bb + v(ii).bytes;
    ss = ss + prod(v(ii).size);

    s{N+4} = sprintf('%13.0f',ss);
    b{N+4} = sprintf('%9.0f',bb);

end


n = cat( 2 ,      n  , bl(ones(N+4,1),ones(1,2)) );
s = cat( 2 , char(s) , bl(ones(N+4,1),ones(1,2)) );
b = cat( 2 , char(b) , bl(ones(N+4,1),ones(1,3)) );


str = cat(2,bl(ones(N+4,1),1),n,s,b,c);

str(  2,:) = '-';
str(N+3,:) = '-';

str = cellstr(str);

if nargout == 0

   clear Msg

   fprintf(1,'\n');
   fprintf(1,'%s\n',str{:});
   fprintf(1,'\n');

elseif Nout == 4

  d.creation = info(file);

end

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
