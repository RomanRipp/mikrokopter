function [msg,val,str,gd] = str2vec(str,varargin)   

% STR2VEC  Converts String into Numeric
%
%  [Msg,V] = STR2VEC( String )
%
%  [Msg,V] = STR2VEC( String , [@Class] , AllowedCharacters , NotAllowedCharacters )
%
%  Empty Value for AllowedCharacters: Allowed Characters will get automaticly
%                                         
%  Empty Value for NotAllowedCharacters: no Character is notallowed
%
%  If the Input '@ClassName' is given, the Values will check with this.
%
%--------------------------------------------------------------
%
%  [Msg,V] = STR2VEC( ['@' FileName '@'] , ... )
%  
%   Reads the String from a File, specified by FileName.
%
%--------------------------------------------------------------
% 
%  [Msg,V,String,AllowedCharacters] = STR2VEC( String , ... );
%


msg = '';
val = [];


Nin = nargin;

%*******************************************************

if Nin == 0
   msg = 'Input String is missing.';
   return
end

%------------------------------------------------------
% Check String

if isempty(str)
   str = '';
end

if ~chkstr(str,0)
   msg = 'Input must be a String.';
   return
end

if isempty(str)
   return
end

%-------------------------------------------------------
% Check GOOD, BAD, Class

cc    = cell(3,1);
cc(:) = { char(zeros(1,0)) };

nv = prod(size(varargin));
ok = zeros(nv,1);

cok = zeros(size(cc,1),1);

for ii = 1 : nv

    v = varargin{ii};

    ok(ii) = chkstr(v,0);

    if ok(ii)
       if isempty(v)
           for jj = [ 1  2 ]
               if ~cok(jj)
                   cok(jj) = 1;
                    cc{jj} = v;
               end
               if cok(jj)
                  break
               end
           end
       else
           for jj = [ 3  1  2 ]
               if ~cok(jj)
                   if strcmp(v(1),'@') | ~( jj == 3 )
                      cok(jj) = 1;
                       cc{jj} = v( (1+(jj==3)) : end );
                   end
                  if cok(jj)
                     break
                  end
              end
           end
       end
   end

end

if ~all(ok)
   msg = 'Following Inputs must be Strings';
   return
end

gd = cc{1};
bd = cc{2};
cl = cc{3};


%*******************************************************
% Check for FileName

n2 = size(str,2);

if ( n2 > 2 ) & strcmp(str([1 n2]),'@@') & ...
    ~any( ( str == 10 ) | ( str == 13 ) |  ...
          ( str ==  9 ) | ( str == 32 )         )

    try
        [msg,str] = loadfile(str(2:n2-1),5*2^20,'char');
    catch
         msg = sprintf('Error using LOADFILE.\n%s\n',lasterr);
    end

    if ~isempty(msg)
        if isempty(str)
           return
        end
        warning(msg);
    end

end

%*******************************************************

if isempty(str)
   return
end

n2 = size(str,2);

%*******************************************************
% Remove NotAllowedCharacters

if ~isempty(cc{2})
    ok = zeros(1,n2);
    for c = cc{2}
        ok = ( ok  |  ( str == c ) );
    end
    if any(ok)
           ok  = find(ok);
       str(ok) = char(32);
    else
           ok  = [];
    end
end

%*******************************************************
% Return if Class CHAR, 
%  remove following NotAllowedCharacters

if strcmp(cl,'char')
   if ~isempty(cc{2})
       if prod(size(ok)) > 1
              ok  = ok(find(diff(ok)==1)+1);
          str(ok) = [];
       end
   end
   val = str;
   return
end

%*******************************************************
% Save NewLineCharacters

is_nl = find( double(str) == 10 );  % NewLineCharacters

   nl = [ 10  32 ];                 % Order to Set NewLineCharacters

%*******************************************************
% Extract AllowedCharacters

is_nan = cat( 2 , findstr(str,'NaN') , findstr(str,'nan') );

low = lower(str);

if ~isempty(cc{1})

   %----------------------------------------------
   % Remove following IMAG and EXP

   if any( ( cc{1} == 'i' ) | ( cc{1} == 'j' ) | ( cc{1} == 'e' ) )
      ii = ( ( low == 'i' ) | ( low == 'j' ) | ( low == 'e' ) );
      if any(ii)
         ii = find(ii);
         jj = find( diff(ii) == 1 );
         str(ii([jj jj+1])) = char(32);
      end
   end

   %----------------------------------------------
   % Check for NaN, Remove following

   is_nan = cat( 2 , findstr(str,'NaN') , findstr(str,'nan') );

   is_nan = sort(is_nan,2);

   if ~isempty(is_nan)
      jj = find( diff(is_nan,1,2) == 3 );
      is_nan([jj jj+1]) = [];
   end

   %----------------------------------------------
   % Ok for allowed Characters

   ok = ( str == char(32) );

   for c = cc{1}
       ok = ( ok  |  ( str == c ) );
   end

   if ~all(ok)
       str(find(~ok)) = char(32);
   end

   %----------------------------------------------
   % Reset NaN

   if ~isempty(is_nan)

      n = size(is_nan,2);

      nn = 'NaN';
      nn = nn(ones(n,1),:);
      nn = permute(nn,[2 1]);
      nn = permute(nn(:),[2 1]);

      in = grp2ind(is_nan,3*ones(1,n));

      str(in) = nn;
       
       ok(in) = 1;

   end

   %----------------------------------------------
   % Special Handling for IMAG: Left AND Right must be Ok

   ii = ( ( low == 'i' ) | ( low == 'j' ) );
   ii = ( ii & ok );

   if any(ii)
      ii = find(ii);
      jj = find( ~( ok(ii-(ii>1)) & ok(ii+(ii<n2)) ) );
      str(ii(jj)) = char(32);
   end

   %----------------------------------------------
   % Special Handling for EXP: Left AND Right must be Ok
   % #. before EXP, +-# behind EXP, # = 0 .. 9 

   ii = ( ( low == 'e' ) & ok );

   if any(ii)

      ii = find(ii);

      jj = ii - ( ii > 1 ); % previous
      kk = ( ( '0' <= str(jj) ) & ( str(jj) <= '9' ) );
      kk = ~( ok(jj) & ( kk  | ( str(jj) == '.' ) ) );
      if any(kk)
         kk = find(kk);
         str(ii(kk)) = char(32);
      end

      jj = ii + ( ii < n2 ); % next
      kk = ( ( '0' <= str(jj) ) & ( str(jj) <= '9' ) );
      kk = ~( ok(jj) & ( kk  | ( str(jj) == '+' ) | ( str(jj) == '-' ) ) );
      if any(kk)
         kk = find(kk);
         str(ii(kk)) = char(32);
      end

   end

   %----------------------------------------------

   nl = nl( 1 : 1+(~any(double(cc{2})==10)) );

end

%*******************************************************
% Try EVAL(str)

for c = nl

    str(is_nl) = char(c);

    [msg,val] = evalstr(str);

    if isempty(msg)
       ok = chksize(val,str);
       if ~ok
           msg = 'Invalid Size of expression.';
       end
    end

    if isempty(msg) & ~isempty(cl)
       ok = strcmp( class(val) , cl );
       if ~ok
           [ok,val] = classchk(val,cl);
       end
       if ~ok
           msg = sprintf('Class "%s" is requested.',cl);
       end
    end
 
    if isempty(msg)
       return
    end

end

if ~isempty(cc{1})
    return
end

%*******************************************************
% Set GOOD Automaticly

good = { '0123456789.ije+-*/^()[]:,;~=' 
         '0123456789.ije+-*/^()[]:,;' 
         '0123456789.ije+-*/^()[]:' 
         '0123456789.ije+-*/^' 
         '0123456789.ije-+' 
         '0123456789.-' 
         '0123456789.' 
         '0123456789'            };
              
str(is_nl) = char(10);

for ii = 1 : size(good,1)

    [msg,val,str,gd] = str2vec(str,['@' cl],good{ii});

    if isempty(msg)
       return
    end

end


%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,val] = evalstr(str)

msg = '';
val = [];

ww = warnstat;
     warning('off')

lastwarn('')

try
   val = eval(cat(2,'[',str,']'));
   msg = lastwarn;
catch
   msg = lasterr;
end

warning(ww)

%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  ok = chksize(val,str);

ok = ( isempty(val) & isempty(str) );
if ok
   return
end

val = prod(size(val));

bl = ( ( str == char(32) ) | ( str == char(10) ) );

if ~any(bl) | all(bl)
    ok = ( ( ( val == 1 ) & ~any(bl) ) | ...
           ( ( val == 0 ) &  all(bl) )       );
    return
end

n2 = size(str,2);

%----------------------------------------------------
% Group of Blanks

[i0,lg] = ind2grp(find(bl));

n = size(i0,1);

i1 = i0(n) + lg(n) - 1;

i01 = ( [ i0(1)  i0(n)+lg(n)-1 ] == [ 1  n2 ] );

n = n - sum(i01);

ok = ( val == n+1 );

if ok
   return
end

%----------------------------------------------------
% Check for Sign +/- if Dot or Number follows
%  and no blank before

ii = ( ( str == '+' ) | ( str == '-' ) );

if any(ii)

   ii = find(ii);

   % Number follows
   jj = ii + ( ii < n2 ); % next
   kk = ( ( '0' <= str(jj) ) & ( str(jj) <= '9' ) );
   kk = ( kk  | ( str(jj) == '.' ) );
 
   % Blank before
   jj = ii - ( ii > 1 ); % previous
   kk = ( kk & bl(jj) );

   if any(kk)         % SignCharacter +/-
      kk = find(kk);
      str(ii(kk)) = char(32);
   end

end

%----------------------------------------------------
% Check for OperatorCharacters between Numbers

op = zeros(size(bl));

for c = '+-*/^'
    op = ( op | ( str == c ) );
end

if ~any(op)
    return
end

bl = bl + 2*op;  % Two if Operator

[i0,lg] = ind2grp(find(bl));

n = size(i0,1);

% Check Sum of Groups with Length 

sg = cumsum(bl,2);
ii = ( 1 : n-1 );
sg = sg(i0+lg-1) - cat( 2 , 0 , sg(i0(ii)+lg(ii)-1) );

op = ~( sg == lg(:)' );  % True if Operator in Group

i01 = ( [ i0(1)  i0(n)+lg(n)-1 ] == [ 1  n2 ] );

op([1 n]) = ( op([1 n]) & ~i01 );
 
n = n - sum(i01) - sum(op);

ok = ( val == n+1 );


%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ii,nn] = grp2ind(i0,l,s);

% GRP2IND  Built IndexVector from StartIndex and Length
%
% [Index,GroupNumber] = GRP2IND( StartIndex , GroupLength )
%
% where LENGTH(GroupNumber) == MAX(Index)
%
% GRP2IND( ... , LowSampleStep  ) LowSamples in Groups
%
% See also: IND2GRP, GRPMEAN, CUMSUM, CUMPROD
%

ii = [];
nn = [];

if isempty(i0);
   return
end

if nargin < 3
   s = 1;
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

jj = ( l == 0 );
if all(jj)
   return
end

if any(jj)
      jj  = find(jj);
   i0(jj) = [];
    l(jj) = [];
end

n = size(l,1);

l = ceil( l / s );

ii = s * ones(sum(l),1);
jj = cumsum( cat(1,1,l) , 1 );

ii(jj(1:n)) = i0;

if n > 1
   ii(jj(2:n)) = ii(jj(2:n))-(i0(1:n-1,1)+s*(l(1:n-1)-1));
end

ii = cumsum(ii,1);

if ( nargout == 2 ) & ~isempty(ii)
   nn = zeros(max(ii),1);
   kk = zeros(size(ii));
   kk(jj(1:n)) = 1;
   kk = cumsum(kk);
   nn(ii) = kk;
end

if n == 1
   return
end

jj = find( si > 1 );

if jj == 1
   return
end
  
perm = cat( 2 , (1:jj-1)+1 , 1 , (jj+1:size(si,2)) );

ii = permute(ii,perm);

if ( nargout == 2 )
   nn = permute(nn,perm);
end

%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [i0,l] = ind2grp(ii);

% IND2GRP  Built StartIndex and Length from IndexVector
%
% [ StartIndex , GroupLength ] = IND2GRP( Index )
%
% See also: GRP2IND, GRPMEAN, CUMSUM, CUMPROD
%

Nout = nargout;

i0 = zeros(0,1);
l  = zeros(0,1);

if isempty(ii);
   return
end

ii = ii(:);
n  = size(ii,1);

i0 = cat( 1 , 1 , find( diff(ii,1,1) > 1 )+1 , n+1 );

l  = diff(i0,1,1);

i0 = ii(i0(1:end-1));

%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,bb] = loadfile(file,varargin);

% LOADFILE  Load binary Data from File, using FREAD
%
% [ Msg , V ] = LOADFILE( FileName , MaxSize , Precision )
%
%   MaxSize     MaximumFileSize [Bytes],   default: 2097152 == 2MBytes
%   Precision   Precision for using FREAD, default: 'char'
%
%  For more Informations  type: >> help fread
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

m = 2*2^20; % 500000;  % MaximumFileSize [Bytes]
p = 'char';

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

    if chkstr(v,1)

       p = v;

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
     msg = [ msg0  'Can''t open File.' ];
     return
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

     ok = ( ( bb ==  9 ) |  ...
            ( bb == 10 ) |  ...
            ( bb == 13 ) |  ...
            (  28 <= bb  &   bb <= 126 ) | ...
            ( 160 <= bb  &   bb <= 255 )        );

     if ~all(ok)

         msg = [ msg0 'Invalid Characters in File.' ];

         if ~any(ok)
             bb = [];
         else
                ok  = find(~ok);
             bb(ok) = 32;
         end

     end

     bb = char(bb(:)');

  end

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ww = warnstat

% WARNSTAT  Returns global WarningStatus
%
%  WARNSTAT returns the Status of WARNING
%
% Matlab R<13   WARNING
% Matlab R>12   WARNING for Identifier ALL
%

ww = warning;

if isstruct(ww)   % New Matlab R>12 Syntax
   try
      id = strcmp({ww.identifier},'all');
      if any(id)
         id = find(id);
         ww = ww(id(1)).state;
      else
         ww = '';
      end
   catch
      ww = '';
   end
elseif ~chkstr(ww)
   ww = '';
end

%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

