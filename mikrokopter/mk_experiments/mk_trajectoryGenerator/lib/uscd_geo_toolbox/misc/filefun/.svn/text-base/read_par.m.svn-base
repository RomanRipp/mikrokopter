function [msg,ini,str] = read_par(file,sp,cl,cm);

% READ_PAR  Read ParameterFile
%
% [Msg,V] = READ_PAR( ParameterFile   ,  Seperator , ContinueMarker , CommentMarker )
% [Msg,V] = READ_PAR( ParameterString ,  Seperator , ContinueMarker , CommentMarker )
%
% Load Parameter-Value-Pairs from an ParameterFile.
% 
% The Parameter-Value-Pairs are seperated by Seperator.
%   A Value is valid until End of Line (NL), 
%     i.e. a Line can contain maximum ONE Parameter-Value-Pair.
%
%   Following Lines, beginning with the ContinueMarker, will added to the Value.
%
%   Seperator can be a 2-Element CellString: { Seperator ParameterMarker }
%
%   If Seperator is empty, only the Parameter, start with ParameterMarker,
%     will read, the ContinueMarker has no effect in that case.
%
% The Comment's are marked in the File with the CommentMarker.
%   A Comment is valid until End of Line (NL).
%
%
% Use "" around a Marker to mask them.
%
% Defaults:
%
%    Seperator      =   { '='  '' }   (no ParameterMarker)
%    ContinueMarker =   '++'
%    CommentMarker  =   '$$'
%
%
%----------------------------------------------------------------------------
% OUTPUT
%
%  Msg contains the ErrorMessage. Msg is empty, if READ_PAR works successful.
%
%  V   returns an [ N by 2 ]-CellArray for N Parameter: { Parameter  Value }
%
% [Msg,V,String] = READ_PAR( ... )
%
%  Returns the String without Parameter-Value-Pairs and Comments.
%
%----------------------------------------------------------------------------
%
% see also:  READ_OBJ, LOAD_OBJ, STR2PAR, LOADFILE
%
%



Nin  = nargin;
Nout = nargout-1;

msg = cell(0,1);

nl  = char(10);

ini = cell(0,2);  % { Parameter Value }
str = '';

quota = '"';

ext = '.par';  % DefaultExtension to use to expand File

%*********************************************************************
% Get Inputs

if Nin < 1
   msg = 'Input File is missing.';
   try
      msg = catmsg(msg,'',i);
   end
   return
end

if Nin < 2
 sp = { '=' '' };
end
if Nin < 3
 cl = '++';
end
if Nin < 4
 cm = '$$';
end

%*********************************************************************
% Check Inputs

%---------------------------------------------------------------------
% File

is_file = 0;

s1 = size(file,1);

if chkcstr(file,1)
   str = sprintf('%s\n',file{:});
elseif ~( ischar(file) & ( ndims(file) == 2 ) )
   msg = cat( 1 , msg , {'Input File must be a CharacterArray or CellArray of Strings.'} );
else
  if isempty(file)
     str = '';
  elseif ( s1 > 1 ) | any( file == 10 ) | any( file == 13 )
     str = file;
     if s1 > 1
        str = cat( 2 , str , char(10*ones(s1,1)) );
        str = permute(str,[2 1]);
        str = str(:);
        str = permute(str,[2 1]);
     end
  else
     is_file = 1;
  end
end
   
%---------------------------------------------------------------------
% Seperator

[ok,sp] = chkcstr(sp);
if ok
   sp = sp(:);
   sp = cat( 1 , sp , {''} );
   sp = sp([1 2]);
   st = rmblank(sp{2},2);    % ParameterStart
   sp = rmblank(sp{1},2);    % ParameterSeperator
   ok = ~( isempty(st) & isempty(sp) );
end

if ~ok
   msg = cat( 1 , msg , ...
              {'Input Seperator must be a CharacterArray or CellArray of Strings.'} , ...
              {' 2 Seperator defines { ParameterSeperator ParameterStart }.' } , ...
              {' where any of them must be nonempty.' }  );
end

%---------------------------------------------------------------------
% ContinueMarker

if ~chkstr(cl)
   msg = cat( 1 , msg , {'Input ContinueMarker must be a String.'} );
end

%---------------------------------------------------------------------
% CommentMarker

if ~chkstr(cm)
   msg = cat( 1 , msg , {'Input CommentMarker must be a String.'} );
end

%---------------------------------------------------------------------

if ~isempty(msg)
   try
      msg = catmsg( 'Invalid Inputs,' , msg , i );
   catch
      msg = sprintf('%s\n','Invalid Inputs.',msg{:});
   end
   return
end

%*********************************************************************
% Read File

if is_file

   [msg,str] = loadfile(file,20*2^20,'char',ext);

   if ~isempty(msg)
       try
          msg = catmsg( 'Invalid ParameterFile.' , msg , i );
       catch
          msg = sprintf('%s\n%s','Invalid ParameterFile.',msg);
       end
       if isempty(str)
          return
       end
       warning(msg);
   end

end

msg = '';

if ( Nout == 0 ) | isempty(str)
   return
end

str = strrep(str,char([13 10]),nl);  % DOS: CRLF --> LF  
str = strrep(str,char(13),nl);       % MAC: CR   --> LF

%*********************************************************************
% Mask

msk = { st  char(1)
        sp  char(2)
        cl  char(3)
        cm  char(4)
      quota char(5)  };

str = mask(str,msk,quota);

%*********************************************************************
% Initialisation for KeyWords
%
% { Parameter Value }
%

[ini,im,fm,lm,nm] = get_par(str,st,sp,cl,cm,nl,msk);

ret = ( ( Nout < 2 ) | isempty(im) );
if ~ret
    ret = all( fm == nm );
end

if ret
   str = mask(str,msk);
   return
end

%*********************************************************************
% Remove Property, Value, Comment

nn = size(im,2);

isn = ( fm == nm );

ok   = ~isn;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Use NotCommentKeys after NewLine too

   ind  = ( 1 : (nn-1) );
ok(ind) = ( ok(ind) | ( ok(ind+1) & ~( fm(ind+1) == (nm-1) ) & ( fm(ind) == nm ) ) );

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ii = find(ok);

%------------------------------------------
% StartIndex & GroupLength
% !!! FM starts with ZERO !!!

i0 = im(ii);
lg = im(ii+1) - im(ii) + lm(fm(ii+1)+1) .* ~( fm(ii+1) == nm );

%------------------------------------------
% Check StartIndex

dd = ( i0 - 1 ) .* ( i0 < 1 );
i0 = i0 - dd;
lg = lg + dd;

%------------------------------------------
% Check EndIndex

ns     = size(str,2);
   jj  = find( i0 > ns );
i0(jj) = [];
lg(jj) = [];

i1 = i0 + lg - 1;
dd = ( ns - i1 ) .* ( ns < i1 );
lg = lg + dd;

%------------------------------------------
% Group ---> Index

ind = grp2ind(i0,lg);

str(ind) = [];

str = mask(str,msk);

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ini,im,fm,lm,nm] = get_par(str,st,sp,cl,cm,nl,msk)

ini = cell(0,2);  % { Parameter Value }

%---------------------------------------------------------------------
% Get MarkerIndize, Flag and Length

no_sep = isempty(sp);   % No Seperator
if no_sep
   cl = '';
end

mark = { st sp cl cm nl };

nm = size(mark,2) - 1;         % !!! Start with ZERO !!!

[im,fm,lm] = get_mark( str , mark , 2-no_sep );

if isempty(im)
   return
end

fm = fm - 1;                   % !!! Start with ZERO !!!

%---------------------------------------------------------------------
% NewLine before Start and End

im = cat( 2 , 0  , im , size(str,2)+1 );
fm = cat( 2 , nm , fm , nm            );
  
%---------------------------------------------------------------------
% Check correct Statement of Marker

%       Flag FlagBefore True/False

cp = nm * ( lm(1) == 0 );      % Check for Seperator: ParStart | NewLine

chk = [  0    nm         1     % NewLine before Start
%%%      0    1         -1     % Seperator after Start
         1    cp         1     % ParStart | NewLine before Seperator
         2    nm         1     % NewLine before ContinueMarker
         3    3          0  ]; % No Comment before Comment


[im,fm] = chk_mark(im,fm,chk,0);

ok = ~isempty(im);
if ok
   ok = any( fm == ( 1 - no_sep ) );
end
if ~ok
    return
end

%---------------------------------------------------------------------
% Get Parameter-Value

[ini,i0,ok] = get_str( str , im , fm , lm , msk , nl , ...
                       [1 2]+[-1 1]*no_sep , ~no_sep + i*no_sep );

if any(~ok)
   jj  = find(~ok);
   if ~( lm(1) == 0 )
       jj = cat(1,jj,jj-1);  % Seperator and Start
   end
      jj  = i0(jj);
   im(jj) = [];
   fm(jj) = [];
end

%---------------------------------------------------------------------
% Check for correct Start

if ~any( lm([1 2]) == 0 ) 
    [im,fm] = chk_mark(im,fm,[0 1 -1],0);  % Seperator after Start
end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [im,fm,lm] = get_mark(str,mark,ret);

% GET_MARK  Returns MarkerIndize, Flag and Length
%

nm = size(mark,2);

im = zeros( 1 , 0  );  % Index for Marker
fm = zeros( 1 , 0  );  % Flag  for Marker
lm = zeros( 1 , nm );  % Length of Marker

%---------------------------------------------------------------------
% Find Marker: Index & Flag

for ii = 1 : nm

    lm(ii) = size( mark{ii} , 2 );

    if ~( lm(ii) == 0 )

        jj = findstr( str , mark{ii} );

        if ~isempty(jj)

            im = cat( 2 , im ,      jj );
            fm = cat( 2 , fm , ii+0*jj );
    
        elseif any( ii == ret )

            return
 
        end

    end

end

%---------------------------------------------------------------------
% Sort

[im,si] = sort(im);

fm = fm(si);

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [im,fm,cc] = chk_mark(im,fm,chk,mode);

% CHK_MARK  Check Correct Statement of Marker

% chk = [ Flag   FlagBefore   True/False ]
%                FlagBehind  -True       ]

cc = ones(size(im));

for ii = 1 : size(chk,1)

    jj = find( fm == chk(ii,1) );

    if ~isempty(jj)

        ck = ~( chk(ii,3) == 0 );         % True or False

        cj = 1 - 2 * ( chk(ii,3) >= 0 );  % Before or After

        ok = ( fm(jj+cj) == chk(ii,2) );

        ok = ( 1 - ck ) + ( 2*ck - 1 ) * ok;

        ok = ~ok;

        if any(ok)

           ok = find(ok);

           jj = jj(ok);

           cc(jj) = 0;

           if mode

              fm(jj) = [];
              im(jj) = [];

           else

              fm(jj) = NaN;

           end

        end

    end

end

if mode
   return
end

jj = find(~cc);

im(jj) = [];
fm(jj) = [];


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ini,ik0,ok] = get_str(str,im,fm,lm,msk,nl,fl,mode);

% GET_STR  Returns CellArray: { KeyWord String }
%
% im     MarkerIndex in str
% fm     MarkerFlag
% lm     MarkerLength
% fl     [ KeyWordFlag  CommentFlag ]
% nl     NewLineMarker                   % for mode = 1
% mode   0  KeyWord   and Comment
%        1  Seperator and ContinueMarker
%

get_value = ( imag(mode) == 0 );

mode = real(mode);

%---------------------------------------------------------------------
% Append KeyWordStart

im = cat( 2 , im , size(str,2)+1 );
fm = cat( 2 , fm , fl(1) );

%---------------------------------------------------------------------

ik0 = find( fm == fl(1) );  % KeyWordStart
icm = find( fm == fl(2) );  % Comment

%---------------------------------------------------------------------

nk = size(ik0,2) - 1;   % KeyWordNumber

ini = cell(nk,1+get_value);

ini(:) = { char(zeros(1,0)) };

ok = zeros(nk,1);

for ii = 1 : nk
  
    %----------------------------------------------------------
    % KeyWord
    % !!! FM starts with ZERO !!!

    i0 = ik0(ii) - mode;
    i1 = ik0(ii) - mode + 1;

    i00 = im(i0) + lm(fm(i0)+1);
    i11 = im(i1) - 1;

    ini{ii,1} = mask( rmblank( str(i00:i11) , 2 ) , msk );

    ok(ii) = ~isempty(ini{ii,1});

    if ok(ii) & get_value

        %----------------------------------------------------------
        % Extract Full String for KeyWord

        i00 = im(i1) + lm(fm(i1)+1) * ( fm(i1) == 2 );

        i1  = ik0(ii+1) - mode;
        i11 = im(i1) - 1;

        s = str(i00:i11);

        %----------------------------------------------------------
        % Find CommentMarker for KeyWord

        ic = icm( find( ( i0+mode < icm ) & ( icm < i1 ) ) );

        if mode
           % Add Seperator
           ic = cat( 2 , ik0(ii) , ic );
        end

        %----------------------------------------------------------
        % Remove Comment for KeyWord
        % Extract from Seperator/ContinueMarker to next Marker==NewLine

        if ~isempty(ic)

             i0  = im(ic) - i00 + 1;   % Start
             i1  = im(ic+1) - im(ic);  % Lenght

             if mode
                jj = ( fm(ic) == fl(2) );  % ContinueMarker
                dd = lm(fm(ic)+1) - jj;
                i0 = i0 + dd;
                i1 = i1 - dd;
                s(i0(find(jj))) = nl;      % Continue --> NewLine
             else
                i1 = i1 + ( im(ic-1) == im(ic)-1 );
             end

             ic  = grp2ind( i0 , i1 );

             if mode
                ic(find(ic>size(s,2))) = [];
                s = s(ic);
             else
                s(ic) = [];
             end

        end

        %----------------------------------------------------------
        % String for KeyWord

        ini{ii,2} = mask( rmblank(s,2) , msk );

    end
 
end

ini = ini(find(ok),:);

ik0 = ik0( 1 : nk );

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = mask(str,msk,quota);

% MASK   Mask Characters in String

rev = ( nargin < 3 );   % Reverse

for ii = 1 : size(msk,1)

    if ~isempty(msk{ii,1}) & ~isempty(msk{ii,2})

        if rev
           str = strrep( str , msk{ii,2} , msk{ii,1} );
        else
           str = strrep( str , cat(2,quota,msk{ii,1},quota) , msk{ii,2} );
        end

    end

end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ii = grp2ind(i0,l,s);

% GRP2IND  Built IndexVector from StartIndex and Length
%
% Index = GRP2IND( StartIndex , GroupLength )
%
% Index = GRP2IND( StartIndex , GroupLength , LowSampleStep  )
%

if isempty(i0);
   ii = [];
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

n = size(l,1);

l = ceil( l / s );

ii = s * ones(sum(l),1);
jj = cumsum( cat(1,1,l) , 1 );

ii(jj(1:n)) = i0;

if n > 1
   ii(jj(2:n)) = ii(jj(2:n))-(i0(1:n-1,1)+s*(l(1:n-1)-1));
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

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

