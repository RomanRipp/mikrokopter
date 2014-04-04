function [msg,head,dat] = read_asc(file,sp,cl,cm);

% READ_ASC   Reads Ascii-DataFile with Header
%
% [ Msg , Header , Data ] = READ_ASC( File )
% [ Msg , Header , Data ] = READ_ASC( String )
% [ Msg , Header , Data ] = READ_ASC( CellString )
%
% Header = { Parameter  String }
% Data   = DataMatrice
%
% In the Header the Parameter-Value-Statement is separated 
%   by the Seperator "="
%
% Comments are marked with the CommentMarker "%"
%
% To redefine this Markers use:
%
% READ_ASC( File , Separator , ContinueMarker , CommentMarker ) 
%
%
% see also: READ_PAR, READ_SBE
%

Nin  = nargin;
Nout = nargout;

head = cell(0,2);
dat  = [];

if Nin < 1
   msg = 'Input File is missing.';
   return
elseif ~( iscellstr(file) | ( ischar(file) & ~isempty(file) & ...
         ( prod(size(file)) == size(file,2) ) ) )
   msg = 'File must be a nonempty String or CellArray of Strings.';
   return
end


if Nin < 2
   sp = '=';
end
if Nin < 3
   cl = '';
end
if Nin < 4
   cm = '%';
end

%**************************************************
% Special Characters

bl  = char(32);      % Blank
nl  = char(10);      % NewLine

nrs = '0123456789';  % NumberCharacters
opt = '.e+-';        % Optional Characters

tab = char(9);       % TAB to replace with Blank

acc = 6;             % Accuracy for Dummy

frm = sprintf( '%s%.0f.%.0ff' , '%' , 2*acc+1 , acc );

dummy_val = 10^(acc) - 10^(-acc);
dummy_str = sprintf( frm , dummy_val );


%**************************************************
% Read Header and Data

spc = char(32*ones(1,6));

isf = chkstr(file,1);
if isf
   isf = all( file > 32 );
end
 
if isf
   fprintf(1,'Read  %s  ...  ',file);
else
   fprintf(1,'Analyse String  ...  ');
end

try
   [msg,head,str] = readpar(file,sp,cl,cm);
catch 
    msg = lasterr;
end

if isempty(msg)
   fprintf(1,'ok\n');
else
   fprintf(1,'error\n');
   msg = sprintf('Error using READPAR.\n%s',msg);
end

if ~isempty(msg) | ( Nout < 3 ) | isempty(str)
    return
end

%**************************************************
% Get Number of Columns

jj = strcmp(lower(head(:,1)),'columns');

if any(jj)
   nc = sum( head{jj,2} == ':' ) + 1;
else
   nc = NaN;
end

%**************************************************
% Check DataString

[dat,msg,nn,nr] = scandat(str,nc,spc);

if ~isempty(msg)

    str = chkdat(str,nl,bl,tab,nrs,opt,dummy_str,spc);

    if isempty(str)
       return
    end

    [dat,msg,nn,nr] = scandat(str,nc,spc);

end

if isempty(dat)
   return
end

if isnan(nr) | isnan(nc)

   [str,nc] = checkcol(str,nc,nl,bl,spc);

   if ~( size(dat,2) == nc )
       [dat,msg,nn,nr] = scandat(str,nc,spc);
   end

   if isempty(dat)
      return
   end

end

dat(find(dat==dummy_val)) = NaN;

%**************************************************
% Check Columns

if ~( size(dat,2) == nc )
    fprintf( 1 , '%sWarning: %s (%.0f) doesn''t match expected %s (%.0f).\n', ...
             spc , 'Number of Datasets' , nn , 'Number of Columns' , nc );
end


%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = chkdat(str,nl,bl,tab,nrs,opt,dummy_str,spc)

if all( ( str == nl ) | ( str == bl ) | ( str == tab ) )
   str = '';
   return
end

%--------------------------------------------------
% TAB --> Blank

str = strrep(str,tab,bl);

%--------------------------------------------------
% Add NewLine at Begin and End

str = lower( cat( 2 , nl , str , nl ) );

%--------------------------------------------------
% Replace NaN's, surrounded by Newline or Blank

nn = 'nan';

n  = size(str,2);

for cc = { [nl nl]  [nl bl]  [bl bl]  [bl nl] }
     c = cc{1};
   str = strrep( str , [ c(1)  nn         c(2) ] , ...
                       [ c(1)  dummy_str  c(2) ] );
end

if ~( n == size(str,2) )

    [dat,wrn] = scandat(str,NaN,spc);

    if isempty(wrn)
       return
    end

end

%--------------------------------------------------
% Check for MonthStatements

str = mon2num(str);

%--------------------------------------------------
% Check for valid Characters: nrs | opt | nl | bl

ok  = ( ( min(nrs) <= str ) & ( str <= max(nrs) ) );

for cc = [ opt nl bl ]
    ok = ( ok | ( str == cc ) );
end

if ~all(ok)
    fprintf(1,'%sWarning: Invalid Characters.\n',spc);
    str(find(~ok)) = bl;
end

%--------------------------------------------------
% Check for Duplicate Characters

for cc = opt

    % No Duplicate Characters
    nn = size(str,2);
    str = remfollw(str,cc);
    mm = size(str,2);
    if ~( nn == mm )
        fprintf(1,'%sWarning: Duplicate "%s"-Characters.\n',spc,cc);
    end

end

%--------------------------------------------------
% Check for Single Characters

for cc = opt

    ii = ( str == cc );

    if any(ii)

       ii = find(ii);

       jj = ( ( ( str(ii-1) == nl ) | ( str(ii-1) == bl ) ) & ...
              ( ( str(ii+1) == nl ) | ( str(ii+1) == bl ) )       );

       if any(jj)
          fprintf(1,'%sWarning: Single "%s"-Characters.\n',spc,cc);
              jj  = ii(find(jj));
          str(jj) = bl;
       end

    end

end

%--------------------------------------------------
% Check for Valid "e"-Statement

cc = 'e';

ii = ( str == cc );

if any(ii)

   ii = find(ii);

   % No "+-" bl nl before "e"

   jj = ( 1 * ( str(ii-1) == nl  ) + 1 * ( str(ii-1) == bl  ) + ...
          2 * ( str(ii-1) == '+' ) + 2 * ( str(ii-1) == '-' ) );

   if any(jj)
      fprintf(1,'%sWarning: Invalid Characters before "%s"-Character.\n',spc,cc);
      str(ii(find(jj))) = bl;       % Remove Character
      str(ii(find(jj==2))-1) = bl;  % Remove "+-"-Character
      ii(find(jj)) = [];
   end

   if ~isempty(ii)

       % No "."  bl nl after  "e"

       jj = ( 1 * ( str(ii+1) == nl  ) + 1 * ( str(ii+1) == bl  ) + ...
              2 * ( str(ii-1) == '.' ) );

       if any(jj)
          fprintf(1,'%sWarning: Invalid Characters after "%s"-Character.\n',spc,cc);
          str(ii(find(jj))) = bl;      % Remove Character
          str(ii(find(jj==2))+1) = bl; % Remove '.'-Character
       end

   end

end

%--------------------------------------------------
% Check for Valid "." - Statement

cc = '.';

ii = ( str == cc );

if any(ii)

   ii = find(ii);

   % No "+-" arround "."
 
   jj = ( ( str(ii-1) == '+' ) | ( str(ii-1) == '-' ) ) & ...
          ( str(ii+1) == '+' ) | ( str(ii+1) == '-' );

   if any(jj)
      fprintf(1,'%sWarning: Invalid Characters arround "%s"-Character.\n',spc,cc);
          jj    = ii(find(jj));
      str(jj+0) = bl;           % Remove Character
      str(jj-1) = bl;           % Remove Character before
   end

end

%--------------------------------------------------
% Check for Valid "+-" - Statement

cc = '+-';

ii = ( ( str == cc(1) ) | ( str == cc(2) ) );

if any(ii)

   ii = find(ii);

   % No "+-" before "+-"

   jj = ( ( str(ii-1) == '+' ) | ( str(ii-1) == '-' ) );
   if any(jj)
      fprintf(1,'%sWarning: Invalid Characters before "%s"-Character.\n',spc,cc);
             jj     = find(jj);
      str(ii(jj)-1) = bl;       % Remove Characters before
          ii(jj)    = [];
   end

   if ~isempty(ii)

       % No Blank after Characters

       jj = ( ( str(ii+1) == bl ) | ( str(ii+1) == nl ) );

       if any(jj)
          fprintf(1,'%sWarning: Blank after "%s"-Character.\n',spc,cc);
              jj  = ii(find(jj));
          str(jj) = bl;            % Remove Characters
       end

   end

end

%--------------------------------------------------
% Check for Single Characters

for cc = opt

    ii = ( str == cc );

    if any(ii)

       ii = find(ii);

       jj = ( ( ( str(ii-1) == nl ) | ( str(ii-1) == bl ) ) & ...
              ( ( str(ii+1) == nl ) | ( str(ii+1) == bl ) )       );

       if any(jj)
          fprintf(1,'%sWarning: Single "%s"-Characters.\n',spc,cc);
              jj  = ii(find(jj));
          str(jj) = bl;
       end

    end

end

%--------------------------------------------------

if all( ( str == nl ) | ( str == bl ) )
   str = '';
end

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [str,nc] = checkcol(str,nc,nl,bl,spc);

%--------------------------------------------------
% Remove following NewLine and Blanks

for cc = [ bl nl ]
    str = remfollw(str,cc);
end

%--------------------------------------------------
% Remove Blanks, surrounding NewLine

str = remfollw(str,[bl nl],nl);

%--------------------------------------------------
% Check Number of Blanks between NewLines

ncl = cumsum( str == bl );
inl =   find( str == nl );

ncl = diff(ncl(inl),1,2) + 1;

if all( ncl == ncl(1) )

   nc = ncl(1);

else

   fprintf(1,'%sWarning: Different Number of Columns.\n',spc);

   if isnan(nc) | all( ncl < nc )
      nc = median(ncl);
   end

   jj = ~( ncl == nc );

   if any(jj)

      jj = find(jj);

      nn = size(inl,2) - 1;

      ii = inl(1:nn) + 1;    % StartIndex
      ll = diff(inl,1,2);    % Length
       
      ii = ii(jj);
      ll = ll(jj);

      mm = size(ll,2);

      fprintf('%sWarning: Remove %.0f of %.0f Lines.\n',spc,mm,nn);

      jj = ones(1,sum(ll));
      kk = cumsum( cat(2,1,ll) , 2 );

      jj(kk(1:mm)) = ii;

      if mm > 1
         jj(kk(2:mm)) = jj(kk(2:mm))-(ii(1:mm-1)+(ll(1:mm-1)-1));
      end

      jj = cumsum(jj);
 
      str(jj) = [];

   end
    
end

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [dat,wrn,nn,nr] = scandat(str,nc,spc)

nr = NaN;

[dat,nn,wrn] = sscanf(str,'%f');

if ~isempty(wrn)
    fprintf(1,'%sWarning: %s\n',spc,wrn);
end

if isempty(dat)
   dat = [];
   nr  = 0;
   return
end

if isnan(nc) | ( nc == 1 )
   nr = nn;
elseif ( mod(nn,nc) == 0 )
   nr = nn / nc;
   dat = reshape(dat,nc,nr);
   dat = permute(dat,[2 1]);
end


%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function x = remfollw(x,z,r)

% REMFOLLW  Remove Following equal Elements from Matrice
%
% X = REMFOLLW( X , Z , R )
%
%  Removes form X the follwing Values, matching the Elements of Z.
%
%  The single Values will replaced by R (optional)
%

Nin = nargin;

%------------------------------------------------------

if Nin < 1
   x = [];
elseif ~( isnumeric(x) | ischar(x) )
   error('X must be a Numeric or Character-Array.');
end

%------------------------------------------------------

if Nin < 2
   z = [];
elseif ~( isnumeric(z) | ischar(z) )
   error('Z must be a Numeric or Character-Array.');
end

%------------------------------------------------------

if Nin < 3
   r = [];
elseif ~( isnumeric(r) | ischar(r) )
   error('R must be a Numeric or Character-Array.');
elseif ~( isempty(r)  |  ( prod(size(r)) == 1 ) )
   error('R must be a single Element or EMPTY.');
end

%------------------------------------------------------

if isempty(x) | isempty(z)
   return
end

%*******************************************************

m = prod(size(z));

jj = zeros(size(x));

for ii = 1 : m
    
    jj = ( jj | ( x == z(ii) ) );

end

if ~any(jj)
    return
end

%------------------------------------------------------
% IND2GRP

jj = jj(:);

jj = find(jj);

n  = size(jj,1);

ii = cat( 1 , 1 , find( diff(jj,1,1) > 1 )+1 , n+1 );
ll = diff(ii,1,1);
ii = jj(ii(1:end-1));

jj = ( ll > 1 );

if ~any(jj)
    return
end

jj = find(jj);

ii = ii(jj);   % StartIndex
ll = ll(jj);   % Length

%------------------------------------------------------
% GRP2IND

ii = ii + 1;
ll = ll - 1;

m = size(ll,1);

jj = ones(sum(ll),1);
kk = cumsum( cat(1,1,ll) , 1 );

jj(kk(1:m)) = ii;

if m > 1
   jj(kk(2:m)) = jj(kk(2:m))-(ii(1:m-1,1)+(ll(1:m-1)-1));
end

jj = cumsum(jj,1);

%------------------------------------------------------
% Replace and Remove

if ~isempty(r)
    x(ii-1) = r;
end

x(jj) = [];


%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,ini,str] = readpar(file,sp,cl,cm);

% READ_PAR  Read ParameterFile
%
% [Msg,V] = READ_PAR( ParameterFile ,  Seperator , ContinueMarker , CommentMarker )
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
% The Comment's are marked in the File with the CommentMarker.
%   A Comment is valid until End of Line (NL).
%
%
% Use "" around a Marker to mask them.
%
% Defaults:
%
%    Seperator      =   { '='  '' }   (no special ParameterMarker)
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

if chkcstr(file,1)
   str = sprintf('%s\n',file{:});
elseif ~( ischar(file) & ( ndims(file) == 2 ) )
   msg = cat( 1 , msg , {'Input File must be a CharacterArray or CellArray of Strings.'} );
else
  if isempty(file)
     str = '';
  elseif ( size(file,1) > 1 ) | any( file == 10 ) | any( file == 13 ) 
     str = cat( 2 , file , char(10*ones(size(file,1),1)) );
     str = permute(str,[2 1]);
     str = str(:);
     str = permute(str,[2 1]);
  else
     is_file = 1;
  end
end
   
%---------------------------------------------------------------------
% Seperator

[ok,sp] = chkcstr(sp);
if ok
   ok = ~isempty(sp{1});
   if ok
      sp = sp(:);
      sp = cat( 1 , sp , {''} );
      sp = sp([1 2]);
      st = rmblank(sp{2},2);    % ParameterStart
      sp = rmblank(sp{1},2);    % ParameterSeperator
      ok = ~isempty(sp);
   end
end

if ~ok
   msg = cat( 1 , msg , ...
              {'Input Seperator must be a CharacterArray or CellArray of Strings.'} , ...
              {' 2 Seperator defines { ParameterSeperator ParameterStart }.' } , ...
              {' ParameterSeperator must be nonempty.' }  );
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

   [msg,str] = loadfile(file,'char',ext);

   if ~isempty(msg)
       try
          msg = catmsg( 'Invalid ParameterFile.' , msg , i );
       catch
          msg = sprintf('%s\n%s','Invalid ParameterFile.',msg);
       end
       return
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

msk = msk(:,[2 1]);  % Flip to remask

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
   return
end

%*********************************************************************
% Remove Property, Value, Comment

nn = size(im,2);

ind = ( 1 : (nn-1) );

ok      = ~( fm == nm );
ok(ind) = ( ok(ind) | ( ok(ind+1) & ~( fm(ind+1) == (nm-1) ) ) );

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


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ini,im,fm,lm,nm] = get_par(str,st,sp,cl,cm,nl,msk)

ini = cell(0,2);  % { Parameter Value }

%---------------------------------------------------------------------
% Get MarkerIndize, Flag and Length

mark = { st sp cl cm nl };

nm = size(mark,2) - 1;         % !!! Start with ZERO !!!

[im,fm,lm] = get_mark(str,mark,2);

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
   ok = any( fm == 1 );
end
if ~ok
   return
end

%---------------------------------------------------------------------
% Get Parameter-Value

[ini,i0,ok] = get_str(str,im,fm,lm,[1 2],1,msk);

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

if ~( lm(1) == 0 )
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

    jj = findstr( str , mark{ii} );

    if ~isempty(jj)

       im = cat( 2 , im ,      jj );
       fm = cat( 2 , fm , ii+0*jj );
    
    elseif any( ii == ret )

       return
 
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

% chk = [ Flag FlagBefore True/False ];
%                         -True == after

cc = ones(size(im));

for ii = 1 : size(chk,1)

    jj = find( fm == chk(ii,1) );

    if ~isempty(jj)

        ck = ~( chk(ii,3) == 0 );         % True or False

        cj = 1 - 2 * ( chk(ii,3) >= 0 );  % Before or After

        ok = ( fm(jj+cj) == chk(ii,2) );

        ok = ( 1 - ck ) + ( 2*ck - 1 ) * ok;

        jj = jj(find(~ok));

        cc(jj) = 0;

        if mode

          fm(jj) = [];
          im(jj) = [];

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

function [ini,ik0,ok] = get_str(str,im,fm,lm,fl,mode,msk);

% GET_STR  Returns CellArray: { KeyWord String }
%
% im     MarkerIndex in str
% fm     MarkerFlag
% lm     MarkerLength
% fl     [ KeyWordFlag  CommentFlag ]
% mode   0  KeyWord   and Comment
%        1  Seperator and ContinueMarker
%

%---------------------------------------------------------------------
% Append KeyWordStart

im = cat( 2 , im , size(str,2)+1 );
fm = cat( 2 , fm , fl(1) );

%---------------------------------------------------------------------

ik0 = find( fm == fl(1) );  % KeyWordStart
icm = find( fm == fl(2) );  % Comment

%---------------------------------------------------------------------

nk = size(ik0,2) - 1;   % KeyWordNumber

ini = cell(nk,2);

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

    if ok(ii)

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

             i0 = i0 + lm(fm(ic)+1) * mode;
             i1 = i1 + ( lm(fm(ic+1)+1) - lm(fm(ic)+1) ) *   mode + ...
                       (    im(ic-1)    ==   im(ic)-1  ) * (~mode);

             i1(end) = i1(end) - lm(fm(ic(end)+1)+1)*mode;

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

if nargin < 3
   quota = '';
end

for ii = 1 : size(msk,1)

    if ~isempty(msk{ii,1}) & ~isempty(msk{ii,2})

       str = strrep( str , cat(2,quota,msk{ii,1},quota) , msk{ii,2} );

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

m   = 60*2^20; % 2 MBytes;  % MaximumFileSize [Bytes]
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

    if ~all( ok )

       bb(find(~ok)) = 32;

       warning('Invalid Characters in File.')

       %%% msg = [ msg0 'Invalid Characters in File.' ];

    end

    bb = char(bb(:)');


  end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = mon2num(str,mon,cc);

% MON2NUM  Replace a MonthString in  String by a Number
%
% STR = MON2NUM(STR)
%
% Replace in the String STR all english and german 
%  abbrevations (3 characters) for months by their Number.
%
% STR = MON2NUM( STR , MonthString , ReplaceCharacter )
%
% MonthString define the strings to replace by the Number
% MonthString can be CharacterArray or CellArray wich 
%  elements are Strings or CellArray of Strings for multiple
%  matches per Number. MON2NUM without any Input returns the
%  default value for MonthString:
%
%  MonthString = MON2NUM
%
% ReplaceCharacters are characters surrounding the strings 
%  from MonthString which will replaced by a Blank.
%  default: '-./+*^\:'   mathematical Operators
%
%
% see also: FINDSTR, STRREP
%

Nin = nargin;

%-------------------------------------------------

if Nin > 1
   if isempty(str)
      return
   elseif ~chkstr(str)
      error('String required.');
   end
end

%-------------------------------------------------

if Nin < 2

   mon = {  { 'JAN'       }
            { 'FEB'       }
            { 'MAR' 'MÄR' }
            { 'APR'       }
            { 'MAY' 'MAI' } 
            { 'JUN'       }
            { 'JUL'       }
            { 'AUG'       }
            { 'SEP'       }
            { 'OCT' 'OKT' }
            { 'NOV'       } 
            { 'DEC' 'DEZ' }  };

else

   if ischar(mon)
      mon = cellstr(mon);
   end

   if ~iscell(mon)
       error('MonthString must be a CharacterArray or CellArray.');
   end

   mon = mon(:);

   n  = size(mon,1);

   for ii = 1 : size(mon,1)
       [ok,mon{ii}] = chkcstr(mon{ii});
       if ~ok
           error('Elements of MonthString must be Strings or CellArrays of Strings.');
       end
       mon{ii} = upper(mon{ii}(:)');  % Row
   end

end

%-------------------------------------------------

if Nin < 3
   cc = '-./+*^\:';
elseif ~isempty(cc)
    if ~ischar(cc)
        error('Invalid Input for ReplaceCharacters.')
    end
    cc = cc(:)';
end

%-------------------------------------------------

if Nin < 1
   if nargout == 0
      fprintf(1,'\n');
      for ii = 1 : size(mon,1)
          len = size(mon{ii}{1},2);
          frm = sprintf('%%%.0f.%.0fd',len,len);
          rep = sprintf(frm,ii);
          fprintf(1,'   %s : {',rep);
          fprintf(1,' ''%s'' ',mon{ii}{:});
          fprintf(1,'}\n');
      end
      fprintf(1,'\n');
      clear str
   else
      str = mon;
   end
   return
end

%*************************************************

% UpperCase

upp = upper(str);

siz = size(str,2);

for ii = 1 : size(mon,1)

    for mm = mon{ii}

        len = max(1,size(mm{1},2));
        frm = sprintf('%%%.0f.%.0fd',len,len);
        rep = sprintf(frm,ii);
        ins = ~( size(rep,2) == len );

        if len == siz

           if strcmp( upp , mm{1} )
              str = rep;
              if ins
                 upp = upper(str);
                 siz = size(str,2);
              end
           end
  
        elseif ( 0 < len ) & ( len < siz )

           kk = findstr(upp,mm{1});
           if ~isempty(kk)
               for jj = kk(end:-1:1)
                   ll =  jj + [ -1  len ];
                   ll =  ll( (1+(ll(1)<1)) : (2-(ll(2)>siz)) );
                   if ~any(isletter(str(ll)));
                       if ~isempty(cc)
                           for c = cc
                               ok = ( str(ll) == c );
                               if any(ok)
                                   ok = find(ok);
                                   str(ll(ok)) = ' ';
                               end
                           end
                       end
                       if ins
                          str = cat( 2 , str(1:ll(1)) , rep , str(ll(2):end) );
                       else
                          str(jj-1+(1:len)) = rep;
                       end
                   end
               end
               if ins
                  upp = upper(str);
                  siz = size(str,2);
               end
           end

        end

    end

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

