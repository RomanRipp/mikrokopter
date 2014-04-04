function [msg,head,data,ini] = read_sbe(varargin)

% READ_SBE  Reads SBE-Ascii-DataFiles with decimal Values (CNV)
%
% [ Msg , Header , Data , HeadINI ] = READ_SBE( FileName )
%
% READ_SBE( FileName , MaxNumber , -LogID , CommentMarker )
%
% MaxNumber:  Maximum Number of HeaderLines
%
% Use a NonZero imaginary Part to return the Header only:
%
% [ Msg , Header , HeadINI ] = READ_SBE( FileName , MaxNumber+i , ... )
%
% In this case MaxNumber can be also ZERO to read the Header complete.
%
% Note: The Time in JulianDays starts with Day 1 at 01.Jan.YYYY 00:00
% 
% see also: READ_PAR, READ_HEX
%

msg  = '';
head = [];
data = [];
ini  = {};

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

% SBE-HeaderInitialisation: { Field  Marker  Seperator }
    
ini = { 'head'   { '**'        ':' }
        'cols'   { '# name'    '=' }
        'sens'   { '# sensor'  '=' }
        'lim'    { '# span'    '=' }
        'par'    { '#'         '=' }
        'proc'   { '*'         '=' }
        'sbe'    { '* SBE'     ''  }
        'info'   { '*'         ':' }
        'misc'   { '*'         ''  }
        'run'    { ''          '=' }  };

% HeaderFields, to search for Columns
%
%  MultipleCharacter Marker: ColumnNumber == ParameterNumber
%    SingleCharacter Marker: Search for ColumnParameter and ColumnSeperator
%

cf  = { 'cols'  'par'  'proc'  'run' };  

cp  = 'columns';  % ColumnParameter
cs  = ':';        % ColumnSeperator

nh  = 512;       % Maximum Number of HeaderLines to read

%****************************************************************
% Special Characters

cl  = '';    % ContinueMarker
cm  = '%';   % CommentMarker

bl  = char(32);      % Blank
nl  = char(10);      % NewLine

nrs = '0123456789';  % NumberCharacters
opt = '.E+-';        % Optional Characters, UPPERCASE !!!
sep = ',;:';         % Valid Seperator

tab = char([9 160]); % TAB to replace with Blank

acc = 6;             % Accuracy for Dummy

frm = sprintf( '%s%.0f.%.0ff' , '%' , 2*acc+1 , acc );

dummy_val = 10^(acc) - 10^(-acc);
dummy_str = sprintf( frm , dummy_val );

spc = char(32*ones(1,6));   % Space for WarningMessage

%****************************************************************

% Check Inputs

[msg,file,siz,lid,nh,nd,cm] = checkin(varargin,nh,cm);

if ~isempty(msg)
    return
end

%****************************************************************
% Open File

fprintf(lid,'\nOpen SBE-File: %s  ... ',file);

fid = fopen(file,'r');
if fid == -1
   fprintf(lid,'error\n');
   msg = sprintf('Can''t open File: %s',file);
   return
end
  
fprintf(lid,'ok\n');

%****************************************************************
% Read Header

fprintf(lid,'Read Header  ... ');

[str,dat,ok,nhl] = readhead(fid,nh);

pos = ftell(fid);

if any( ok == -1 )
   fclose(fid);
end

if ok == 0
   fprintf(lid,'Warning: Maximum Number %.0f of HeaderLines reached',nh);
end

if isempty(str)
   fprintf(lid,'; Warning: No Header');
end

fprintf(lid,'\n');

%****************************************************************
% Analyze Header

nc = NaN;   % Number of Columns from Header

if isempty(str)
   
   head = permute(ini,[2 1]);
   head(2,:) = {{{}}};

   head = struct(head{:});

else

   fprintf(lid,'Analyze Header  ... ');

   [wrn,head,ini] = gethead(str,ini,cl,cm);

   if ~isempty(wrn)
       fprintf(lid,'error');
   else
       fprintf(lid,'ok');
   end

   %-----------------------------------
   % Get Number of Columns from Header

   if ~isempty(head) & isa(head,'struct')
       for ff = cf
           hv = getfield(head,ff{1});  % Value
           if ~isempty(hv)
               hm = getfield(ini,ff{1});  % { Marker Seperator }
               if size(hm{1},2) > 1
                  nc = size(hv,1);
                  break
               else
                  jj = strcmp(lower(hv(:,1)),cp);
                  if any(jj)
                     jj = max(find(jj));
                     nc = sum( hv{jj,2} == cs ) + 1;
                     break
                  end
               end
           end
       end
   end

   if ~isnan(nc)
       fprintf(lid,';  %.0f Columns',nc);
   end

   %-----------------------------------

   fprintf(lid,'\n');

end

head.lines = nhl;
head.bytes = pos;

%****************************************************************

if nd | ( nargout < 3 )

   data = ini;

   fclose(fid);

   return

end

%****************************************************************
% Read Data

if ok == -1

   fprintf(lid,'Warning: No Data\n');

else

   %--------------------------------------
   % Get Number of Columns from FirstLine

   nd = NaN;

   if ~isempty(dat)

       [str,dat,wrn] = chkdat(lid,dat,cm,nl,bl,tab,sep,nrs,opt,dummy_str,spc);

       if isempty(wrn)
          nd = size(dat,1);
       else
          [str,nd] = chkcol(lid,str,NaN,nl,bl,spc);
       end

       if isnan(nc)
          fprintf(lid,'Number of Columns in 1. DataLine:  %.0f\n',nd);
       elseif ~( nc == nd )
          fprintf(lid,'Warning: Number of Columns in 1. DataLine:  %.0f\n',nd);
       end

   end

   %-----------------------------------
   % Scan Data, Build Matrice

   if ~( isnan(nc) & isnan(nd) )

       fprintf(lid,'Scan DataFile  ... ');

       [ok,data] = scanfile(fid,siz,sep);

       if ~ok
           fprintf(lid,'error\n');
       else
           fprintf(lid,'ok');
           if isempty(data)
              fprintf(lid,';  Warning: No Data\n');
           else
              fprintf(lid,'\n');
              [ok,data] = data2mat(lid,data,nc,nd);
           end
       end

   end

   %-----------------------------------
   % Read DataString
 
   if ~ok
       fprintf(lid,'Read DataString  ... ');
       data = fread(fid,'char');
       if isequal(data,-1)
          fprintf(lid,'error\n');
          data = [];
       else
          fprintf(lid,'ok\n');
          data = char(permute(data,[2 1]));
       end
   end

   %-----------------------------------

   fclose(fid);

   %-----------------------------------
   % CheckDataString

   if ~ok & ~isempty(data)

       fprintf(lid,'Check DataString  ...\n');

       [str,data,wrn] = chkdat(lid,data,cm,nl,bl,tab,sep,nrs,opt,dummy_str,spc);

       ok = ( isempty(wrn) & ~( isnan(nc) & isnan(nd) ) );
       if ok
          [ok,data] = data2mat(lid,data,nc,nd);
       end

       if ~ok
           fprintf(lid,'Check Columns  ...\n');
           [str,nc] = chkcol(lid,str,nc,nl,bl,spc);
           [ok,data] = data2mat( lid , sscanf(str,'%f') , nc );
       end

   end

end

fprintf(lid,'\n');

%****************************************************************

if ~isempty(data)

    data( find( data == dummy_val ) ) = NaN;

end

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,file,siz,lid,nh,nd,cm] = checkin(vin,nh,cm);

msg  = '';
file = '';
siz  = 0;
lid  = 1;
nd   = 0;  % No Data

nv = prod(size(vin));

if nv == 0
   msg = 'Input FileName is missing.';
   return
end

msg = cell(0,1);

%************************************************************
% Check File

file = vin{1};

if ~( ischar(file) & ~isempty(file) & ...
      ( prod(size(file)) == size(file,2) ) )

    msg = cat( 1 , msg , {'FileName must be a nonempty String.'} );

elseif ~( exist(file,'file') == 2 )

    msg = cat( 1 , msg , {sprintf('File doesn''t exist: %s',file)} );

else

    d = dir(file);
    if isempty(d)
       f = which(file);
       if ~isempty(f)
           file = f;
           d    = dir(file);
       end
    end

    if isempty(d)
       msg = cat( 1 , msg , {sprintf('Can''t get DIR(''%s'').',file)} );
    else
        siz = getfield(d,'bytes');
        if siz == 0
           msg = cat( 1 , msg , {sprintf('Empty File: %s',file)} );
        end
    end

end

%************************************************************
% Check for LogID, MaxNumber, CommentMarker

for ii = 2 : nv

    v = vin{ii};

    if ischar(v)
       if isempty(v)
          cm = '';
       elseif ( prod(size(v)) == size(v,2) )
          cm = v;
       else
          msg = cat(1,msg,{'CharacterInput must be a String for CommentMarker.'});
       end
    elseif isnumeric(v) & ( prod(size(v)) <= 1 )
       if isempty(v)
          lid = 0;
       elseif ( real(v) <= 0 ) & ( imag(v) == 0 )
          lid = abs(v);
          if isempty(fopen(lid))
             msg = cat(1,msg,{'Invalid LogID.'});
             lid = 0;
          end
       elseif ( mod(real(v),1) == 0 )
          nd = ~( imag(v) == 0 );
           v = real(v);
          if ~( ( v == 0 ) & nd )
              nh = v;
          end
       else
          msg = cat(1,msg,{'MaxNumber must be a Integer.'});
       end
    else
       msg = cat(1,msg,{'Following Inputs must be single Numerics or Strings.'});
    end

end

%************************************************************

if isempty(msg)
   msg = '';
else
   msg = sprintf('%s\n',msg{:});
end

 
%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [str,dat,ok,z] = readhead(fid,mm)

str = '';
dat = '';   % First DataLine

val = '0123456789';  % Valid NumberCharacters for Start of DataLine
opt = '.';           % Valid optional Character, "val" has to follow
sgn = '+-';          % Valid SignCharacter, "opt" or "val" has to follow

nl = char(10);

z   = 0;
ok  = 0;

while ~ok & ( z < mm )

       z = z + 1;

       pos = ftell(fid);

       bb = fgetl(fid);
       if isequal(bb,-1)
          ok = -1;
          break
       end

       bb = rmblank(bb,2);

       if ~isempty(bb)

           ok = any( bb(1) == val ) ;

           if ~ok
               si = size(bb,2);
               ok = ( any( bb(1) == opt ) & ( si > 1 ) );
               if ok
                  ok = any( bb(2) == val );
               else
                  ok = ( any( bb(1) == sgn ) & ( si > 1 ) );
                  if ok
                     ok = any( bb(2) == val );
                     if ~ok
                         ok = ( any( bb(2) == opt ) & ( si > 2 ) );
                         if ok
                            ok = any( bb(3) == val );
                         end
                     end
                  end
               end
           end  

           if ok
              dat = bb;
              act = ftell(fid);
              fseek(fid,pos-act,'cof');
              break
           end

           str = cat(2,str,nl,bb);

       end

end

z = z - ok;

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,dat,ini] = gethead(str,ini,cl,cm)

msg = '';

%****************************************************************
% Initialisation: { Field  Marker  Seperator }
    
n = size(ini,1);

dat = permute(ini(:,[1 1]),[2 1]);

dat(2,:) = {{cell(0,2)}};   % { Parameter  Value }

dat = struct(dat{:});

msg    = cell(n,1);
msg(:) = {''};

ok = zeros(n,1);

for ii = 1 : n

    sep = ini{ii,2};

    ini{ii,2} = {sep};

    [msg{ii},d,str] = read_par(str,sep([2 1]),cl,cm);

    ok(ii) = isempty(msg{ii});

    if ok(ii)
       dat = setfield(dat,ini{ii,1},d);
    end

end

ini = permute(ini,[2 1]);
ini = struct(ini{:});

if all(ok)
   msg = '';
   return
end

ok = find(~ok);

msg = sprintf('%s\n',msg{ok}); 


%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,data] = scanfile(fid,siz,sep);


sep = cat( 1 , {''} , cellstr(sep(:)) );
sep = permute(sep,[2 1]);

pos = ftell(fid);

%****************************************************************
% FSCANF

frm = '%f';

for sp = sep

    data = fscanf( fid , cat(2,frm,sp{1}) );

     act = ftell(fid);

    fseek( fid , pos-act , 'cof' );

    ok = ( act == siz );

    if ok
       break
    end

end

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,data] = data2mat(lid,data,nc,nd)


fprintf(lid,'Reshape Data  ... ');

%****************************************************************
% Check numeric Data

nn = size(data,1);

ok = ( nargin == 3 );

if ok

   mm = nc * floor(nn/nc);

   if mm == 0
      data = zeros(0,nc);
   else
      data = data( 1 : mm );
   end

   if ~( mm == nn )
       fprintf(lid,'Warning: Cut  %.0f of %.0f  Values  ... ',nn-mm,nn);
   end

   if mm == 0
      fprintf(lid,'ok;  %.0f Columns\n',size(data,2));
      return
   end

   nn = mm;
   mm = nc;

else

   for mm = [ nc  nd ]
       ok = ( mod(nn,mm) == 0 );
       if ok
          break
       end
   end

end

if ok

   data = reshape(data,mm,nn/mm);
   data = permute(data,[2 1]);

   fprintf(lid,'ok;  %.0f Columns\n',size(data,2));

else

   fprintf(lid,'error\n');

end

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [str,dat,wrn] = chkdat(lid,str,cm,nl,bl,tab,sep,nrs,opt,dummy,spc)

dat = [];
wrn = '';

%**************************************************

str = strrep(str,char([13 10]),nl);  % DOS: CRLF --> LF  
str = strrep(str,char(13),nl);       % MAC: CR   --> LF

%--------------------------------------------------
% Sep, TAB --> Blank

for cc = [ tab sep ]
    str( find( str == cc ) ) = bl;
end

%--------------------------------------------------

if all( ( str == nl ) | ( str == bl )  )
   str = '';
   return
end

%--------------------------------------------------
% Add NewLine at Begin and End, make UPPERCASE !!!

str = upper( cat( 2 , nl , str , nl ) );

%--------------------------------------------------
% Replace NaN's, surrounded by Newline or Blank

nn = 'NAN';

n  = size(str,2);

for cc = { [nl nl]  [nl bl]  [bl bl]  [bl nl] }
     c = cc{1};
   str = strrep( str , [ c(1)  nn     c(2) ] , ...
                       [ c(1)  dummy  c(2) ] );
end

if ~( n == size(str,2) )

    [dat,nn,wrn] = sscanf(str,'%f');

    if isempty(wrn)
       return
    end

end

%--------------------------------------------------
% Check for Month

str = mon2num(str);

%--------------------------------------------------
% Check for Comments

[msg,comm,str] = read_par(str,{'' cm},'','');

%--------------------------------------------------
% Check for valid Characters: nrs | opt | nl | bl

ok  = ( ( min(nrs) <= str ) & ( str <= max(nrs) ) );

for cc = [ opt nl bl ]
    ok = ( ok | ( str == cc ) );
end

if ~all(ok)
    fprintf(lid,'%sWarning: Invalid Characters.\n',spc);
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
        fprintf(lid,'%sWarning: Duplicate "%s"-Characters.\n',spc,cc);
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
          fprintf(lid,'%sWarning: Single "%s"-Characters.\n',spc,cc);
              jj  = ii(find(jj));
          str(jj) = bl;
       end

    end

end

%--------------------------------------------------
% Check for Valid "E"-Statement

cc = 'E';

ii = ( str == cc );

if any(ii)

   ii = find(ii);

   % No "+-" bl nl before "E"

   jj = ( 1 * ( str(ii-1) == nl  ) + 1 * ( str(ii-1) == bl  ) + ...
          2 * ( str(ii-1) == '+' ) + 2 * ( str(ii-1) == '-' ) );

   if any(jj)
      fprintf(lid,'%sWarning: Invalid Characters before "%s"-Character.\n',spc,cc);
      str(ii(find(jj))) = bl;       % Remove Character
      str(ii(find(jj==2))-1) = bl;  % Remove "+-"-Character
      ii(find(jj)) = [];
   end

   if ~isempty(ii)

       % No "."  bl nl after  "E"

       jj = ( 1 * ( str(ii+1) == nl  ) + 1 * ( str(ii+1) == bl  ) + ...
              2 * ( str(ii-1) == '.' ) );

       if any(jj)
          fprintf(lid,'%sWarning: Invalid Characters after "%s"-Character.\n',spc,cc);
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
      fprintf(lid,'%sWarning: Invalid Characters arround "%s"-Character.\n',spc,cc);
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
      fprintf(lid,'%sWarning: Invalid Characters before "%s"-Character.\n',spc,cc);
             jj     = find(jj);
      str(ii(jj)-1) = bl;       % Remove Characters before
          ii(jj)    = [];
   end

   if ~isempty(ii)

       % No Blank after Characters

       jj = ( ( str(ii+1) == bl ) | ( str(ii+1) == nl ) );

       if any(jj)
          fprintf(lid,'%sWarning: Blank after "%s"-Character.\n',spc,cc);
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
          fprintf(lid,'%sWarning: Single "%s"-Characters.\n',spc,cc);
              jj  = ii(find(jj));
          str(jj) = bl;
       end

    end

end

%--------------------------------------------------

if all( ( str == nl ) | ( str == bl ) )
   str = '';
else
   [dat,nn,wrn] = sscanf(str,'%f');
end

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [str,nc] = chkcol(lid,str,nc,nl,bl,spc);

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

   fprintf(lid,'%sWarning: Different Number of Columns.\n',spc);

   if isnan(nc) | ~any( ncl == nc )
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

      fprintf(lid,'%sWarning: Remove %.0f of %.0f Lines.\n',spc,mm,nn);

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

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = mon2num(str);

mon = {    'Jan'
           'Feb'
         { 'Mar' 'Mär' }
           'Apr' 
         { 'May' 'Mai' } 
           'Jun' 
           'Jul' 
           'Aug' 
           'Sep' 
         { 'Oct' 'Okt' }
           'Nov' 
         { 'Dec' 'Dez' }  };


for ii = 1 : size(mon,1)

    mn = mon{ii};

    if ischar(mn);
       mn = {mn};
    end

    for mm = mn
        str = strrep( str , upper(mm{1}) , sprintf('%3.2d',ii) );
        str = strrep( str , sprintf(' %s ',upper(mm{1}([2 3]))) , ...
                            sprintf('%3.2d ',ii) );
    end

end

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
