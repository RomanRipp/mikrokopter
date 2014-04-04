function [z,x,y,z0,pc,c] = read_hgt(file,repl,sub);

% READ_HGT  Reads Topographie from SRTM-HGT-Files
%
% [Z,Lon,Lat,NN,PPH,Mask] = READ_HGT( File )
% [Z,Lon,Lat,NN,PPH,Mask] = READ_HGT( [ LAT LON ] )
%
% A single HGT-Files covers a [ 1° x 1° ] - Tile with
% a Resolution of [ 1201 x 1201 ] (3arcsec). 
%
% The Coordinates [LAT,LON] refers to the SW-Corner of the Tile,
% the FileName is build like "#LAT*LON.hgt", 
%  where "#" is "N" or "S", "*" is "E" or "W".
%
% The Data are in INT16-Format with big endian order (ieee-be).
%
% READ_HGT( ... , N ) calls REPL_NAN( Z , N ) to Replace NaN-Values
%                      default: 0 (No Replace)
%
% READ_HGT( ... , N , SUB ) 
%
% Subsamples Resolution by a negative IntegerValue of SUB
%  or increase the Resolution by a positive Value of SUB.
%  
% If exist for the specific Tile a a Gray-Scaled-PNG ShoreImage 
%  from SRTM_SHORE: #LAT*LON.png; or from SWBD_SHORE: #lon*lat_res.png,
%  with the ColorValues: Land == 000, Ocean == 255, 000 < Lakes < 255;
% the Values on Ocean are set to NaN, the Values on Lakes will averaged.
%
% Special Values for returned Gray-Mask: 
%
%  Red ---> 77 / Green ---> 88 / Blue ---> 99
%
% READ_HGT( ... , N-i ) does not read/use the ShoreImage
%
% READ_HGT( ... , N+1i ) does not average the Lake-Z-Values
% READ_HGT( ... , N+2i ) does not set Ocean-Z-Values
% READ_HGT( ... , N+3i ) does not modify any Z-Values by C
%
% NN   = Mean(Z==Ocean) + i * Median(Z==Ocean)
% PPH  = Percent of Ocean
% Mask = Mask of ShoreImage
%
% see also: SRTM_SHORE, SWBD_SHORE
%

if nargin < 2
   repl = 0;
end

if nargin < 3
   sub = 0;
end

mode = imag(repl);
repl = real(repl);

z0 = NaN;
pc = NaN;

 c = [];

set_lake  = any( mode == [ 0  2 ] );
set_ocean = any( mode == [ 0  1 ] );

if isnumeric(file) & ( prod(size(file)) == 2 )

   ew = 'EW';
   ns = 'NS';

   ext = '.hgt';

   y = file(1);
   x = file(2);

   ns = ns(1+(y<0));
   ew = ew(1+(x<0));

   file = sprintf('%s%2.2d%s%3.3d%s',ns,abs(y),ew,abs(x),ext);

elseif ~chkstr(file,1)

   error('File must be a String or [ Lat Lon ].');

end

%*********************************************
% Check for ZIP-File

file0 = file;

[p,n,e] = fileparts(file);
is_zip = strcmp(e,'.zip');
if ~( exist(file,'file') == 2 ) 
    zf = file;
    if ~is_zip
        zf = cat(2,file,'.zip');
    end
    is_zip = ( exist(zf,'file') == 2 ); 
    if ~is_zip
        [n,p,f] = whichfile(zf);
        is_zip = ~isempty(n);
        if is_zip
           ok = strcmp(n,zf);
           is_zip = any(ok);
           if is_zip
              ok = find(ok);
              zf = f{ok(1)};
           end
        end
    end
    if is_zip
       file = zf;
    end
end

%---------------------------------------------

if ~is_zip & ~( exist(file,'file') == 2 )
    error(sprintf('File "%s" doesn''t exist.',file));
end

f = which(file);
if ~isempty(f)
    file = f;
end

%*********************************************
% Unzip Archive

pfd = '';

if is_zip

   [pfd,tmp] = fileparts(tempname);
    pfd = fullfile(pfd,sprintf('hgt_%s',tmp));
    while ( exist(pfd,'file') == 2 ) | ( exist(pfd,'dir') == 7 )
          pfd = sprintf('%s0',pfd);
    end

    ok = mkdir(pfd);
    if ~ok
        error(sprintf('Cann''t create TempDir "%s" for UNZIP.',pfd));
    end

    cmd = sprintf('unzip -d "%s" "%s"',pfd,file);

    [s,w] = unix(cmd);

    if ~( s == 0 )

        m = sprintf('Error call UNZIP for UNIX: %s',cmd);
        if ~isempty(w)
            m = sprintf('%s\n%s',m,w);
        end

        [s,w] = unix(sprintf('rm -r "%s"',pfd));
        if exist(pfd,'dir') == 7
           m = sprintf('%s\nCann''t remove TempDir: "%s".',m,pfd);
           if ~isempty(w)
               m = sprintf('%s\n%s',m,w);
           end
        end     

        error(m)

    end

    [p,n,e] = fileparts(file);
       file = fullfile(pfd,n);

end

%*********************************************
% Read Z

fid = fopen(file,'r','ieee-be');

if fid == -1
   if ~isempty(pfd), [s,w] = unix(sprintf('rm -r "%s"',pfd)); end
   error(sprintf('Can''t open File "%s".',file));
end

z = fread(fid,'int16');

fclose(fid);

%---------------------------------------------------------------
% Check Bytes of File

d = dir(file);
if ~( prod(size(d)) == 1 )
    d = dir(which(file));
end

%---------------------------------------------------------------
% Clean Up

if ~isempty(pfd), [s,w] = unix(sprintf('rm -r "%s"',pfd)); end

%---------------------------------------------------------------

if ( prod(size(d)) == 1 )
   d = d.bytes;
else
   d = 0;
end


if d == 0
   d = 2*prod(size(z));
   warning(sprintf('Can''t estimate Resolution, use: %.0f',sqrt(d/2)));
end

res = sqrt(d/2);         % 2 Bytes, Square
if ~( mod(res,1) == 0 )
     res = floor(res);
     z = z( 1 : (res*res) );
     warning(sprintf('Resolution not square-defined, use %.0f.',res));
end

 
%*********************************************
% Transform Z

z = reshape(z,res,res);

z = permute(z,[2 1]);

z = z(end:-1:1,:);

z(find(z==-32768)) = NaN;

%*********************************************

rs = res - 1;

xi = ( 0 : rs );

if abs(sub) > 1

   if sub < 0
       ii = ( 1 : -sub : res );
       xi = xi(ii);
        z = z(ii,ii);
   else
       rn = sub * rs;
       xj = linspace(0,rs,rn+1);
        z = interp2( xi ,  xi(:) , z , xj , xj(:) , 'cubic' );  
       xi = (0:rn);
       rs = rn;
      res = rs + 1;
   end

end


%*********************************************
% X and Y

[m,v] = str2vec(file0);

x = v(2) * ( 1 - 2 * any(file=='W') );
y = v(1) * ( 1 - 2 * any(file=='S') );

x = x + xi / rs;
y = y + xi / rs;

%*********************************************
if mode < 0

   c = uint8(zeros(size(z)));

   return

end

%*********************************************
% Check for ShoreImage

[p,n,e] = fileparts(file0);
while ~isempty(e)
       [q,n,e] = fileparts(n);
end

[m,ll] = str2vec(n);

c   = [];
msg = '';

%---------------------------------------------
% Try SWBD-Image
%
% e###n##_**.png / ** == RES[sec] 
%

if isequal(size(ll),[1 2])
   ll = ll([2 1]);  % [ Lon Lat ]
   ew = 'ew';
   ns = 'ns';
   ff = '%s%3.3d%s%2.2d_%2.2d.png';
   f = sprintf( ff , ew(1+any(n=='W')) , abs(ll(1)) , ...
                     ns(1+any(n=='S')) , abs(ll(2)) , round(3600/(res-1)) );
   if exist(f,'file') == 2 
      try
         c = imread(f);
      catch
         c = [];
      end
   end
end

%---------------------------------------------
% Try GMT-Image
%
% n##e###.png 
%

if isempty(c)

   f = sprintf('%s.png',fullfile(p,n));

   if ~exist(f,'file') == 2 
       f = lower(f);
   end

   if exist(f,'file') == 2 
      try
         c = imread(f);
      catch
         msg = sprintf('Cann''t read ShoreImage "%s".\n%s',f,lasterr);
      end
   else
      msg = sprintf('ShoreImage "%s" doesn''t exist.',f);
   end

end

if isempty(msg)
   sc = size(c); sz = size(z);
   if ~isequal( sc([1 2]), sz([1 2]) )
       msg = sprintf('ShoreImage "%s" has %.0f x %.0f, Data requires  %.0f x %.0f.' , ...
                     f , sc([2 1]) , sz([2 1]));
   end
end

if ~isempty(msg)
    warning(msg);
    if repl > 0 
       z = repl_nan(z,repl);
    end
    c = uint8(zeros(size(z)));
    return
end

if size(c,3) == 1
   c = c(:,:,[1 1 1]);
end

% !!! Blue Rivers in SWBD-Images !!!

isb = find( ( c(:,:,1) == 000 ) & ...
            ( c(:,:,2) == 000 ) & ...
            ( c(:,:,3) == 255 )       );

isg = find( ( c(:,:,1) == 000 ) & ...
            ( c(:,:,2) == 255 ) & ...
            ( c(:,:,3) == 000 )       );

isr = find( ( c(:,:,1) == 255 ) & ...
            ( c(:,:,2) == 000 ) & ...
            ( c(:,:,3) == 000 )       );

c = c(:,:,1);

c(isb) = uint8(99);
c(isb) = uint8(88);
c(isb) = uint8(77);

c = c(end:-1:1,:);

if set_lake

   ii = find( ( 0 < c ) & ( c < 255 ) );  % Lakes

   if ~isempty(ii)

       v = sort(double(c(ii)));

       v(find(diff(v,1,1)==0)+1) = [];

       for jj  = 1 : size(v,1)

           kk  = find( c == v(jj) );

  %%% [ floor(mednan(z(kk))) meannan(z(kk)) min(z(kk)) max(z(kk)) ]

          z(kk) = floor(mednan(z(kk)));

       end

   end

end

   kk  = find( c == 255 );

   if isempty(kk)
      pc = 0;
      return
   end

   z0 = meannan(z(kk)) + i * floor(mednan(z(kk)));

   pc = 100 * prod(size(kk)) / prod(size(z));

   if set_ocean
      z(kk) = imag(z0);
   end

   if repl > 0 
      z = repl_nan(z,repl);
   end

   if set_ocean
      z(kk) = NaN;   % Water
   end


%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

    [msg,str] = loadfile(str(2:n2-1),5*2^20,'char');

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

