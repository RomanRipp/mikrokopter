function [Msg,ccm,ctxt] = wrtxpm(c,m,file);

% WRTXPM  Writes XPM-ImageFile (indexed) from ImageData
%
% Indexed   ImageData:  Msg = WRTXPM( CData , ColorMap , file );
%  
% TrueColor ImageData:  Msg = WRTXPM( CData , file );
%
% CData can be a Indexed ColorMatrice, refering to ColorMap
%   with NaN's for Color None; or a TrueColor-Matrice:
%     - a 3-dimensional RGB-ColorMatrice or 
%     - a 2-dimensional GrayScaled Matrice
%    of class:  UINT8, UINT16 or DOUBLE ( 0 <= CData <= 1 ) 
%
% Note: The XPM-Format is recommended for small indexed Images like Icons etc.
%       The conversion of truecolor to indexed CData can spend a lot of Time !!!
%
%   Multiple Frames can be added to the 3. (4.) Dimension of CData.
%
%
% ColorMap  [ R G B ], 0 <= ColorMap <= 1
%
% [Msg,CharacterColorMap,ColorMapText] = WRTXPM( ... )
%
%    returns the Colormap, translated into Characters (without Color "None").
% 
%    CharacterColorMap = [ size(ColorMap,1) by CharacterPerPixel ]
%    ColorMapText      = [ size(ColorMap,1) by NString           ]
%
% WRTXPM uses follwing special Characters:
%
%     '.'   Color NaN    None        comes first in List !!!
%     'o'   Color White  #FFFFFF
%     '#'   Color Black  #000000
%
%   The Characters ' "%\' will not used.
%
% 
%  see also:  READXPM, WRTPPM, READPPM, IMREAD, IMWRITE
%


Nin = nargin;

Msg  = '';

ccm  = [];
ctxt = '';

nl   = char(10);

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
Msg0 = sprintf('%s: ',upper(fcn));

%*****************************************************************

ColorDepth = 2;   % [ 1  |  2  |  4 ]

ColorScale = 16^ColorDepth - 1;

% Special Characters

ColorNaN   = '.';
ColorWhite = 'o';
ColorBlack = '#';

% Forbidden Characters

ColorForbid = ' "%\';

% Comment Characters

cm0 = '/*';
cm1 = '*/';

%*****************************************************************
% Check Inputs

if Nin < 2
  Msg = sprintf('%sNot enough InputArguments.',Msg0);
  return
end

Msg = cell(0,1);

%**************************************************
% Check File

if ischar(m)
   file = m;
   m    = [];
   Nin  = 3;
elseif Nin < 3
   Msg = cat(1,Msg,{'Input File is missing.'});
end

if Nin == 3 
   if ~( ischar(file) & ~isempty(file) & ...
         ( prod(size(file)) == size(file,2) ) );
       Msg = cat(1,Msg,{'File must be a String.'});
   end
end

%**************************************************
% Check CData and ColorMap

[mc,c,m] = chkcdata(c,m);

if ~isempty(mc)
    Msg = cat(1,Msg,{mc});
end

if ~isempty(Msg)
    Msg = sprintf('%s\n',Msg{:});
    Msg = sprintf('%sInvalid Inputs.\n%s',Msg0,Msg);
    return
end

Msg = '';

if isempty(c)
   return
end

%-------------------------------------------------------
% Check for NaN-Colors in ColorMatrice

c(find(isnan(c))) = 0;

is_nan = any( c(:) == 0 );

if is_nan
   c = c + 1;
   m = cat( 1 , NaN*ones(1,3) , m );
end

nc = size(m,1);


%-------------------------------------------------------
% ColorCharacters

%                 a ..  z       A .. Z       0 .. 9

ccm = cat( 2 , ( 97 : 122 ) , ( 65 : 90 ) , ( 48 : 57 ) , ...
               ( 33 :  47 ) , ( 58 : 64 ) , ( 91 : 94 ) , ...
               (123 : 126 ) );
  
% abcdefghijklmnopqrstuvwxyz ABCDEFGHIJKLMNOPQRSTUVWXYZ 0123456789
% !"#$%&'()*+,-./   :;<=>?@   [\]^   {|}~

ccm = char(ccm);

ccm = ccm(:);


% Remove Special and Forbidden Colors

for cf = [ ColorNaN ColorWhite ColorBlack ColorForbid ]

  ccm( find( double(ccm) == double(cf) ) ) = [];

end

nini = size(ccm,1);

ccm = cat( 1 , ColorNaN  , ColorWhite , ColorBlack , ccm );

%--------------------------------
% Check for Special Color

in  = find( isnan(sum(m,2)) );
isn = ~isempty(in);

m(in,:) = NaN;
 
iw  = find( sum(m,2) == 3 );
isw = ~isempty(iw);

ik  = find( sum(m,2) == 0 );
isk = ~isempty(ik);

isp = [ prod(size(in)) prod(size(iw)) prod(size(ik)) ];

sp  = sum(isp);

%--------------------------------
% Get Char_Per_Pixel

cpp = 1;
while  nc-sp > nini^cpp 
  cpp = cpp + 1; 
end

%--------------------------------
% ColorCharacters

ok     = ones(nc,1);
ok(in) = 0;
ok(iw) = 0;
ok(ik) = 0;
ok     = find(ok);

ii     = zeros(nc,cpp);

ii(ok,:) = comb(nini,cpp,nc-sp) + 3;

ii(in,:) = 1;
ii(iw,:) = 2;
ii(ik,:) = 3;

ccm     = ccm(ii);


%--------------------------------------------------------  
% Transform ColorMap to Character

if isn
   m(in,:) = 1;   %  NaN  --->  1   for DEC2HEX
end

m = round([ ones(1,3) ; m ] * ColorScale );
m = m(:);

m = sprintf(['%0' int2str(ColorDepth) 'X'],m);  

% cmap = [ 3*(nc+1) by ColorDepth ]
% cmap = [ R1 ; R2 ; ... ; G1 ; G2 ; ... ; B1 ; ... BN ] 

m = permute(m,[2 1]);  

% [ R1  .. RN , G1 .. GN , B1 .. BN ]
m = reshape(m,ColorDepth,nc+1,3);
m = permute(m,[2 1 3]);

m = reshape(m,nc+1,3*ColorDepth);

m = m(2:nc+1,:);   % Remove First Color, added before dec2hex

 
i1   = ones(nc,1);

c1 = double('"');
c2 = [  9       double('c #') ];
c3 = [ double('",')  10        ];

m = cat( 2 , char(i1*c1)  , ...
             ccm          , ...
             char(i1*c2)  , ...
             m            , ...
             char(i1*c3)         );      

% NaN-Color  -->  None

if isn

   cn = 'None';

   cn = cat( 2 , cn , char(32*ones(1,3*ColorDepth-size(cn,2)+1)) );

   i1 = ones( prod(size(in)) , 1 );

   cn = cat( 2 , char(i1*c1)   , ...
                 ccm(in,:)     , ...
                 char(i1*c2(1:end-1))  , ...
                 cn(i1,:)              , ...
                 char(i1*c3)                 );      

   m(in,:) = cn;

end
 
ctxt = m;

m = permute(m,[2 1]);
m = permute(m(:),[2 1]);               


%****************************************************
% Transform Image to Character

% Matrice for CharacterIndize

nf = size(c,3);   % Number of Frames

sc = size(c);

%-----------------------------------------------
% Multiple Frames

if nf > 1

   c = reshape( permute( c , [ 2 1 3 ] ) , sc(2) , sc(1)*nf );

   c = permute( c , [ 2  1 ] );

end

%-----------------------------------------------

cc = zeros(sc(1)*nf,sc(2)*cpp);

ccm = permute(ccm,[2 1]);   % [ cpp by  nc(+1) ] 

for ii = 1 : nc

 [i1,i2] = find( c == ii );

 if ~isempty(i1)

    jj = ones(size(i1,1),1);
    ip = ones(1,cpp);


   % Indize in Row of cc
    ic = ((i2-1)*cpp) * ip  +  jj * ( 1 : cpp );

   % Indize in cc
    ic = (ic-1)*sc(1)*nf + i1 * ip;

   % Index of Color ii in ccm
    ind = (ii-1)*cpp + (1:cpp);

    cc(ic) = jj * ind;

 end
end


% Matrice of ColorCharacters

cc = ccm(cc);


if ~isequal( size(cc) , [ sc(1)*nf  sc(2)*cpp ] )
  cc = permute(cc,[2 1]);
  if ~isequal( size(cc) , [ sc(1)*nf  sc(2)*cpp ] )
    Msg = [ Msg0 'Can''t determine CharacterMatrice.' nl ...
                 'File: '  file  '  is empty.'               ];
    fclose(fid);
    return
  end
end

i1   = ones(sc(1)*nf,1);

cc = cat( 2 , char(i1*c1)  , ...
              cc           , ...
              char(i1*c3)         );      

%-------------------------------------
% CharacterColorMap

ccm = permute(ccm,[2 1]);

if is_nan

   % Remove added NaN-Color

    ccm(1,:) = [];  
   ctxt(1,:) = [];

end

%-------------------------------------
% Add FrameNumber as Comment

if nf > 1

   ind = ones( (sc(1)+1)*nf , 1 );

    ic = ( 1 : sc(1)+1 : (sc(1)+1)*nf-sc(1) );  % Index for added Commentary

   ind( ic + 1 ) = 0;

   ind = cumsum(ind);

   cc = cc(ind,:);

   cm = cell(nf,2);    % [ Dummy  True ]
   s2 = size(cc,2)-1;  % Without NewLine

   for ii = 1 : nf

       cm{ii,1} = sprintf('%%%3.3d',ii);
       cm{ii,1} = cat( 2 , cm{ii,1} , char( 32*ones(1,s2-size(cm{ii,1},2))) );

       cm{ii,2} = sprintf( '%s Frame %3.0f of %.0f %s' , cm0 , ii , nf ,cm1 );      

       cc(ic(ii),1:s2) = cm{ii,1};

   end

end

%-------------------------------------
% Build String

cc = permute(cc,[2 1]);
cc = permute(cc(:),[2 1]);               

%-------------------------------------
% Replace DummyComment for FrameNumber

if nf > 1

   for ii = 1 : nf

       cc = strrep(cc,cm{ii,1},cm{ii,2});

   end

end

%****************************************************
% Write File

% File starts with:

% /* XPM */
% /* DATE */
% /* MFile VERSION Computer */
% static char *go_left[] = {
% /* width height Ncolor chars_per_pixel Nframe*/
% "    70    70   7       1              1",
% /* colors */
% ". c None",
% "a c #003910",
% /* pixels */

%-------------------------------------

dt = clock;
dt = datenum(dt(1),dt(2),dt(3),dt(4),dt(5),dt(6));
dt = datestr(dt,0);

fid = fopen(file,'wt');

if fid == -1
  Msg = [ Msg0 'Can''t open file: '  file '   for writing.' ];
  return
end

[pfad,name]=fileparts(file);


fprintf(fid,'%s %s %s\n',cm0,'XPM',cm1);
fprintf(fid,'%s %s %s\n',cm0,dt,cm1);
fprintf(fid,'%s %s: Matlab %s %s %s\n',cm0,fcn,version,computer,cm1);
fprintf(fid,'static char *%s[] = {\n',name);
fprintf(fid,'%s width height Ncolor chars_per_pixel Nframe %s\n',cm0,cm1);
fprintf(fid,'"%s ",\n',sprintf('   %3.0f',[ sc([2 1])  nc  cpp  nf]));
fprintf(fid,'%s %s %s\n',cm0,'colors',cm1);
fprintf(fid,'%s',m);
fprintf(fid,'%s %s %s\n',cm0,'pixels',cm1);

fprintf( fid ,  cc(1:end-2) );

fprintf(fid,'\n};\n');

fclose(fid);



%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function c = comb(n,p,l)

% C = COMB(N,P,L);
%
% Combinations of N Elements to P
%
% Returns C with Max. Length L
%

 c = zeros(n,n^(p-1),p);

 c(:,1,1) = (1:n)';
 
 for ii = 2 : p
    z = n^(ii-2);  % RowNumber of existing Combinations
                   %  c(:,1:z,:)
   jj = ( 1 : z*n );

   jj = jj - z*(ceil(jj/z)-1); % [ 1 .. z   1 .. z   1 .. z  ... ]

   c( : , 1 : (z*n) , 1 : (ii-1) ) = c(:,jj,1:(ii-1));

   jj                  =  ones(n,z*n);
   jj(1,:)             = 0;
   jj(1,z*(0:(n-1))+1) = 1;
   jj(1,:)             = cumsum(jj(1,:),2);
   jj                  = cumsum(jj,1);
   jj                  = jj - n * ( jj > n );

   c(:,1:z*n,ii)       = jj;

 end

 c = permute(c,[3 2 1]);  % [ p   by  n^(p-1)  by  n ]
 c = reshape(c,p,n^p);    % [ p   by  n^p ]
 c = permute(c,[2 1]);    % [ n^p by  p   ]

 c = c( 1:min(n^p,l) , [ 1  (p:-1:2) ] );


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,c,m] = chkcdata(c,m,grey)

% CHKCDATA  Checks Input CData (and ColorMap) for valid ImageData
%
% returns RGB-ColorMatrice or indexed ColorMatrice and ColorMap
%
%-----------------------------------------------------------------------
%     RGB-Input,     RGB-Output: [Msg,CData] = CHKCDATA(CData)
% Indexed-Input,     RGB-Output: [Msg,CData] = CHKCDATA(CData,ColorMap)
%
%     RGB-Input, Indexed-Output: [Msg,CData,ColorMap] = CHKCDATA(CData)
% Indexed-Input, Indexed-Output: [Msg,CData,ColorMap] = CHKCDATA(CData,ColorMap)
%
%-----------------------------------------------------------------------
% RGB-Input
%
%   The RGB-Input CData can be:
%     - a 3-dimensional RGB-ColorMatrice or 
%     - a 2-dimensional GrayScaled Matrice
%    of class:  UINT8, UINT16 or DOUBLE ( 0 <= CData <= 1 ) 
%
%   Multiple Frames (Nf) can be added to the 4. (3.) Dimension of CData.
%
%   The Dimension of the Input RGB-CData can be:
%
%     [ Ny  by  Nx  by  3   by  Nf ] RGB-ColorMatrice
%     [ Ny  by  Nx  by  1   by  Nf ] GrayScaledMatrice
%     [ Ny  by  Nx  by  Nf  by  1  ] GrayScaledMatrice, ( Nf ~= 3 )
%
%-----------------------------------------------------------------------
% Indexed-Input
%
%   The Indexed-Input CData should be a 2-dimensional indexed Matrice,
%    refering to ColorMap, which is a 3-Column RGB-Matrice
%    of class:  UINT8, UINT16 or DOUBLE ( 0 <= ColorMap <= 1 ).
%
%   The Class of the indexed CData can be  UINT8, UINT16 or DOUBLE.
%
%   Multiple Frames (Nf) can be added to the 3. or 4. Dimension of CData.
%
%   The Dimension of the Input indexed CData can be:
%
%     [ Ny  by  Nx  by  Nf  by  1  ]  indexed
%     [ Ny  by  Nx  by  1   by  Nf ]  indexed
%   
%-----------------------------------------------------------------------
% RGB-Output
%
%   The RGB-Output returns a 3-dimensional UINT8-ColorMatrice,
%    contains the RGB-ColorValues ( 0 <= C <= 255 ).
%
%   Multiple Frames (Nf) are added to the 4. Dimension:
%
%     [ Ny  by  Nx  by  3   by  Nf ]  UINT8
%
%-----------------------------------------------------------------------
% Indexed-Output
%
%   The Indexed Output returns a 2-dimensional indexed Matrice,
%    refering to ColorMap, which is a 3-Column RGB-Matrice
%    of class DOUBLE ( 0 <= ColorMap <= 1 ).
%
%   Multiple Frames (Nf) are added to the 3. Dimension:
%
%     [ Ny  by  Nx  by  Nf ]  DOUBLE
%
%-----------------------------------------------------------------------
% NaN-Colors
%
% The optional GREY defines the Color to use for NaN-Colors.
%
%      CHKCDATA( ... , GREY )
%
%  GREY must be a single UINT8-Value.
%
%  Colors, which have in all Values of RGB a NaN will set to GREY.
%  Colors, which have any but not all Values NaN will set to ZERO in NaN-Values.
%
%  default for RGB-Output: GREY = uint8(214)
%
%-----------------------------------------------------------------------
%
% See also: IND2RGB, MAT2IND
%

msg = '';

Nin  = nargin;
Nout = nargout;

if Nin == 0
   c = [];
   m = zeros(0,3);
   return
end

if Nin < 2
   m = [];
end

if Nin < 3
   grey = [];
end

msg = cell(0,1);

%**************************************************
% Check GreyColor

if strcmp(class(m),'uint8') & ( prod(size(m)) == 1 )
   grey = m;
   m    = [];
end

is_grey = ~isempty(grey);

if ~is_grey
    grey = uint8(214);
elseif ~( strcmp(class(grey),'uint8') & ( prod(size(grey)) == 1 ) )
    msg = cat(1,msg,{'Value for GREY must be a single UINT8.'});
end
   
%**************************************************
% Check ColorMap

is_map = ~isempty(m);

cm = class(m);

if any(strcmp(cm,{'uint8' 'uint16'}))
   m = double(m);
   if is_map
      p = 8 * ( 1 + strcmp(cm,'uint16') );
      m = m / ( 2^p - 1 );
   end
   cm = 'double';
end

if ~strcmp(cm,'double')
    msg = cat(1,msg,{'ColorMap must be of class DOUBLE, UINT16 or UINT8.'});
elseif ~is_map
    m = zeros(0,3);
elseif ~( ( ndims(m) == 2 ) & ( size(m,2) == 3 ) & all( abs(m(:)-0.5) <= 0.5 ) )
    msg = cat(1,msg,{'ColorMap must have 3 Columns with Values between 0 and 1.'});
end

%**************************************************
% Check CData

cl = class(c);

if ~any(strcmp(cl,{'uint8' 'uint16' 'double'}))
    msg = cat(1,msg,{'CData must be of class DOUBLE, UINT16 or UINT8.'});
elseif ~isempty(c)
    s3 = size(c,3);
    s4 = size(c,4);
    nd = ndims(c);
    if ( nd > 4 ) | ( ( s4 > 1 ) & ~any( s3 == [ 1  3 ] ) )
        str = sprintf('%s\n%s','CData must have max. 4 Dimensions:', ...
                '2(3|4) Indexed | GrayScale  or  3(4) TrueColor (RGB).');
        msg = cat(1,msg,{str});
    elseif strcmp(cl,'double')
        if ~all( isfinite(c(:)) | isnan(c(:)) )
            msg = cat(1,msg,{'Values of CData must be finite or NaN.'});
        end
    end
end

%**************************************************

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    return
end

msg = '';

if isempty(c) | ( Nout == 1 )
   return
end

%*************************************************
% Input: Indexed CData

s31 = ( s3 == 1 );
s41 = ( s4 == 1 );

if ( is_map & ( s31 | s41 ) )

   is_nan = [];

   if ~strcmp(cl,'double');
      c = double(c) + 1;
   else
      is_nan = find(isnan(c));
   end

   nc = size(m,1);

   c = min(max(round(c),1),nc);

   %----------------------------------------------
   % Output: Indexed CData

   if Nout == 3

      if ~isempty(is_nan)
          if is_grey
             m = cat( 1 , m , double(grey([1 1 1]))/255 );
             c(is_nan) = nc+1;
          else
             c(is_nan) = NaN;
          end
      end

      if ~s41
          c = permute(c,[1 2 4 3]);
      end

      return

   end

   %----------------------------------------------
   % Output: TrueColor CData

   m = uint8(round(m*255));

   if ~isempty(is_nan)
       m  = cat(1,m,grey(ones(1,3)));
       nc = nc+1;
       c(is_nan) = nc;
   end

   if ~s31
       c = permute(c,[1 2 4 3]);
   end

   c = m(cat(3,c,c+nc,c+2*nc));

   return

end

%*************************************************
% Check Class of CData

if strcmp(cl,'uint16')
   c  = uint8( round( 255 * double(c) / (  2^16 - 1 ) ) );
   cl = 'uint8';
elseif strcmp(cl,'double')
   if ~all( ( ( 0 <= c(:) ) & ( c(:) <= 1 ) ) | isnan(c(:)) )
       msg = 'TrueColor CData of class DOUBLE must have Values between 0 and 1.';
       return
   end
end

%*************************************************
% Input: GrayScale CData

if ( s4 == 1 ) & ~any( s3 == [ 1  3 ] )
   c = permute(c,[1 2 4 3]);
   s4 = s3;
   s3 =  1;
end

is_gray = ( s3 == 1 );

if is_gray
   c = c(:,:,ones(1,3),:);
end

%*************************************************
% Output: TrueColor CData

if Nout == 2

   if strcmp(cl,'double')
      is_nan = isnan(c);
      c = uint8(round(c*255));
      if any(is_nan(:))
         if is_gray
              is_nan  = find(is_nan);
            c(is_nan) = grey;
         else
                 nan3 = ( sum(is_nan,3) == 3 );
              is_nan  = find(is_nan);
            c(is_nan) = uint8(0);
            if any(nan3(:))
                 is_nan  = find(nan3(:,:,[1 1 1],:));
               c(is_nan) = grey;
            end    
         end
      end
   end

   return

end

%*************************************************
% Output: Indexed CData

if strcmp(cl,'uint8')
   c = double(c) / 255;
end

s1 = size(c,1);
s2 = size(c,2);

c = permute(c,[1 2 4 3]);
c = reshape(c,s1,s2*s4,3);

[c,m] = mat2ind(c);

is_nan = isnan(m);

m = min(max(m,0),1);

c = reshape(c,s1,s2,s4);

if any(is_nan(:))

   m(find(is_nan)) = 0;

   nan3 = ( sum(is_nan,2) == 3 );

   if any(nan3)

      jj = find(nan3);

      if is_grey
         m(jj,:) = double(grey)/255;
      else
         m(jj,:) = [];
         c(find(c==jj)) = NaN;
         c = c - ( c > jj );
      end
 
   end

end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [c,m] = mat2ind(c)

% MAT2IND  Converts Matrice to Indexed
%
% [ Y , Map ] = MAT2IND( X )
%
%   X   numeric Matrice with max. 3 Dimensions
%
%   Y   Indexed Matrice, refering to Values of Map
%
%  Map  Unique Values of X: [ N  by  size(X,3) ] 
%
% Convert back:
%
%   s1 = size(Y,1); s2 = size(Y,2); s3 = size(Map,2);
%   X  = feval(class(Map),zeros(s1,s2,s3));
%   for ii = 1 : s3
%       X(:,:,ii) = reshape(Map(Y,ii),s1,s2);
%   end
%

if ~( isnumeric(c) & ( ndims(c) <= 3 ) )
    error('Input must be numeric with max. 3 Dimensions.');
elseif isempty(c)
    c = zeros(size(c));
    m = zeros(0,size(c,3));
    return
end

cl = class(c);
if ~strcmp(cl,'double')
    c = double(c);
end

if ~all( isfinite(c(:)) | isnan(c(:)) )
    error('Input must have finite Values or NaN.');
end

s1 = size(c,1);
s2 = size(c,2);
s3 = size(c,3);

n = s1 * s2;

c = reshape(c,n,s3);

si = ( 1 : n )';

for ii = s3 : -1 : 1
    [m,jj] = sort(c(si,ii));
       si  = si(jj);
end

c = c(si,:);

is_nan = find(isnan(c));
if ~isempty(is_nan)
    nn = ceil(max(c(:))) + 1;
    c(is_nan) = nn+1;
    is_nan = 1;
end

c(2:n,:) = c(2:n,:) - c(1:(n-1),:);

if is_nan    
    n = ~( sum( abs(c) < 1e3*eps , 2 ) == s3 );
else
    n = ~( sum(     c  < 1e3*eps , 2 ) == s3 );
end

n(1) = 1;

m = cumsum(c(find(n),:),1);

if is_nan
   m( find( m > nn ) ) = NaN;
end

n = cumsum(n);

c = reshape(n,s1,s2);

c(si) = c;

if ~strcmp(cl,'double')
    m = feval(cl,m);
end

