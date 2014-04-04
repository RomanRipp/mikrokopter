function [c,msg] = cutimage(c,varargin);

% CUTIMAGE removes Margin from Image
%
% CUTIMAGE( CData , ColorLimit , i*PixelOffset )
%
% CData can be:
%   - a TrueColor-CDataMatrice, 
%   - a Structure, returned from GETFRAME,
%   - a FigureHandle, which Image will captured using GETFRAME, or 
%   - a Name of an ImageFile, which ImageData will read using IMREAD.
%
% ColorLimit defines the surrounding Colors to remove, Range: 0 .. 255
%
% CUTIMAGE removes surrounding Rows/Columns of Image, in which
%          ALL Pixels have a mean ColorValue 
%           - larger  if positive ColorLimit, or
%           - smaller if negative ColorLimit 
%          then ColorLimit 
%
%   In case ColorLimit is a imaginary Number with a NonZero Part:
%
%      R + i*ColorLimit
%
%   all surrounding Rows/Columns of Image will removed, in which 
%   ALL Pixels have the same mean ColorValue then ColorLimit.
%
% PixelOffset sets the Number of Rows/Columns, which will added
%              to the cutted Image and will not removed!
%   Give this Input as imaginary Number, use NaN for ZERO  PixelOffset.
%
%
% [C,Msg] = CUTIMAGE( ... ) returns the cutted Image in a UINT8-CDataMatrice,
%                             and ErrorMessages.
%
%-------------------------------------------------------------
% Save the cutted Image into File
%
% CUTIMAGE( ... , OutFile , '.Format' , Options , ... )
%
% Saves the cutted Image to OutFile, using IMWRITE:
%
%   IMWRITE( C , OutFile , Format , Options{:} )
%
% Format can be any of:
%
%    'tif' or 'tiff' Tagged Image File Format (TIFF)
%    'jpg' or 'jpeg' Joint Photographic Experts Group (JPEG)
%    'bmp'           Windows Bitmap (BMP)
%    'png'           Portable Network Graphics (PNG)
%    'hdf'           Hierarchical Data Format (HDF)
%    'pcx'           Windows Paintbrush (PCX)
%    'xwd'           X Window Dump (XWD)
%
%    default: 'png'
%
% Options is a CellArray with Options for IMWRITE.
%  Type  >> help imwrite 
%  for more Informations about valid Options.
%
% Example:
%
%    figure('color','w'), plot(1:10)
%
%    cutimage(gcf,250,5*i,'test.jpg','.jpg',{'Quality' 90})
%
%  saves the current Figure with 5-Pixel white Border
%    to File "test.jpg" with Quality 90%.
%
%-------------------------------------------------------------
%
% see also: IMREAD, IMWRITE, SAVEIMG
%

pix = 10;
val = 250;
img = '';
fmt = '';
opt = cell(1,0);

Nin  = nargin;
Nout = nargout;

if Nin < 1
   error('Not enough InputArguments.');
end

%**********************************************************
% Check Inputs

msg = cell(0,1);

if isempty(c)
   return
end

ok = ischar(c);

if ok  

   ok = ( ~isempty(c) & ( prod(size(c)) == size(c,2) )  );

   if ok
      ok = ( exist(c,'file') == 2 );
      if ~ok
          msg = cat(1,msg,{sprintf('File "%s" doesn''t exist.',c)});
      else
          try
             c = imread(c);
          catch
             ok  = 0;
          end
          if ~ok
              str = sprintf('Error call IMREAD( %s ).\n%s',c,lasterr);
              msg = cat(1,msg,{str});
          end
      end
   end

else

  ok = isnumeric(c);

  if ok

     if prod(size(c)) == 1
        ok = ishandle(c);
        if ok
           ok = strcmp(get(c,'type'),'figure');
        end
        if ok
           figure(c); drawnow; c = getframe(c); 
        else
           msg = cat(1,msg,{'Single Input must be a FigureHandle.'});
        end
      end

   else 

      ok = ( isstruct(c) & ( prod(size(c)) == 1 ) );

   end

end

%-------------------------------------------------------------
% CData

if ~ok

    msg = cat(1,msg,{'First Input must be a nonempty String, Numeric or Structure.'});

else

   [m,c] = chkcdata(c);

   if ~isempty(m)
       msg = cat(1,msg,{sprintf('Invalid CData.\n%s',m)});
   end

end

%-------------------------------------------------------------
% Parameter

for ii = 1 : Nin-1
    v = varargin{ii};
    if chkstr(v,1)
       if v(1) == '.'
          fmt = v(2:end);
       else
          img = v;
       end
    elseif iscell(v)
       opt = cat(2,opt,v(:)');
    elseif isnumeric(v) & ( prod(size(v)) == 1 )
       if isnan(v)
          pix = 0;
       elseif ( real(v) == 0 ) & ~( imag(v) == 0 )
          pix = imag(v);
       else
          val = v;
       end
    else
       str = sprintf('Unrecognized %.0f. Input.',ii+1);
       msg = cat(1,msg,{str});
    end
end

%-------------------------------------------------------------
% Check File

if ~isempty(img)

    if isempty(fmt)
       fmt = 'png';
    end

    [m,name,img] = chk_name(img,'',fmt,'','none');

    if ~isempty(m)
        msg = cat(1,msg,{sprintf('Invalid FileName: %s\n%s',img,m)});
    end

end

%-------------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    if Nout < 2
       error(msg);
    end
    return
end

%**********************************************************

s = sum(double(c),3) / 3;

s1 = size(s,1);
s2 = size(s,2);

if ( imag(val) == 0 )
   s = ( sign(val)*s > val );
else
   s = ( s == imag(val) );
end

i2  = ( sum(s,1) == s1  );
i1  = ( sum(s,2) == s2  );

i2 = [  sum(cumprod(i2))-1  s2-sum(cumprod(i2(s2:-1:1))) ];
i2 = i2 + [ -1  1 ] * pix;

i2(1) = max(i2(1),1);
i2(2) = min(i2(2),s2);

i1 = [  sum(cumprod(i1))-1  s1-sum(cumprod(i1(s1:-1:1))) ];
i1 = i1 + [ -1  1 ] * pix;

i1(1) = max(i1(1),1);
i1(2) = min(i1(2),s1);

c = c( i1(1) : i1(2) , i2(1) : i2(2) , : );

if isempty(img)
   return
end

if isempty(c)
   fprintf(1,'%s\nWarning: Image empty, no File saved.\n',char(7));
   return
end

%-------------------------------------------------------

msg = '';

fprintf(1,'\nWrite %s ... ',img);

try
   imwrite(c,img,fmt,opt{:});
catch
   msg = lasterr;
end

if isempty(msg)
   fprintf(1,'ok\n');
   return
end

fprintf(1,'error\n');

msg = sprintf('Error call IMWRITE( %s ).\n%s',img,msg);

if Nout < 2
   error(msg);
end

%********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
% Check Structure

if isstruct(c)

   f = fieldnames(c);
   l = lower(f);

   ic = strcmp(l,'cdata');
   if ~any(ic)
       ic  = [];
       msg = cat(1,msg,{'StructureInput must have the FieldName "cdata".'});
   else
       ic = find(ic);
       ic = ic(1);
   end

   im = strcmp(l,'colormap');
   if ~any(im)
       im = strcmp(l,'cmap');
   end

   if ~any(im)
       im = [];
   else
       im = find(im);
       im = im(1);
   end

   if ~isempty(im)
        m = getfield(c,f{im});
   end

   if ~isempty(ic)
        c = getfield(c,f{ic});
   end

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

   m = uint8(m*255);

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
   c  = uint8( 255 * double(c) / (  2^16 - 1 ) );
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
      c = uint8(c*255);
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

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,name,file] = chk_name(name,pfad,ext,pre,chk);

% CHK_NAME   Checks FileName of valid Characters, Convention and Directory
%
% [ Msg , Name , FullName ] = CHK_NAME( Name , Pfad , Extension , Prefix , Check )
%
%    Check = 'new'    Warning, if File or Directory exists 
%            'none'   no check for existing File/Directory
%            'file'   Warning, if File not exist 
%            'dir'    Warning, if File not exist 
%


Msg = '';

nl = char(10);

file = '';

%-------------------------------------------------------------

if nargin < 2
  pfad = '';
end

if nargin < 3
  ext = '';
end

if nargin < 4
  pre = '';
end

if nargin < 5
 chk = 'none';
end

%-------------------------------------------------------------

if ~( ischar(pfad)  &  ( prod(size(pfad)) == size(pfad,2) )  )
   Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
           'Input Pfad must be a String.'      ];
end

if ~( ischar(name)  &  ( prod(size(name)) == size(name,2) )  )
   Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
           'Input Name must be a String.'      ];
end

if ~( ischar(ext)  &  ( prod(size(ext)) == size(ext,2) )  )
   Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
           'Input Extension must be a String.'      ];
end

if ~( ischar(pre)  &  ( prod(size(pre)) == size(pre,2) )  )
   Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
           'Input Prefix must be a String.'      ];
end

if ~( ischar(chk)  &  ( prod(size(chk)) == size(chk,2) )  )
   Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
           'Input Check must be a String.'      ];
end

if ~isempty(Msg)
   return
end

%-------------------------------------------------------------

name = rmblank(name,2);

  fs = filesep;

if ~isempty(pfad)
    pfad = pfad( 1 : (end-strcmp(pfad(end),fs)) );
end


if ~isempty(name)

  jj = find( name == fs );

  if ~isempty(jj)

      jj = max(jj);

      if isempty(pfad)
         pfad = name( 1 : (jj-1) );
      else
         if ~strcmp(pfad,name(1:(jj-1)));
             Msg = sprintf('Directory "%s" is not conform with "%s".',pfad,name(1:(jj-1)));
             return
         end
      end
      
      name = name( (jj+1) : end );

   end

end


if ~isempty(name)

  %---------------------------------------------------
  % Check for valid Characters
  % 0 .. 9  |  A .. Z   |  a .. z  |  .  |   _  |  FileSeparator

  name = double(name); 
    
    fs = double(fs);

  ii = find( ~( (  48 <= name  &  name <= 57  )  | ...
                (  65 <= name  &  name <= 90  )  | ...
                (  97 <= name  &  name <= 122 )  | ...
                  name == 46   |  name == 95     |  name == 126  | ...
                ( name == fs   &  strcmp(chk,'dir') )   ) );

  name = char(name);

  if ~isempty(ii)
     Msg = 'Invalid Character''s in Name.';
     return
  end

  %---------------------------------------------------
  % Append pre

  if ~isempty(pre)

    s2 = size(name,2);

    sp = size(pre ,2);

    if ~strcmp(lower(pre),lower(name(1:min([sp s2]))))
       name = cat( 2 , pre , name );
    end

  end

  %---------------------------------------------------
  % Append ext

  if ~isempty(ext)


    f  = '.';

    ext = cat( 2 , f(1:(end*(~strcmp(ext(1),f)))) , ext );

    s2 = size(name,2);

    se = size(ext,2);

    if ~strcmp(lower(ext),lower(name(s2-min([se s2])+1:s2)))
       name = cat( 2 , name , ext( (1+strcmp(name(s2),f)) : se ) );
    end
    
  end
  
end 


%---------------------------------------------------
% Append pfad

if ~isempty(pfad)

   pfad = cat( 2 , pfad , fs );

   if isempty(name)
      file = pfad;
   else
      file = cat( 2 , pfad , name );
   end

else

   file = name;

end


%-----------------------------------------------------
% Check

isfile = ( exist(file,'file') == 2 );
isdir  = ( exist(file,'dir')  == 7 );

chk = lower(chk);

switch chk

   %---------------------------------------------------
   case 'new'

     if isfile
        Msg = [ 'File allready exist: ' file ];
     elseif isdir
        Msg = [ 'Directory allready exist: ' file ];
     end
   
   %---------------------------------------------------
   case 'file'

     if ~isfile
        Msg = [ 'File does not exist: ' file ];
     end

     if isdir
        Msg = [ 'Directory exist: ' file ];
     end

   %---------------------------------------------------
   case 'dir' 

     if ~isdir
        Msg = [ 'Directory does not exist: ' file ];
     end

     if isfile
        Msg = [ 'File exist: ' file ];
     end

end

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

