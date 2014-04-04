function Msg = wrtppm(c,m,file);

% WRTPPM  Writes PPM-ImageFile (truecolor) from ImageData
%
% Indexed   ImageData:  Msg = WRTPPM( CData , ColorMap , file );
%  
% TrueColor ImageData:  Msg = WRTPPM( CData , file );
%
% Multiple Frames can be added 
%    to the 3. Dimension of CData in case of indexed   ImageData
%    to the 4. Dimension of CData in case of truecolor ImageData 
%
%    to the 3. or 4. Dimension of CData in case of grayscaled ImageData 
%
% Description:  Writes image in PPM format, a widely-used, simple, 
%      portable, but non-compressed format.  PPM images can be converted
%      to gif, jpg, tif, bmp, pict, and nearly every other image format
%      know to man (or nerd).   Look for the 'netpbm' distribution on 
%      the internet.
%
%
%  see also:  READPPM, WRTXPM, READXPM, IMREAD, IMWRITE
%
%

Nin = nargin;

Msg = '';

nl  = char(10);

grey = uint8(214);     % NaN-Color

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
Msg0 = sprintf('%s: ',upper(fcn));

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

[mc,c] = chkcdata(c,m);

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

%**************************************************

w = size(c,2);
h = size(c,1);

c = double(permute(c,[3 2 1 4]));

%**************************************************
% Write file

dt = clock;
dt = datenum(dt(1),dt(2),dt(3),dt(4),dt(5),dt(6));
dt = datestr(dt,0);

fid=fopen(file,'w');

if fid == -1
   Msg = sprintf('%sError open File "%s" for writing.',Msg0,file);
   return
end

fprintf(fid,'P6\n');
fprintf(fid,'# %s\n',dt);
fprintf(fid,'# %s: Matlab %s %s\n',fcn,version,computer);

nf = size(c,4);
if nf > 1
   fprintf(fid,'# %.0f Frames\n',nf);
end

fprintf(fid,'%4d %4d\n',w,h);
fprintf(fid,'%4d\n',255);

fwrite(fid,c,'uchar');

fclose(fid);


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
