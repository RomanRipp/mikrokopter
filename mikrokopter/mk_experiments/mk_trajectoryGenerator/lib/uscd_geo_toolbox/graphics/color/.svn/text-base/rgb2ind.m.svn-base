function [X,map] = rgb2ind(varargin),

%RGB2IND Convert RGB image to indexed image.
%   RGB2IND converts RGB images to indexed images using one of three
%   different methods: uniform quantization, minimum variance quantization,
%   and colormap approximation. RGB2IND dithers the image unless you specify
%   'nodither' for DITHER_OPTION.
%
%   [X,MAP] = RGB2IND(RGB,N) converts the RGB image to an indexed image X
%   using minimum variance quantization. MAP contains at most N colors.  N
%   must be <= 65536.
%
%   X = RGB2IND(RGB,MAP) converts the RGB image to an indexed image X with
%   colormap MAP by matching colors in RGB with the nearest color in the
%   colormap MAP.  SIZE(MAP,1) must be <= 65536.
%
%   [X,MAP] = RGB2IND(RGB,TOL) converts the RGB image to an indexed image X
%   using uniform quantization. MAP contains at most (FLOOR(1/TOL)+1)^3
%   colors. TOL must be between 0.0 and 1.0.
%
%   [...] = RGB2IND(...,DITHER_OPTION) enables or disables
%   dithering. DITHER_OPTION is a string that can have one of these values:
%
%       'dither'   dithers, if necessary, to achieve better color
%                  resolution at the expense of spatial
%                  resolution (default)
%
%       'nodither' maps each color in the original image to the
%                  closest color in the new map. No dithering is
%                  performed.
%
%   Class Support
%   -------------
%   The input image can be of class uint8, uint16, or double. The output
%   image is of class uint8 if the length of MAP is less than or equal to
%   256, or uint16 otherwise.
%
%   Example
%   -------
%       RGB = imread('flowers.tif');
%       [X,map] = rgb2ind(RGB,128);
%       imshow(X,map)
%
%   See also CMUNIQUE, DITHER, IMAPPROX, IND2RGB, RGB2GRAY.

% Grandfathered syntax
% --------------------
%   [X,MAP] = RGB2IND(RGB) converts the RGB image in the array
%   RGB to an indexed image X with colormap MAP using direct
%   translation. The resulting colormap may be very long, as it
%   has one entry for each pixel in RGB. Do not set DITHER_OPTION
%   if you use this method.
%

%   Copyright 1993-2002 The MathWorks, Inc.  
%   $Revision: 5.21 $  $Date: 2002/03/15 15:29:05 $

[RGB,m,dith] = parse_inputs(varargin{:});

[so(1) so(2) thirdD] = size(RGB);

% Converts depending on what is m:
if isempty(m),% Convert RGB image to an indexed image.
    X = reshape([1:so(1)*so(2)],so(1),so(2));
    if so(1)*so(2) <= 256
        X = uint8(X-1);
    elseif so(1)*so(2) <= 65536
        X = uint16(X-1);
    end
    map = im2double(reshape(RGB,so(1)*so(2),3));

elseif length(m)==1,% TOL or N is given
    RGB = im2uint8(RGB);

    if m<1,% tol is given. Use uniform quantization
        max_colors = 65536;
        max_N = floor(max_colors^(1/3)) - 1;
        N = round(1/m);
        if (N > max_N)
            N = max_N;
            warning(sprintf('Too many colors; increasing tolerance to %g',...
                            1/N));
        end
        
        [x,y,z] = meshgrid([0:N]/N);
        map = [x(:),y(:),z(:)];

        if dith(1) == 'n'; 
            RGB = round(im2double(RGB)*N);
            X = RGB(:,:,3)*((N+1)^2)+RGB(:,:,1)*(N+1)+RGB(:,:,2)+1;
        else
            X = dither(RGB,map);
        end
        [X,map] = cmunique(X,map);

    else % N is given. Use variance minimization quantization
        [map,X] = cq(RGB,m);
        map = double(map) / 255;
        if dith(1)=='d',% Use standalone dither if map is an approximation.
            X = dither(RGB,map);
        end
    end

else % MAP is given
    RGB = im2uint8(RGB);

    map = m;
    if dith(1)=='n', % 'nodither'
        X = dither(RGB,map,5,4); % Use dither to do inverse colormap mapping.
    else
        X = dither(RGB,map);
    end
end

if isa(map, 'uint8'),    % Make sure that the colormap is doubles
    map = double(map)/255;
elseif isa(map, 'uint16')
    map = double(map)/65535;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function: parse_inputs
%

function [RGB,m,dith] = parse_inputs(varargin)
% Outputs:  RGB     image
%           m       colormap
%           dith    dithering option
% Defaults:
dith = 'dither';
m = [];

error(nargchk(1,5,nargin));
switch nargin
case 1,               % rgb2ind(RGB)
  RGB = varargin{1};
  warning(sprintf('%s\n%s', 'RGB2IND(RGB) is an obsolete syntax.', ...
                  'Specify number of colors, tolerance, or colormap.'));
case 2,               % rgb2ind(RGB,x) where x = MAP | N | TOL
  RGB = varargin{1};
  m = varargin{2};
case 3,               % rgb2ind(R,G,B) OBSOLETE
  if isequal(size(varargin{1}),size(varargin{2}),size(varargin{3})),
    warning(['RGB2IND(R,G,B) is an obsolete syntax. ',...
    'Use a three-dimensional array to represent RGB image.']);
    RGB = cat(3,varargin{1},varargin{2},varargin{3});
  else                % rgb2ind(RGB,x,DITHER_OPTION)
    RGB = varargin{1};  %              where x = MAP | N | TOL
    m = varargin{2};
    dith = varargin{3};
  end;
case 4,               % rgb2ind(R,G,B,x) OBSOLETE
  warning(['RGB2IND(R,G,B,x) is an obsolete syntax. ',...
  'Use a three-dimensional array to represent RGB image.']);
  if isequal(size(varargin{1}),size(varargin{2}),size(varargin{3})),
    RGB = cat(3,varargin{1},varargin{2},varargin{3});
  else
    error('R,G,B arrays must be of equal size.');
  end;
  m = varargin{4};
case 5,               % rgb2ind(R,G,B,x,DITHER_OPTION) OBSOLETE
  warning(['RGB2IND(R,G,B,x) is an obsolete syntax. ',...
  'Use a three-dimensional array to represent RGB image.']);
  if isequal(size(varargin{1}),size(varargin{2}),size(varargin{3})),
    RGB = cat(3,varargin{1},varargin{2},varargin{3});
  else
    error('R,G,B arrays must be of equal size.');
  end;
  m = varargin{4};
  dith = varargin{5}; 
otherwise,
  error('Invalid input arguments.');
end

% Check validity of the input parameters
if (ndims(RGB)==3)&(size(RGB,3) ~= 3)|(ndims(RGB)>3),
  error('RGB image has to be M-by-N-by-3 array.');
end;

if any(m(:)<0),
  error('Colormap, Number of colors, or Tolerance have to be non-negative.');
elseif (length(m)==1)&(m~=round(m))&(m>1),
  error('Number of colors in the colormap has to be a non-negative integer.');
elseif (length(m)>1)&(ndims(m)==2)&(size(m,1)==1)&(size(m,2)~=3),% MAP
  error('Input colormap has to be a 2D array with at least 2 rows and exactly 3 columns.');
elseif (length(m)>1)&(max(m(:))>1),%
  error('All colormap intensities must be between 0 and 1.');
end;

if ischar(dith),% dither option
  strings = {'dither','nodither'};
  idx = strmatch(lower(dith),strings);
  if isempty(idx),
    error(sprintf('Unknown dither option: %s',dith));
  elseif length(idx)>1,
    error(sprintf('Ambiguous dither option: %s',dith));
  else
    dith = strings{idx};
  end  
else
  error(sprintf('Dither option has to be a string.'));  
end;

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function d = im2double(img, typestr)
%IM2DOUBLE Convert image to double precision.
%   IM2DOUBLE takes an image as input, and returns an image of
%   class double.  If the input image is of class double, the
%   output image is identical to it.  If the input image is of
%   class logical, uint8 or uint16, im2double returns the 
%   equivalent image of class double, rescaling or offsetting
%   the data as necessary.
%
%   I2 = IM2DOUBLE(I1) converts the intensity image I1 to double
%   precision, rescaling the data if necessary.
%
%   RGB2 = IM2DOUBLE(RGB1) converts the truecolor image RGB1 to
%   double precision, rescaling the data if necessary.
%
%   BW2 = IM2DOUBLE(BW1) converts the binary image BW1 to double
%   precision.
%
%   X2 = IM2DOUBLE(X1,'indexed') converts the indexed image X1 to
%   double precision, offsetting the data if necessary.
% 
%   See also DOUBLE, IM2UINT8, UINT8.

%   Copyright 1993-2002 The MathWorks, Inc.  
%   $Revision: 1.16 $  $Date: 2002/03/15 15:27:36 $

checknargin(1,2,nargin,mfilename);
checkinput(img,{'double','logical','uint8','uint16'},{},mfilename,'Image',1);

if isa(img, 'double')
   d = img;
elseif isa(img, 'logical')
   d = double(img);
elseif isa(img, 'uint8') | isa(img, 'uint16')
   if nargin==1
      if isa(img, 'uint8')
          d = double(img)/255;
      else
          d = double(img)/65535;
      end
   elseif nargin==2
      checkstrs(typestr,{'indexed'},mfilename,'type',2);
      d = double(img)+1;
   end
end

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function u = im2uint8(varargin)
%IM2UINT8 Convert image to eight-bit unsigned integers.
%   IM2UINT8 takes an image as input, and returns an image of
%   class uint8.  If the input image is of class uint8, the
%   output image is identical to it.  If the input image is of
%   class logical, double or uint16, im2uint8 returns the
%   equivalent image of class uint8, rescaling or offsetting 
%   the data as necessary.
%
%   I2 = IM2UINT8(I1) converts the intensity image I1 to uint8,
%   rescaling the data if necessary.
%
%   RGB2 = IM2UINT8(RGB1) converts the truecolor image RGB1 to
%   uint8, rescaling the data if necessary.
%
%   BW2 = IM2UINT8(BW1) converts the binary image BW1 to uint8.
%
%   X2 = IM2UINT8(X1,'indexed') converts the indexed image X1 to
%   uint8, offsetting the data if necessary. Note that it is
%   not always possible to convert an indexed image to
%   uint8. If X1 is double, then the maximum value of X1 must
%   be 256 or less.  If X1 is uint16, the maximum value of X1
%   must be 255 or less.
%
%   See also DOUBLE, IM2DOUBLE, IM2UINT16, UINT8, UINT16, IMAPPROX.

%   Copyright 1993-2002 The MathWorks, Inc.  
%   $Revision: 1.20 $  $Date: 2002/03/15 15:27:39 $

checknargin(1,2,nargin,mfilename);

img = varargin{1};
checkinput(img,{'double','logical','uint8','uint16'},{},mfilename,'Image',1);

if isa(img, 'uint8')
    u = img; 
elseif isa(img, 'logical')
    u=uint8(img);
    u(img)=255;
elseif isa(img, 'double') | isa(img, 'uint16')
    if nargin==1
         % intensity image; call MEX-file
         u = grayto8(img);
    elseif nargin==2
       typestr = varargin{2};
       checkstrs(typestr,{'indexed'},mfilename,'type',2);
       if (isa(img, 'uint16'))
            if (max(img(:)) > 255)
                msg = 'Too many colors for 8-bit integer storage.';
                eid = sprintf('Images:%s:tooManyColorsFor8bitStorage');
                error(eid,msg);
            else
                u = uint8(img);
            end
        else
            % img is double
            if max(img(:))>=257 
                msg = 'Too many colors for 8-bit integer storage.';
                eid = sprintf('Images:%s:tooManyColorsFor8bitStorage');
                error(eid,msg);
            elseif min(img(:))<1
                msg = 'Invalid indexed image: an index was less than 1.';
                eid = sprintf('Images:%s:invalidIndexedImage');
                error(eid,msg);
            else
                u = uint8(img-1);
            end
        end
    end
end

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [c,map] = cmunique(varargin)
%CMUNIQUE Find unique colormap colors and corresponding image.
%   [Y,NEWMAP] = CMUNIQUE(X,MAP) returns the indexed image Y and 
%   associated colormap NEWMAP that produce the same image as
%   (X,MAP) but with the smallest possible colormap. CMUNIQUE
%   removes duplicate rows from the colormap and adjusts the
%   indices in the image matrix accordingly.
%
%   [Y,NEWMAP] = CMUNIQUE(RGB) converts the truecolor image RGB
%   to the indexed image Y and its associated colormap
%   NEWMAP. NEWMAP is the smallest possible colormap for the
%   image, containing one entry for each unique color in
%   RGB. (Note that NEWMAP may be very large, as much as P-by-3
%   where P is the number of pixels in RGB.) 
%
%   [Y,NEWMAP] = CMUNIQUE(I) converts the intensity image I to an
%   indexed image Y and its associated colormap NEWMAP. NEWMAP is
%   the smallest possible colormap for the image, containing one
%   entry for each unique intensity level in I. 
%
%   Class Support
%   -------------
%   The input image can be of class uint8, uint16, or double. 
%   The class of the output image Y is uint8 if the length of 
%   NEWMAP is less than or equal to 256. If the length of 
%   NEWMAP is greater than 256, Y is of class double.
%
%   See also RGB2IND, GRAY2IND.

%   Copyright 1993-2002 The MathWorks, Inc.  
%   $Revision: 5.18 $  $Date: 2002/03/15 15:26:49 $

%   I/O Spec
%   ========
%   IN
%      X      - image of class uint8, uint16, or double
%      MAP    - M-by-3 array of doubles (colormap)
%   OUT
%      Y      - uint8 if NEWMAP has <= 256 entries, double 
%               if NEWMAP has > 256 entries.
%      NEWMAP - M-by-3 array of doubles (colormap)

checknargin(1,3,nargin,mfilename);
checkinput(varargin{1},{'double' 'uint8' 'uint16'},{'real' 'nonsparse'}, mfilename,'X',1);

% Convert all possible input arguments to an indexed image.
if nargin==1, % cmunique(I) or cmunique(RGB)
    arg1 = varargin{1};
    if ndims(arg1)==3, % cmunique(RGB)
        [c,map] = rgb2ind(arg1);
    else % cmunique(I)
        [c,map] = rgb2ind(arg1,arg1,arg1);
    end
elseif nargin==2, % cmunique(a,cm)
    c = varargin{1}; map = varargin{2};
elseif nargin==3, % cmunique(r,g,b)
    warning('Images:cmunique:obsoleteSyntax',['CMUNIQUE(r,g,b) is an obsolete syntax.',...
    'Use a three dimensional array to represent RGB image.']);
    [c,map] = rgb2ind(varargin{1},varargin{2},varargin{3});
end

if ~isa(c, 'double')    % The promotion is necessary for the indexing into
    c = im2double(c, 'indexed');  % pos below --  ...loc(pos(c))...
end

tol = 1/1024;

% Quantize colormap entries to help matching below.
map = round(map/tol)*tol;

%  
% Remove matching entries from colormap
%
 
% Sort colormap entries
[dum,ndx1] = sort(map(:,1));
[dum,ndx2] = sort(map(ndx1,2));
[dum,ndx3] = sort(map(ndx1(ndx2),3));
                % ndx maps from sorted cm to original cm
ndx = ndx1(ndx2(ndx3));
                % pos maps from original cm to sorted cm
pos = zeros(size(ndx)); pos(ndx) = 1:length(ndx);

% Find matching entries
                % d indicates the location of matching entries
d = all(abs(diff(map(ndx,:)))'<tol)';

% Mapping from full cm to compressed cm
                % loc maps from sorted cm to compressed cm
loc = [1:length(ndx)]' - [0;cumsum(d)]; 
c(:) = loc(pos(c));

% Remove matching entries (compress cm)
ndx(find(d)) = [];
map = map(ndx,:);

%
% Remove colormap entries that are not used in c
%
[n,x] = imhist(c,map);
d = (n==0); % Find unused colormap entries
loc = [1:length(map)]' - cumsum(d);

% Update image values
c(:) = loc(c);

% Remove unused entries (compress cm)
map(find(d),:) = [];

if max(c(:))<=256    % Output a uint8 array if we can
    c = uint8(c-1);
end

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function varargout=dither(varargin)
%DITHER Convert image using dithering.
%   X = DITHER(RGB,MAP) creates an indexed image approximation of
%   the RGB image in the array RGB by dithering the colors in
%   colormap MAP.  MAP cannot have more than 65536 colors.
%
%   X = DITHER(RGB,MAP,Qm,Qe) creates an indexed image from RGB,
%   specifying the parameters Qm and Qe. Qm specifies the number
%   of quantization bits to use along each color axis for the
%   inverse color map, and Qe specifies the number of
%   quantization bits to use for the color space error
%   calculations.  If Qe < Qm, dithering cannot be performed and
%   an undithered indexed image is returned in X.  If you omit
%   these parameters, DITHER uses the default values Qm = 5, 
%   Qe = 8.
%
%   BW = DITHER(I) converts the intensity image in the matrix I
%   to the binary image BW by dithering.
%
%   Class Support
%   -------------
%   The input image (RGB or I) can be of class uint8, uint16, or
%   double. All other input arguments must be of class
%   double. The output image (X or BW) is of class logical if it is
%   a binary image or uint8 if it is an indexed image with 256 or fewer
%   colors; otherwise its class is uint16.
%
%   See also RGB2IND.

%   Copyright 1993-2002 The MathWorks, Inc.  
%   $Revision: 5.24 $  $Date: 2002/03/15 15:27:09 $

%   References: 
%      R. W. Floyd and L. Steinberg, "An Adaptive Algorithm for
%         Spatial Gray Scale," International Symposium Digest of Technical
%         Papers, Society for Information Displays, 36.  1975.
%      Spencer W. Thomas, "Efficient Inverse Color Map Computation",
%         Graphics Gems II, (ed. James Arvo), Academic Press: Boston.
%         1991. (includes source code)

[X,m,qm,qe] = parse_inputs_dither(varargin{:});

if ndims(X)==2,% Convert intensity image to binary by dithering
  im = logical(ditherc(X,m,qm,qe));
  map = 2;
else % Create an indexed image from RGB 
  im = ditherc(X,m,qm,qe);
  map = m;
end

if nargout==0,
  imshow(im,map);
else 
  varargout{1} = im;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function: parse_inputs_dither
%

function [X,m,qm,qe] = parse_inputs_dither(varargin)
% Outputs:  X  the input RGB (3D) or intensity image (2D)
%           m  colormap (:,3)
%           qm number of quantization bits for colormap
%           qe number of quantization bits for errors, qe>qm

checknargin(1,6,nargin,mfilename);

% Default values:
qm = 5;
qe = 8;

switch nargin
case 1,                        % dither(I)
  X = varargin{1};
  m = gray(2); 
case 2,                        % dither(RGB,m)
  X = varargin{1};
  m = varargin{2};
case 4,
  if length(varargin{4})==1,   % dither(RGB,m,qm,qe)
    X = varargin{1};
    m = varargin{2};
    qm = varargin{3};
    qe = varargin{4};
  else                        % dither(R,G,B,m) OBSOLETE
    eid = sprintf('Images:%s:obsoleteSyntax',mfilename);
    msg = ['DITHER(R,G,B,m) is an obsolete syntax. ',...
           'Use a three-dimensional array to represent RGB image.'];
    warning(eid,msg);
    X = cat(3,varargin{1},varargin{2},varargin{3});
    m = varargin{4};
  end;
case 6,                       % dither(R,G,B,m,qm,qe) OBSOLETE
  eid = sprintf('Images:%s:obsoleteSyntax',mfilename);
  msg = ['DITHER(R,G,B,m,qm,qe) is an obsolete syntax. ',...
         'Use a three-dimensional array to represent RGB image.']
  warning(eid, msg);
  X = cat(3,varargin{1},varargin{2},varargin{3});
  m = varargin{4};
  qm = varargin{5};
  qe = varargin{6};
otherwise,
  eid = sprintf('Images:%s:invalidInput',mfilename);
  error(eid,'Invalid input arguments in function %s.',mfilename);
end

% Check validity of the input parameters 
if (ndims(X)==3)&(nargin==1),
  eid = sprintf('Images:%s:imageMustBe2D',mfilename);  
  error(eid,'DITHER(I): the intensity image I has to be a two-dimensional array.');
elseif (ndims(X)==2)&(nargin==2),
  eid = sprintf('Images:%s:imageMustBe3D',mfilename);  
  error(eid,'DITHER(RGB,map): the RGB image has to be a three-dimensional array.');
end;

X = im2uint8(X);
 
if (size(m,2) ~= 3)|(size(m,1)==1)|ndims(m)>2,
  eid = sprintf('Images:%s:colormapMustBe2D',mfilename);  
  error(eid,['In function %s, input colormap has to be a ',...
             '2D array with at least 2 rows and exactly 3 columns.'], mfilename);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function checkinput(A, classes, attributes, function_name, ...
                    variable_name, argument_position)
%CHECKINPUT Check validity of array.
%   CHECKINPUT(A,CLASSES,ATTRIBUTES,FUNCTION_NAME,VARIABLE_NAME, ...
%   ARGUMENT_POSITION) checks the validity of the array A and issues a
%   formatted error message if it is invalid.
%
%   CLASSES is either a space separated string or a cell array of strings
%   containing the set of classes that A is expected to belong to.  For
%   example, CLASSES could be 'uint8 double' if A can be either uint8 or
%   double.  CLASSES could be {'logical' 'cell'} if A can be either logical
%   or cell.  The string 'numeric' is interpreted as an abbreviation for all
%   the numeric classes.
%
%   ATTRIBUTES is either a space separated string or a cell array of strings
%   containing the set of attributes that A must satisfy.  To see the list of
%   valid attributes, see the subfunction init_table below.  For example, if
%   ATTRIBUTES is {'real' 'nonempty' 'finite'}, then A must be real and
%   nonempty, and it must contain only finite values.
%
%   FUNCTION_NAME is a string containing the function name to be used in the
%   formatted error message.
%
%   VARIABLE_NAME is a string containing the documented variable name to be
%   used in the formatted error message.
%
%   ARGUMENT_POSITION is a positive integer indicating which input argument
%   is being checked; it is also used in the formatted error message.

%   Copyright 1993-2002 The MathWorks, Inc.
%   $Revision: 1.5 $  $Date: 2002/05/20 20:27:03 $

% Input arguments are not checked for validity.

check_classes(A, classes, function_name, variable_name, argument_position);

check_attributes(A, attributes, function_name, variable_name, ...
                 argument_position);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = is_numeric(A)

numeric_classes = {'double' 'uint8' 'uint16' 'uint32' 'int8' ...
                   'int16' 'int32' 'single'};

tf = false;
for p = 1:length(numeric_classes)
    if isa(A, numeric_classes{p})
        tf = true;
        break;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = expand_numeric(in)
% Converts the string 'numeric' to the equivalent cell array containing the
% names of the numeric classes.

out = in(:);

idx = strmatch('numeric', out, 'exact');
if (length(idx) == 1)
    out(idx) = [];
    numeric_classes = {'uint8', 'int8', 'uint16', 'int16', ...
                       'uint32', 'int32', 'single', 'double'}';
    out = [out; numeric_classes];
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_classes(A, classes, function_name, ...
                       variable_name, argument_position)

if isempty(classes)
    return
end

if ischar(classes)
    if isempty(classes)
        % Work around bug in strread.
        classes = {};
    else
        classes = strread(classes, '%s');
    end
end

is_valid_type = false;
for k = 1:length(classes)
    if strcmp(classes{k}, 'numeric') && is_numeric(A)
        is_valid_type = true;
        break;

    else        
        if isa(A, classes{k})
            is_valid_type = true;
            break;
        end
    end
end

if ~is_valid_type
    messageId = sprintf('Images:%s:%s', function_name, 'invalidType');
    classes = expand_numeric(classes);
    validTypes = '';
    for k = 1:length(classes)
        validTypes = [validTypes, classes{k}, ', '];
    end
    validTypes(end-1:end) = [];
    message1 = sprintf('Function %s expected its %s input argument, %s,', ...
                       function_name, ...
                       num2ordinal(argument_position), ...
                       variable_name);
    message2 = 'to be one of these types:';
    message3 = sprintf('Instead its type was %s.', class(A));
    error(messageId, '%s\n%s\n\n  %s\n\n%s', message1, message2, validTypes, ...
          message3);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_attributes(A, attributes, function_name, ...
                          variable_name, argument_position)

if ischar(attributes)
    if isempty(attributes)
        % Work around bug in strread.
        attributes = {};
    else
        attributes = strread(attributes, '%s');
    end
end

table = init_table;

for k = 1:length(attributes)
    if strcmp(attributes{k}, '2d')
        tableEntry = table.twod;
    else
        tableEntry = table.(attributes{k});
    end
    
    if ~feval(tableEntry.checkFunction, A)
        messageId = sprintf('Images:%s:%s', function_name, ...
                            tableEntry.mnemonic);
        message1 = sprintf('Function %s expected its %s input argument, %s,', ...
                           function_name, num2ordinal(argument_position), ...
                           variable_name);
        message2 = sprintf('to be %s.', tableEntry.endOfMessage);
        error(messageId, '%s\n%s', message1, message2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_real(A)

try
    tf = isreal(A);
catch
    tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_even(A)

try
    tf = ~any(rem(double(A(:)),2));
catch
    tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_vector(A)

try
    tf = (ndims(A) == 2) & (any(size(A) == 1) | all(size(A) == 0));
catch
    tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_row(A)

try
    tf = (ndims(A) == 2) & ((size(A,1) == 1) | isequal(size(A), [0 0]));
catch
    tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_column(A)

try
    tf = (ndims(A) == 2) & ((size(A,2) == 1) | isequal(size(A), [0 0]));
catch
    tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_scalar(A)

try
    tf = all(size(A) == 1);
catch
    tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_2d(A)

try
    tf = ndims(A) == 2;
catch
    tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_nonsparse(A)

try
    tf = ~issparse(A);
catch
    tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_nonempty(A)

try
    tf = ~isempty(A);
catch
    tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_integer(A)

try
    A = A(:);
    switch class(A)

      case {'double','single'}
        tf = all(floor(A) == A) & all(isfinite(A));

      case {'uint8','int8','uint16','int16','uint32','int32','logical'}
        tf = true;

      otherwise
        tf = false;
    end

catch
    tf = false;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_nonnegative(A)

try
    tf = all(A(:) >= 0);
catch
    tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_positive(A)

try
    tf = all(A(:) > 0);
catch
    tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_nonnan(A)

try
    tf = ~any(isnan(A(:)));
catch
    % if isnan isn't defined for the class of A,
    % then we'll end up here.  If isnan isn't
    % defined then we'll assume that A can't
    % contain NaNs.
    tf = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_finite(A)

try
    tf = all(isfinite(A(:)));
catch
    % if isfinite isn't defined for the class of A,
    % then we'll end up here.  If isfinite isn't
    % defined then we'll assume that A is finite.
    tf = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = check_nonzero(A)

try
    tf = ~all(A(:) == 0);
catch
    tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = init_table

persistent table

if isempty(table)
    table.real.checkFunction        = @check_real;
    table.real.mnemonic             = 'expectedReal';
    table.real.endOfMessage         = 'real';
    
    table.vector.checkFunction      = @check_vector;
    table.vector.mnemonic           = 'expectedVector';
    table.vector.endOfMessage       = 'a vector';
    
    table.row.checkFunction         = @check_row;
    table.row.mnemonic              = 'expectedRow';
    table.row.endOfMessage          = 'a row vector';
    
    table.column.checkFunction      = @check_column;
    table.column.mnemonic           = 'expectedColumn';
    table.column.endOfMessage       = 'a column vector';
    
    table.scalar.checkFunction      = @check_scalar;
    table.scalar.mnemonic           = 'expectedScalar';
    table.scalar.endOfMessage       = 'a scalar';
    
    table.twod.checkFunction        = @check_2d;
    table.twod.mnemonic             = 'expected2D';
    table.twod.endOfMessage         = 'two-dimensional';
    
    table.nonsparse.checkFunction   = @check_nonsparse;
    table.nonsparse.mnemonic        = 'expectedNonsparse';
    table.nonsparse.endOfMessage    = 'nonsparse';
    
    table.nonempty.checkFunction    = @check_nonempty;
    table.nonempty.mnemonic         = 'expectedNonempty';
    table.nonempty.endOfMessage     = 'nonempty';
    
    table.integer.checkFunction     = @check_integer;
    table.integer.mnemonic          = 'expectedInteger';
    table.integer.endOfMessage      = 'integer-valued';
        
    table.nonnegative.checkFunction = @check_nonnegative;
    table.nonnegative.mnemonic      = 'expectedNonnegative';
    table.nonnegative.endOfMessage  = 'nonnegative';
    
    table.positive.checkFunction    = @check_positive;
    table.positive.mnemonic         = 'expectedPositive';
    table.positive.endOfMessage     = 'positive';
    
    table.nonnan.checkFunction      = @check_nonnan;
    table.nonnan.mnemonic           = 'expectedNonNaN';
    table.nonnan.endOfMessage       = 'non-NaN';
    
    table.finite.checkFunction      = @check_finite;
    table.finite.mnemonic           = 'expectedFinite';
    table.finite.endOfMessage       = 'finite';
    
    table.nonzero.checkFunction     = @check_nonzero;
    table.nonzero.mnemonic          = 'expectedNonZero';
    table.nonzero.endOfMessage      = 'non-zero';
    
    table.even.checkFunction        = @check_even;
    table.even.mnemonic             = 'expectedEven';
    table.even.endOfMessage         = 'even';
    
end

out = table;

%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function checknargin(low, high, numInputs, function_name)
%CHECKNARGIN Check number of input arguments.
%   CHECKNARGIN(LOW,HIGH,NUM_INPUTS,FUNCTION_NAME) checks whether NUM_INPUTS
%   is in the range indicated by LOW and HIGH.  If not, CHECKNARGIN issues a
%   formatted error message using the string in FUNCTION_NAME.
%
%   LOW should be a scalar nonnegative integer.
%
%   HIGH should be a scalar nonnegative integer or Inf.
%
%   FUNCTION_NAME should be a string.

%   Copyright 1993-2002 The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 2002/03/15 15:57:05 $

% Input arguments are not checked for validity.

if numInputs < low
  msgId = sprintf('Images:%s:tooFewInputs', function_name);
  if low == 1
    msg1 = sprintf('Function %s expected at least 1 input argument', ...
                   function_name);
  else
    msg1 = sprintf('Function %s expected at least %d input arguments', ...
                   function_name, low);
  end
  
  if numInputs == 1
    msg2 = 'but was called instead with 1 input argument.';
  else
    msg2 = sprintf('but was called instead with %d input arguments.', ...
                   numInputs);
  end
  
  error(msgId, '%s\n%s', msg1, msg2);
  
elseif numInputs > high
  msgId = sprintf('Images:%s:tooManyInputs', function_name);

  if high == 1
    msg1 = sprintf('Function %s expected at most 1 input argument', ...
                   function_name);
  else
    msg1 = sprintf('Function %s expected at most %d input arguments', ...
                   function_name, high);
  end
  
  if numInputs == 1
    msg2 = 'but was called instead with 1 input argument.';
  else
    msg2 = sprintf('but was called instead with %d input arguments.', ...
                   numInputs);
  end
  
  error(msgId, '%s\n%s', msg1, msg2);
end
    
%**************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function out = checkstrs(in, valid_strings, function_name, ...
                         variable_name, argument_position)
%CHECKSTRS Check validity of option string.
%   OUT = CHECKSTRS(IN,VALID_STRINGS,FUNCTION_NAME,VARIABLE_NAME, ...
%   ARGUMENT_POSITION) checks the validity of the option string IN.  It
%   returns the matching string in VALID_STRINGS in OUT.  CHECKSTRS looks
%   for a case-insensitive nonambiguous match between IN and the strings
%   in VALID_STRINGS.
%
%   VALID_STRINGS is a cell array containing strings.
%
%   FUNCTION_NAME is a string containing the function name to be used in the
%   formatted error message.
%
%   VARIABLE_NAME is a string containing the documented variable name to be
%   used in the formatted error message.
%
%   ARGUMENT_POSITION is a positive integer indicating which input argument
%   is being checked; it is also used in the formatted error message.

%   Copyright 1993-2002 The MathWorks, Inc.
%   $Revision: 1.3 $  $Date: 2002/03/15 15:57:05 $

% Except for IN, input arguments are not checked for validity.

checkinput(in, 'char', 'row', function_name, variable_name, argument_position);

idx = strmatch(lower(in), valid_strings);

num_matches = prod(size(idx));

if num_matches == 1
    out = valid_strings{idx};

else
    % Convert valid_strings to a single string containing a space-separated list
    % of valid strings.
    list = '';
    for k = 1:length(valid_strings)
        list = [list ', ' valid_strings{k}];
    end
    list(1:2) = [];

    msg1 = sprintf('Function %s expected its %s input argument, %s,', ...
                   function_name, num2ordinal(argument_position), ...
                   variable_name);
    msg2 = 'to match one of these strings:';

    if num_matches == 0
        msg3 = sprintf('The input, ''%s'', did not match any of the valid strings.', in);
        id = sprintf('Images:%s:unrecognizedStringChoice', function_name);

    else
        msg3 = sprintf('The input, ''%s'', matched more than one valid string.', in);
        id = sprintf('Images:%s:ambiguousStringChoice', function_name);
    end

    error(id,'%s\n%s\n\n  %s\n\n%s', msg1, msg2, list, msg3);
end

