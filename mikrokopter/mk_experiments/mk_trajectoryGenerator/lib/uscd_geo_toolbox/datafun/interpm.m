function vi = interpm(varargin)

% INTERPM N-D interpolation (table lookup).
%   VI = INTERPM(X1,X2,X3,...,V,Y1,Y2,Y3,...) interpolates to find VI,
%   the values of the underlying N-D function V at the points in
%   arrays Y1,Y2,Y3,etc.  For an N-D V, INTERPM should be called with
%   2*N+1 arguments.  Arrays X1,X2,X3,etc. specify the points at which
%   the data V is given.  Out of range values are returned as NaN's.
%   Y1,Y2,Y3,etc. must be arrays of the same size or vectors.  Vector
%   arguments that are not the same size are passed through NDGRID to
%   create the Y1,Y2,Y3,etc. arrays.  INTERPM works for all N-D arrays
%   with 2 or more dimensions.
%
%   VI = INTERPM(V,Y1,Y2,Y3,...) assumes X1=1:SIZE(V,1),X2=1:SIZE(V,2),etc.
%   VI = INTERPM(V,NTIMES) expands V by interleaving interpolates between
%   every element, working recursively for NTIMES.  
%   VI = INTERPM(V) is the same as INTERPM(V,1).
%
%   VI = INTERPM(...,'method') specifies alternate methods.  The default
%   is linear interpolation.  Available methods are:
%
%     'linear'  - linear interpolation
%     'cubic'   - cubic interpolation
%     'nearest' - nearest neighbor interpolation
%     'spline'  - spline interpolation
%   
%   INTERPM requires that X1,X2,X3,etc. be monotonic and plaid (as if
%   they were created using NDGRID).  X1,X2,X3,etc. can be non-uniformly
%   spaced. For faster interpolation when X1,X2,Y3,etc. are known to be
%   equally spaced and monotonic, use the methods '*linear', '*cubic', or
%   '*nearest'.
%
%   See also INTERP1, INTERP2, INTERP3, NDGRID.

%   Clay M. Thompson 8-1-94
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 1.26 $  $Date: 1997/11/21 23:40:50 $

bypass = 0;
uniform = 1;
if isstr(varargin{end}),
  narg = nargin-1;
  method = [varargin{end} '    ']; % Protect against short string.
  if method(1)=='*', % Direct call bypass.
    if method(2)=='l', % linear interpolation.
      vi = linear(varargin{1:end-1});
      return
    elseif method(2)=='c', % cubic interpolation
      vi = cubic(varargin{1:end-1});
      return
    elseif method(2)=='n', % Nearest neighbor interpolation
      vi = nearest(varargin{1:end-1});
      return
    elseif method(2)=='s', % spline interpolation
      method = 'spline'; bypass = 1;
    else
      error([deblank(method),' is an invalid method.']);
    end
  elseif method(1)=='s', % Spline interpolation
    method = 'spline'; bypass = 1;
  end
else
  narg = nargin;
  method = 'linear';
end

if narg==0,
  error('Not enough input arguments.');

elseif narg<=2. % interp(v) or interpm(v,n), Expand V ntimes
  if narg==1, ntimes = 1; else ntimes = floor(varargin{2}); end
  siz = size(varargin{1});
  x = cell(size(siz)); y = cell(size(siz));
  for i=1:length(siz),
    x{i} = 1:siz(i);
    y{i} = 1:1/(2^ntimes):siz(i);
  end
  [msg,x,v,y] = xnchk(x,varargin{1},y);

elseif all(size(varargin{1}) >= 2) & (narg == ndims(varargin{1})+1),
   % If every dimension of V has size >=2 and nargin = ndims + 1,
   % then handle interpm(v,y1,y2,y3,...).
  v = varargin{1};
  if ndims(v)~=narg-1, error('Wrong number of input arguments.'); end
  siz = size(v);
  y = varargin(2:narg);
  x = cell(size(siz));
  for i=1:length(siz), x{i} = 1:siz(i); end
  [msg,x,v,y] = xnchk(x,v,y);

elseif (rem(narg,2)==1) & all(size(varargin{floor((narg+1)/2)}) >= 2) & ... 
       (narg == 2*ndims(varargin{(narg+1)/2})+1),
  % If the number of input arguments is odd and every dimension of V
  % has size >= 2 and nargin = 2*ndims+1,
  % then handle interpm(x1,x2,x3,...,v,y1,y2,y3,...)
  v = varargin{(narg+1)/2};
  siz = size(v);
  x = varargin(1:ndims(v));
  y = varargin(ndims(v)+2:narg);
  [msg,x,v,y] = xnchk(x,v,y);
else
  error(...
    'Wrong number of input arguments or some dimension of V is less than 2.');
end

if ~isempty(msg), error(msg); end


is_flip = zeros(1,length(x));

if ~bypass,

  % Create xx cell array containing the vectors from each Xi array
  xx = cell(size(x));
  for i=1:length(x),
    ind(1:length(x)) = {1}; 
    ind{i} = ':';
    xx{i} = x{i}(ind{:});
    xx{i} = xx{i}(:); % Make sure its a column
  end
  
  %
  % Check for non-equally spaced data.  If so, map(x1,x2,x3,...) and
  % (y1,y2,y3,...) to array index coordinate system.
  %
  for i=1:prod(size(xx));
    dx{i} = diff(xx{i});
    xxmax(i) = max(abs(xx{i}));
    if prod(size(dx{i}))==1,
      maxabs(i) = 0;
    else 
      maxabs(i) = max(abs(diff(dx{i})));
    end
  end
  if any(maxabs > eps*xxmax), % If data is not equally spaced,
    for i=1:length(x),
      % Flip orientation of data so that x{i} is increasing
      is_flip(i) = any( dx{i} < 0 );
      if is_flip(i)
        for j=1:length(x), x{j} = flipdim(x{j},i); end
        for j=1:length(y), y{j} = flipdim(y{j},i); end
        v = flipdim(v,i);
        xx{i} =  flipud(xx{i});
        dx{i} = -flipud(dx{i});
      end
    end
  
    for i=1:length(dx),
      if any(dx{i}<=0), 
        error('X1,X2,X3, etc. must be monotonic vectors or arrays produced by NDGRID.'); 
      end
    end
  
    % Bypass mapping code for cubic
    if method(1)~='c', 
      % Map values in y to values in si via linear interpolation (one for
      % each input argument.
      for k=1:length(y),
        zz = xx{k}; zi = y{k}; dz = dx{k}; si = [];
        
        % Determine the nearest location of zi in zz
        [zzi,j] = sort(zi(:));
        [dum,i] = sort([zz;zzi]);
        si(i) = (1:length(i));
        si = (si(length(zz)+1:end)-(1:length(zzi)))';
        si(j) = si;
        
        % Map values in zi to index offset (si) tia linear interpolation
        si(si<1) = 1;
        si(si>length(zz)-1) = length(zz)-1;
        si = si + (zi(:)-zz(si))./(zz(si+1)-zz(si));
        
        y{k}(:) = si; % Replace value in y
      end
      
      % Change x1,x2,x3,... to be index cordinates.
      for i=1:length(x),
        x{i} = 1:siz(i);
      end
      [x{1:end}] = ndgrid(x{:});
    else
      uniform = 0;
    end
  end
end

% Now do the interpolation based on method.
method = [lower(method),'   ']; % Protect against short string

if method(1)=='l', % linear interpolation.
  vi = linear(x{:},v,y{:});

elseif method(1)=='c', % cubic interpolation
  if uniform
    vi = cubic(x{:},v,y{:});
  else
    d = zeros(size(y{1}));
    for i=1:length(y)
      d = d | y{i} < min(x{i}(:)) | y{i} > max(x{i}(:));
    end
    d = find(d);
    vi = splinen(x,v,y);
    vi(d) = NaN;
  end

elseif method(1)=='n', % Nearest neighbor interpolation
  vi = nearest(x{:},v,y{:});

elseif method(1)=='s', % Spline interpolation
  vi = splinen(x,v,y);

else
  error([deblank(method),' is an invalid method.']);

end


% Flip orientation of data back

for i = find(is_flip)
  vi = flipdim(vi,i);
end

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function F = linear(varargin)

% LINEAR N-D linear data interpolation.
%   VI = LINEAR(X1,X2,X3,...,V,Y1,Y2,Y3,...) interpolates to find VI,
%   the values of the underlying N-D function V, at the points in arrays
%   Y1,Y2,Y3, etc. via linear interpolation.  For an N-D V, LINEAR
%   should be called with 2*N+1 arguments.  Arrays X1,X2,X3,... specify
%   the points at which the data V is given.  X1,X2,X3, etc. can also be
%   vectors specifying the abscissae for the matrix V as for NDGRID. In
%   both cases, X1,X2,X3, etc. must be equally spaced and monotonic.
%   Out of range values are returned as NaN.
%
%   VI = LINEAR(V,Y1,Y2,Y3,...) assumes X1=1:SIZE(V,1), X2=1:SIZE(V,2),
%   X3=1:SIZE(V,3), etc.
%   VI = LINEAR(V,NTIMES) expands V by interleaving interpolates between
%   every element, working recursively for NTIMES.  LINEAR(V) is the
%   same as LINEAR(V,1).
%
%   See also INTERPM.

%   Clay M. Thompson 8-2-94

if nargin==0,
  error('Not enough input arguments.');

elseif nargin<=2, % linear(v) or linear(v,n), Expand V n times
  if nargin==1, ntimes = 1; else ntimes = floor(varargin{2}); end
  v = varargin{1};
  siz = size(v);
  s = cell(size(siz));
  for i=1:length(s),
    s{i} = 1:1/(2^ntimes):siz(i);
  end
  [s{1:end}] = ndgrid(s{:});

elseif nargin>2 & (rem(nargin,2)==0 | nargin == ndims(varargin{1})+1),
  % linear(v,y1,y2,y3,...)
  v = varargin{1};
  if ndims(v)~=nargin-1, 
    error('Wrong number of input arguments.');
  end
  siz = size(v);
  s = varargin(2:nargin);

elseif nargin>2 & rem(nargin,2)==1 &  ...
       (nargin == 2*ndims(varargin{(nargin+1)/2})+1),
  % linear(x1,x2,x3,...,v,y1,y2,y3,...)
  v = varargin{(nargin+1)/2};
  if nargin~=2*ndims(v)+1, error('Wrong number of input arguments.'); end
  siz = size(v);
  x = varargin(1:ndims(v));
  y = varargin(ndims(v)+2:nargin);
  [msg,x,v,y] = xnchk(x,v,y);
  for i=length(x):-1:1,
    xsiz{i} = size(x{i});
    xprodsiz(i) = prod(size(x{i}));
  end
  if ~isequal(xprodsiz,siz) & ~isequal(siz,xsiz{:}),
    error('The lengths of the X1,X2,X3,... vectors must match V.');
  end
  s = cell(size(y));
  for i=1:length(s),
    s{i} = 1 + (y{i}-x{i}(1))/(x{i}(xprodsiz(i))-x{i}(1))*(siz(i)-1);
  end
end

if any(siz<2*ones(size(siz))),
  error('V must be at least 2-by-2-by-2-by-...');
end

for i=length(s):-1:1, ssiz{i} = size(s{i}); end
if ~isequal(ssiz{:}),
  error('Y1,Y2,Y3, etc. must be the same size.');
end

% Check for out of range values of s and set to 1
sout = cell(size(s));
for i=1:length(s),
  sout{i} = find((s{i}<1) | (s{i}>siz(i)));
  if ~isempty(sout{i}), s{i}(sout{i}) = ones(size(sout{i})); end
end

% Matrix element indexing
offset = cumprod([1 siz(1:end-1)]);
ndx = 1;
for i=1:length(s),
  ndx = ndx + offset(i)*floor(s{i}-1);
end

% Compute intepolation parameters, check for boundary value.
for i=1:length(s),
  if isempty(s{i}), d = s{i}; else d = find(s{i}==siz(i)); end
  s{i} = s{i}-floor(s{i});
  if ~isempty(d), s{i}(d) = s{i}(d)+1; ndx(d) = ndx(d)-offset(i); end
end
d = []; % Reclaim memory

% Create index arrays, iw.
iw = cell(size(s));
[iw{1:end}] = ndgrid(0:1);

% Reshape each iw{i} to a column and then arrange into a matrix
iwcol = [prod(size(iw{1})) 1];
for i=1:prod(size(iw)), iw{i} = reshape(iw{i},iwcol); end 
iw = cat(2,iw{:}); % Arrange columns into a matrix

% Do the linear interpolation: f = v1*(1-s) + v2*s along each direction
F = 0;
for i=1:size(iw,1),
  vv = v(ndx + offset*iw(i,:)');
  for j=1:size(iw,2),
    switch iw(i,j)
    case 0 % Interpolation function (1-s)
      vv = vv.*(1-s{j});
    case 1 % Interpolation function (s)
      vv = vv.*s{j};
    end
  end
  F = F + vv;
end

% Now set out of range values to NaN.
for i=1:length(sout),
  if ~isempty(sout{i}), F(sout{i}) = NaN; end
end


%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function F = cubic(varargin)

% CUBIC Cubic data interpolation.
%   CUBIC(...) is the same as LINEAR(....) except that it uses
%   cubic interpolation.
%   
%   See also INTERPM.

%   Clay M. Thompson 8-1-94
%   Based on "Cubic Convolution Interpolation for Digital Image
%   Processing", Robert G. Keys, IEEE Trans. on Acoustics, Speech, and
%   Signal Processing, Vol. 29, No. 6, Dec. 1981, pp. 1153-1160.

if nargin==0,
  error('Not enough input arguments.');

elseif nargin<=2, % cubic(v) or cubic(v,n), Expand V n times
  if nargin==1, ntimes = 1; else ntimes = floor(varargin{2}); end
  v = varargin{1};
  siz = size(v);
  s = cell(size(siz));
  for i=1:length(s),
    s{i} = 1:1/(2^ntimes):siz(i);
  end
  [s{1:end}] = ndgrid(s{:});

elseif nargin>2 & (rem(nargin,2)==0 | nargin == ndims(varargin{1})+1),
  % cubic(v,y1,y2,y3,...)
  v = varargin{1};
  if ndims(v)~=nargin-1, 
    error('Wrong number of input arguments.');
  end
  siz = size(v);
  s = varargin(2:nargin);

elseif nargin>2 & rem(nargin,2)==1 &  ...
       (nargin == 2*ndims(varargin{(nargin+1)/2})+1),
  % cubic(x1,x2,x3,...,v,y1,y2,y3,...)
  v = varargin{(nargin+1)/2};
  if nargin~=2*ndims(v)+1, error('Wrong number of input arguments.'); end
  siz = size(v);
  x = varargin(1:ndims(v));
  y = varargin(ndims(v)+2:nargin);
  [msg,x,v,y] = xnchk(x,v,y);
  for i=length(x):-1:1,
    xsiz{i} = size(x{i});
    xprodsiz(i) = prod(size(x{i}));
  end
  if ~isequal(xprodsiz,siz) & ~isequal(siz,xsiz{:}),
    error('The lengths of the X1,X2,X3,... vectors must match V.');
  end
  s = cell(size(y));
  for i=1:length(s),
    s{i} = 1 + (y{i}-x{i}(1))/(x{i}(xprodsiz(i))-x{i}(1))*(siz(i)-1);
  end
end

if any(siz<3*ones(size(siz))),
  error('V must be at least 3-by-3-by-3-by-...');
end

for i=length(s):-1:1, ssiz{i} = size(s{i}); end
if ~isequal(ssiz{:}),
  error('Y1,Y2,Y3, etc. must be the same size.');
end

% Check for out of range values of s and set to 1
sout = cell(size(s));
for i=1:length(s),
  sout{i} = find((s{i}<1) | (s{i}>siz(i)));
  if ~isempty(sout{i}), s{i}(sout{i}) = ones(size(sout{i})); end
end

% Matrix element indexing
offset = cumprod([1 siz(1:end-1)+2]);
ndx = 1;
for i=1:length(s),
  ndx = ndx + offset(i)*floor(s{i}-1);
end

% Compute intepolation parameters, check for boundary value.
for i=1:length(s),
  if isempty(s{i}), d = s{i}; else d = find(s{i}==siz(i)); end
  s{i} = s{i}-floor(s{i});
  if ~isempty(d), s{i}(d) = s{i}(d)+1; ndx(d) = ndx(d)-offset(i); end
end
d = []; % Reclaim memory


% Expand v so interpolation is valid at the boundaries.  

vv = zeros(siz+2);
for i=length(siz):-1:1,
  ind{i} = 2:siz(i)+1;
end
vv(ind{:}) = v;

% Insert values on the boundary that allow the use of the same
% interpolation kernel everywhere.  These values are computed
% computed from the 3 strips nearest the boundary:
%   v(boundary,:) = 3*v(1,:) - 3*v(2,:) + v(3,:);
% Use comma separated list equivalence to index into vv.  That is,
% vv(edges{1,:}) is the same as vv(edges{1,1},edges{1,2},edges{1,3},...).

for i=length(siz):-1:1,
  ind{i} = 1:siz(i)+2; % vv(ind{:}) is the same as vv(:,:,...,:).
end

for i=1:length(s),
  edges = ind(ones(4,1),:); % Reinitialize edges

  % Set virtual right edge
  edges(:,i) = {1;2;3;4};
  vv(edges{1,:}) = 3*vv(edges{2,:}) - 3*vv(edges{3,:}) + vv(edges{4,:});

  % Set virtual left edge
  edges(:,i) = num2cell(siz(i)+[2 1 0 -1]');
  vv(edges{1,:}) = 3*vv(edges{2,:}) - 3*vv(edges{3,:}) + vv(edges{4,:});
end
siz = siz + 2;

% Create index arrays, iw.
iw = cell(size(s));
[iw{1:end}] = ndgrid(0:3);

% Reshape each iw{i} to a column and then arrange into a matrix
iwcol = [prod(size(iw{1})) 1];
for i=1:prod(size(iw)), iw{i} = reshape(iw{i},iwcol); end 
iw = cat(2,iw{:});

% Do the cubic interpolation: 
% f = v1*((2-s)*s-1)*s     + v2*((3*s-5)*s^2+2) + ...
%     v3*(((4-3*s)*s+1)*s) + v4*((s-1)*s^2) along each direction
F = 0;
for i=1:size(iw,1),
  vvv = vv(ndx + offset*iw(i,:)');
  for j=1:size(iw,2),
    switch iw(i,j)
    case 0 % Interpolation function ((2-s)*s-1)*s
      vvv = vvv.*( ((2-s{j}).*s{j}-1).*s{j} );
    case 1 % Interpolation function (3*s-5)*s*s+2
      vvv = vvv.*( (3*s{j}-5).*s{j}.*s{j}+2 );
    case 2 % Interpolation function ((4-3*s)*s+1)*s
      vvv = vvv.*( ((4-3*s{j}).*s{j}+1).*s{j} );
    case 3 % Interpolation function (s-1)*s*s
      vvv = vvv.*( (s{j}-1).*s{j}.*s{j} );
    end
  end
  F = F + vvv;
end
F(:) = F/pow2(length(s));

% Now set out of range values to NaN.
for i=1:length(sout),
  if ~isempty(sout{i}), F(sout{i}) = NaN; end
end


%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function F = nearest(varargin)

% NEAREST 3-D Nearest neighbor interpolation.
%   NEAREST(...) is the same as LINEAR(...) except that it uses
%   nearest neighbor interpolation.
%
%   See also INTERPM.

%   Clay M. Thompson 1-31-94

if nargin==0,
  error('Not enough input arguments.');

elseif nargin<=2, % nearest(v) or nearest(v,n), Expand V n times
  if nargin==1, ntimes = 1; else ntimes = floor(varargin{2}); end
  v = varargin{1};
  siz = size(v);
  s = cell(size(siz));
  for i=1:length(s),
    s{i} = 1:1/(2^ntimes):siz(i);
  end
  [s{1:end}] = ndgrid(s{:});

elseif nargin>2 & (rem(nargin,2)==0 | nargin == ndims(varargin{1})+1),
  % nearest(v,y1,y2,y3,...)
  v = varargin{1};
  if ndims(v)~=nargin-1, 
    error('Wrong number of input arguments.');
  end
  siz = size(v);
  s = varargin(2:nargin);

elseif nargin>2 & rem(nargin,2)==1 &  ...
       (nargin == 2*ndims(varargin{(nargin+1)/2})+1),
  % nearest(x1,x2,x3,...,v,y1,y2,y3,...)
  v = varargin{(nargin+1)/2};
  if nargin~=2*ndims(v)+1, error('Wrong number of input arguments.'); end
  siz = size(v);
  x = varargin(1:ndims(v));
  y = varargin(ndims(v)+2:nargin);
  [msg,x,v,y] = xnchk(x,v,y);
  for i=length(x):-1:1,
    xsiz{i} = size(x{i});
    xprodsiz(i) = prod(size(x{i}));
  end
  if ~isequal(xprodsiz,siz) & ~isequal(siz,xsiz{:}),
    error('The lengths of the X1,X2,X3,... vectors must match V.');
  end
  s = cell(size(y));
  for i=1:length(s),
    s{i} = 1 + (y{i}-x{i}(1))/(x{i}(xprodsiz(i))-x{i}(1))*(siz(i)-1);
  end
end

for i=length(s):-1:1, ssiz{i} = size(s{i}); end
if ~isequal(ssiz{:}),
  error('Y1,Y2,Y3, etc. must be the same size.');
end

% Check for out of range values of s and set to 1
sout = cell(size(s));
for i=1:length(s),
  sout{i} = find((s{i}<.5) | (s{i}>siz(i)+.5));
  if ~isempty(sout{i}), s{i}(sout{i}) = ones(size(sout{i})); end
end

% Matrix element indexing
offset = cumprod([1 siz(1:end-1)]);
ndx = 1;
for i=1:length(s),
  ndx = ndx + offset(i)*round(s{i}-1);
end

% Now interpolate
F = v(ndx);

% Now set out of range values to NaN.
for i=1:length(sout),
  if ~isempty(sout{i}), F(sout{i}) = NaN; end
end

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,x,v,y] = xnchk(x,v,y)

% XNCHK Check arguments to N-D data routines.
%   [MSG,X,V,Y] = XNCHK(X,V,Y), checks the input aguments and
%   returns either an error message in MSG or valid X,V, and Y
%   data.  X and Y are cell array lists that contain vectors or
%   matrices.


if nargin~=3, error('Wrong number of input arguments.'); end
msg = [];

% Check to make sure the number of dimensions of V matches the
% number of x and y arguments.
if length(x) ==length(y) & ndims(v) ~= length(x),
  msg = sprintf('V must be a %d-D array.',length(x));
  return
end

siz = size(v);
isvec = zeros(size(x));
for i=length(x):-1:1,
  xsiz{i} = size(x{i});
  isvec(i) = isvector(x{i});
end

if ~isvector(v), % v is not a vector or scalar
  % Convert x,y,z to row, column, and page matrices if necessary.
  if all(isvec),
    [x{1:end}] = ndgrid(x{:}); % Grid the vectors
    for i=length(x):-1:1, xsiz{i} = size(x{i}); end
    if ~isequal(siz,xsiz{:}),
      msg = 'The lengths of X1,X2,X3,etc. must match the size of V.';
      return
    end
  elseif any(isvec),
    msg = 'X1,X2,X3, etc. must all be vectors or all be arrays.';
    return
  else
    if ~isequal(siz,xsiz{:}),
      msg = 'Arrays X1,X2,X3, etc. must be the same size as V.';
      return
    end
  end
elseif isvector(v) % v is a vector
  for i=length(x):-1:1, xlength(i) = length(x{i}); end
  if any(~isvec),
    msg = 'X1,X2,X3, etc. must be vectors when V is.';
    return
  elseif ~isequal(length(v),xlength{:}),
    msg = 'X1,X2,X3, etc. must be the same length as V.';
    return
  end
end

isvec = zeros(size(y));
for i=length(y):-1:1, 
  ysiz{i} = size(y{i}); 
  isvec(i) =isvector(y{i});
end

% If y1,y2,y3, etc. don't all have the same orientation, then
% build y1,y2,y3, etc. arrays.
if ismesh(y{:})
  [y{1:end}] = ndgrid(y{:});
elseif ~isequal(ysiz{:}),
  msg = sprintf('Y1,Y2,Y3, etc. must all the same size or vectors of different orientations.');
end

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function tf = isvector(x)

% ISVECTOR True if x has only one non-singleton dimension.

tf = length(x)==prod(size(x));

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function F = splinen(x,v,xi)

% N-D spline interpolation

% Determine abscissa vectors 
for i=1:length(x);
  ind(1:length(x)) = {1}; ind{i} = ':';
  x{i} = reshape(x{i}(ind{:}),1,size(x{i},i));
end

%
% Check for plaid data.
%
isplaid = 1;
for i=1:length(xi),
  ind(1:length(xi)) = {1}; ind{i} = ':';
  siz = size(xi{i}); siz(i) = 1;
  xxi{i} = xi{i}(ind{:});
  if ~isequal(repmat(xxi{i},siz),xi{i})
     isplaid = 0;
  end
  xxi{i} = xxi{i}(:).'; % Make it a row
end

if ~isplaid
  F = splncore(x,v,xi);
else
  F = splncore(x,v,xxi,'gridded');
end

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function tf = ismesh(varargin)

%ISMESH True if the inputs should be automatically meshgridded.
%    ISMESH(X,Y) returns true if X and Y are vectors of
%    different orientations.
%
%    ISMESH(X,Y,Z) returns true if X,Y,Z are vectors of
%    different orientations.
%
%    ISMESH(...) returns true if all the inputs are vectors of
%    different orientations.

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%    $Revision: 1.4 $ $Date: 1997/11/21 23:40:30 $

for i=1:length(varargin)
  ns{i} = size(varargin{i})~=1; % Location of non-singleton dimensions
  isvec(i) = sum(ns{i})<=1;     % Is vector.
  nd(i) = ndims(varargin{i});    % Number of dimensions.
end

% True if inputs are 2-D, all vectors, and their non-singleton
% dimensions aren't along the same dimension.
tf = all(isvec) & ~isequal(ns{:});
