function [xi,yi,zi,em] = objmap(xb,yb,zb,xi,yi,fun,p,err,n,arg10)
% OBJMAP Objective mapping interpolation.
%	[ZI,EM] = OBJMAP( ... ,[LX LY],E) Specifies
%	lengthscales LX, LY and relative error E
%	(0 < E < 1) for "classical" gaussian correlation 
%	function 
%	C(x,y) = E*D(x,y)+(1-E)*exp(-(x/LX)^2-(y/LY)^2),
%	where D - Dirac delta function.
%       The LengthScale could be a 3-element Vector:
%        [LX LY Phi], in this case you get a rotated Elipse:
%        C(x,y) = E*D(x,y)+(1-E)* ...
%         exp( -((x*cos(phi)+y*sin(phi))/LX)^2 ...
%              -((y*cos(phi)-x*sin(phi))/LY)^2 )
%
%       OBJMAP(...,'notrend') didn't remove Trend, only Mean.
%
%	OBJMAP(...,'FUN',P,E)  Allows to specify the
%	correlation function 'FUN' as a string (expression
%	or function name) depending on x, y, r(r^2=x^2+y^2)
%	and parameters (such as lengthscales) in the vector
%	P.
%
%	For large number of known points it divides
%	the domain into subdomains using the adaptive
%	QUADTREE procedure.
%
%	OBJMAP(...,[LX LY],E,OPT) or 
%	OBJMAP(...,'FUN',P,E,OPT) also allows OPT vector
%	argument to specify several parameters of quadtree
%	division:
%	OPT = [NB ND NMAX PB VERBOSE], where
%	NB   - max. number of points in one block
%	ND   - "depopulation" threshold - if number of
%	  points in a block is less than ND is is
%	  considered "depopulated" and its own "secondary"
%	  neighbours are used for interpolation.
%	NMAX - max. number of points below which domain
%	  is not divided into blocks.
%	PB - X,Y scales as part of average block sizes,
%	VERBOSE - verbosity (1 - display some values and
%	  number of processed blocks.
%	Default values [NB ND NMAX PB V] = [32 8 500 1/3 1].
%
%	For more information about quadtree division and
%	options see TREEINFO and TREEDEMO.
%
%	XI can be a row vector, in which case it specifies
%	a matrix with constant columns. Similarly, YI can
%	be a column vector and it specifies a matrix with 
%	constant rows.
%	[XI,YI,ZI] = OBJMAP(X,Y,Z,XI,YI, ...) also returns 
%	matrices XI, YI formed from input vectors XI,YI 
%	in the way described above.
%       When the Rows of Z have the same number like Elements 
%        in X,Y , the Columns of Z will be formed to 
%        the 3. Dimension of ZI (only in Matlab5.#). 
%
%	See also GRIDDATA, MINCURVI, KRIGING, QUADTREE.

%  Copyright (c) 1995 by Kirill K. Pankratov
%	kirill@plume.mit.edu
%	05/25/95, 05/31/95

 % Defaults and parameters .............................
n_dflt = [32 8];  % Default for max. and min number of 
                  %  points in elementary block
nmax = 500;       % Default for max number of points
                  %  to be treated as one block
part_blk = 1/3;   % Default for x and y lengthscales
                  %  as parts of average block sizes
                  % when LY, LY not specified.
verbose = 1;      % Verbosity (shows some parameters and 
                  % number of blocks processed).


 % Handle input ........................................
if nargin==0, help objmap, return, end
if nargin==1
 if strcmp(xb,'info'), more on, type mapinfo, more off, return, end
end
if nargin<5
  error('  Not enough input arguments')
end

no_trend = 0;                                                                  
 if nargin > 5, no_trend = [ strcmp(lower(fun),'notrend') & isstr(fun) ]; end
 if nargin > 6, no_trend = [ strcmp(lower(p),'notrend') & isstr(p) ]; end      
 if nargin > 7, no_trend = [ strcmp(lower(err),'notrend') & isstr(err) ]; end  
 if nargin > 8, no_trend = [ strcmp(lower(n),'notrend') & isstr(n) ]; end      
 if nargin > 9, no_trend = [ strcmp(lower(arg10),'notrend') & isstr(arg10) ];
  end

 % Insert defaults for missing arguments
if nargin<9, n = 0; end
if nargin<8, err = 0; end
if nargin<7, p = 0; end
if nargin<6, fun = 0; end
is_call = [ isstr(fun) & ~strcmp(fun,'notrend') ];
  % If function name or expression

if is_call   % Function name or expression
  call = callchk(fun,'r');
else         % Lengthscales for gaussian form
  if ~isstr(err), n = err; else, n=0; end
  if ~isstr(p) , err = p; else, err = 0; end
  if ~isstr(fun) , sc = fun; else, sc = []; end
  call = '';
end


if ~length(sc), sc = [0 0]; end
if length(sc)==1, sc = sc([1 1]); end
if length(sc)>=2,sc(1:2) = max(sc(1:2),0); end 

 % Options (vector n)
if n==0 | isempty(n), n = n_dflt; end
len = length(n);
if len>=3, nmax = n(3); end
if len>=4, part_blk = n(4); end
if len>=5, verbose = n(5); end

 % Relative error:
if isempty(err), err=0; end  % Not empty
err = min(err,1);       % Not larger than 1

 % If error map is needed
is_err = (nargout==2) | (nargout==4);

 % Check input coordinates and make matrices XI, YI
 % if necessary
if prod(size(xb)) == prod(size(yb)) & ...
   prod(size(xb)) ==  size(zb,1)
 [msg,xb,yb,dd,xi,yi] = xyzchk(xb,yb,zb(:,1),xi,yi); 
else
 [msg,xb,yb,zb,xi,yi] = xyzchk(xb,yb,zb,xi,yi);
end
if length(msg)>0, error(msg); end


 % Calculate quadtree divison into blocks and
 % related objects
[xx,yy,s,nb,lb,ind,bx,by,Na,Lx,Ly,Sp,npb] = ...
           mkblocks(xb,yb,xi,yi,n(1),nmax);

iszero = 0;
if ~isempty(sc)
 iszeros = any(~sc(1:2));
end
if iszero | isempty(sc)
  if lb==1, part_blk = part_blk/sqrt(nb/n(1)); end
  sc(1:2) = max(sc(1:2),[mean(diff(Lx')), mean(diff(Ly'))]*part_blk);
else
  if length(sc) < 3
   sc(3) = 0;         % Phi
  end
end

 % Process basis values vector ZB
szb = size(zb);
is_vec = 0; err_only = 0;
if any(szb==1) | prod(szb)==nb
  zb = zb(:);
elseif min(szb)>1
  if szb(2)==nb & szb(1)~=nb , zb = zb'; end
  is_vec = 1;
elseif zb==[]
  err_only = 1;
end
szb = size(zb);
oz = ones(1,szb(2));

 % Mask for underpopulated blocks .......
isup = npb<n(2);

 % Remove mean and trend ................
if 1
if no_trend
 % 'OBJMAP: no_trend'
 g = mean(zb);
 zb = zb - ones(szb(1),1)*g;
else
 for jj=1:szb(2)
   [zb(:,jj),rc,gc] = detrend2(xb,yb,zb(:,jj));
   r0(jj,:) = rc;
   g(jj,:) = gc;
 end
end
end

 % Initialize output and weights
zi = zeros(prod(size(xi)),szb(2));
wi = zeros(prod(size(xi)),1);
er1 = 1-err;

if is_err, em = zeros(size(xi)); end

 % Tell some parameters and number of blocks ......
if verbose
  fprintf('\n')
  if ~is_call
    fprintf('Scales: %g, %g   ',sc(1),sc(2))
  end
  fprintf('Relative error: %g   ',err) 
  fprintf('Number of blocks: %g\nProcessing: ',lb) 
end

 % Process each block sequentially **********************
for jj=1:lb   % Begin blocks ```````````````````````````0

  [ibp,iip,w,x,y] = ptsinblk(jj,xx,yy,Na,Sp,...
                             Lx,Ly,nb,isup);
  ob = ones(size(ibp));
  oi = ones(size(iip));

  % Green's function (covariance) matrix for basis points
  if is_call
    r = sqrt(x.^2+y.^2);
    eval(call);
  else
    r = exp(-((x*cos(sc(3))+y*sin(sc(3)))/sc(1)).^2 - ...
             ((y*cos(sc(3))-x*sin(sc(3)))/sc(2)).^2);
  end
  r = er1*r+err*eye(size(r));    % Add errors
  A = r;
  if ~err_only, v = r\zb(ibp,:); end

  % Interpolation points ...................
  len = length(iip);
  n_chunk = ceil(len/nmax);
  for j1=1:n_chunk    % Process each chunk

    % Get current chunk
    nn = min(len,j1*nmax);
    i_ch = (j1-1)*nmax+1:nn;
    w_ch = w(i_ch);
    i_ch = iip(i_ch);
    och = ones(size(i_ch));
    x = xx(i_ch,ob)-xx(ibp,och)';
    y = yy(i_ch,ob)-yy(ibp,och)';

    % Calculate covariance matrix
    if is_call
      r = sqrt(x.^2+y.^2);
      eval(call);
    else
      r = exp(-((x*cos(sc(3))+y*sin(sc(3)))/sc(1)).^2 - ...
               ((y*cos(sc(3))-x*sin(sc(3)))/sc(2)).^2);
    end
    r = r*er1;

    % Interpolation itself ...................
    i_ch = i_ch-nb;
    if ~err_only
      zi(i_ch,:) = zi(i_ch,:)+w_ch(:,oz).*(r*v);
    end
    wi(i_ch) = wi(i_ch)+w_ch;

    % Calculate errors
    if is_err & err>0
      r = r';
      em(i_ch) = em(i_ch)+w_ch./(1-sum(r.*(A\r))/er1)';
    end

  end

  % Verbose mode, tell that a block is processed
  if verbose, fprintf('.'); end

end  % End all blocks ''''''''''''''''''''''''''''''''''0

if verbose, fprintf(' All done.\n'); end

zi = zi./wi(:,oz);     % Divide by weights

sz = size(xi);
wi = reshape(wi,sz(1),sz(2));
if is_err & err>0
  em = wi./em;     % Invert error map
end

 % Put mean and slope back ...............
if 1
if no_trend
  zi = zi + ones(size(zi,1),1)*g;
else
 for jj = 1:szb(2)
  zi(:,jj) = zi(:,jj)+g(jj,1)+(xi(:)-r0(jj,1))*g(jj,2);
  zi(:,jj) = zi(:,jj)+(yi(:)-r0(jj,2))*g(jj,3);
 end
end
end

 % Reshape if scalar z value
v=version;
if v(1)=='5'
 if ~err_only
   zi = reshape(zi,sz(1),sz(2),szb(2));
 end
else
 if ~is_vec & ~err_only                                                      
   zi = reshape(zi,sz(1),sz(2)); 
 end
end

 % If XI, YI are not needed
if nargout<=2,
   xi = zi;
  if nargout == 2
    yi = em;
  end
end

