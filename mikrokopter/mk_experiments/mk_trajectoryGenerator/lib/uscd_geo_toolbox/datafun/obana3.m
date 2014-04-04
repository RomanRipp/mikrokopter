function zi = obana3(z,x,y,xi,yi,xr,xc,yr,yc)

% OBANA3 fast data gridding using Gaussian weights. 
% 
%   ZI = obana3(Z,X,Y,XI,YI,XR,XC,YR,YC)
%
%   where  ZI          : gridded matrix
%          Z           : data values to be gridded
%          X,Y         : position of data values
%          XI,YI       : new grid points
%          XR        : influence radius in x-direction
%          XC        : cut-off radius in x-direction
%         [YR]       : influence radius in y-direction
%         [YC]       : cut-off radius in y-direction
%         
%   USES: obana2.{m,mexlx,mexaxp,mexsol or other}
%
%	  yr = 'geo' uses geographical distances
%
%   NOTE: The mex-file version is about 4 times faster than
%         the m-file equivalent.
% 
%   See also OBANA, GRIDDATA, OBJMAP, KRIGING.

% adapted from the original Matlab m-file obana.m by Martin Visbeck
% with the additional extensions of obana2.m by Gerd Krahmann
%
% Matlab 4.2c mex-file 
% by Ulf Garternicht in Jun 97 at IfM Kiel
%
% Version 1.1 (Aug 97)

% check number of I/O parameters

Nin = nargin;

if ~any( Nin == [ 7 8 9 ] );
    error('USAGE: zi=obana3(z,x,y,xi,yi,xr,xc,yr,yr)');
end

% make isotropic if no y values are specified

if Nin == 7
   yr=xr;
   yc=xc;
end

if Nin > 7
   if ischar(yr)
     if strcmp(yr,'geo')
       yr=-1;
       yc=0;
     else
       error('Invalid Option.');
     end
   end
end
     
%*******************************************************
% make column vector of input data

x = x(:); 
y = y(:); 
z = z(:);

if ~isequal(size(x),size(y),size(z))
    error('Size of X, Y and Z must be agree.');
end

%*******************************************************
% Check Size of Inputs

sx = size(xi);  px = prod(sx);  vx = ( px == max(sx) );
sy = size(yi);  py = prod(sy);  vy = ( py == max(sy) );

if ~isequal(sx,sy)
    if     vx & vy
       xi = ones(py,1) * xi(:)';
       yi = yi(:) * ones(1,px);
    elseif vx & ( px == sy(2) ) & ( py == sy(1)*sy(2) )
       xi = ones(sy(1),1) * xi(:)';
    elseif vy & ( py == sx(1) ) & ( px == sx(1)*sx(2) )
       yi = yi(:) * ones(1,sx(2));
    else
       error('Size of XI and YI must be agree.');
    end
end

sz = size(xi);

xi = xi(:);
yi = yi(:);

nz = prod(sz);

%*******************************************************
% Check Size of Radien

msg = cell(0,1);

n  = { 'xr' 'xc' 'yr' 'yc' };
ok = ones(size(n));

for ii = 1 : size(n,2)

    zi = eval(n{ii});

    if strcmp(n{ii},'yr') & isequal(zi,-1)
       break 
    end

    s = size(zi);

    ok(ii) = isequal(s,sz)

    if ~ok(ii)
        p = prod(s); v = ( p == max(s) );
        ok(ii) = ( v & any( v == sz ) );
        if ok(ii)
           if v == sz(1)
              zi = zi(:) * ones(1,sz(2));
           else
              zi = ones(sz(1),1) * zi(:)';
           end
        else
           ok(ii) = ( p == 1 );
           if ok(ii)
              zi = zi* ones(sz);
           end
        end
    end

    if ok(ii)
       eval([n{ii} '= zi;']);
    end
end

if any(~ok)
  ii = find(~ok)
  m  = sprintf(1,'%s ',n{ii});
  m = sprintf('Invalid Size of %s.',m);
  error(m)
end

%*******************************************************
% Remove NaN values in Positions

zi = ~( isnan(z) | isnan(x) | isnan(y) );
if any(zi)
   if all(zi)
      zi = NaN * zeros(si);
      return
   else
      zi = find(zi);
      x = x(zi);
      y = y(zi);
      z = z(zi);
   end
end

%*******************************************************

% do the actual function call

if exist('obana2','file') == 3
   zi=obana2(z,x,y,xi,yi,xr,xc,yr,yc);
   return
end


xr=xr(:);
xc=xc(:);
yr=yr(:);
yc=yc(:);


% reset output

zi  = NaN * zeros(1,nz);

% loop over each output value

pc = 0;

fprintf(1,'\nOBANA: ');

for ii = 1 : nz

  % display process
  if ( ii/nz > pc ) 
     fprintf(1,'\rOBANA: %3.0f%%',pc*100);
     pc = pc+0.1; 
  end

  % positions difference
  dx = xi(ii) - x;
  dy = yi(ii) - y;

  % norm with cutoff radius
  jj = sqrt( (dx/xc(ii)).^2 + (dy/yc(ii)).^2 );

  % select only values within cutoff radius
  jj = find( abs(jj) < 1 );

  if ~isempty(jj)

    % norm with inflence radius
    d = sqrt( (dx(jj)/xr(ii)).^2 + (dy(jj)/yr(ii)).^2 );

    % get factors using a gauss distribution
    d = gauss(d,1);

    % sum up values
    s = sumnan(d);

    if s > 0
       zi(ii) = ( d * z(jj) ) / s;
    end

  end

end

fprintf(1,'\n\n');

% reshape to GridSize

zi = reshape(zi,sz);

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [m,varargout] = chkrad(mode,sz,varargin)

Nin  = nargin - 2;
Nout = nargout - 1;

varargout = cell(1,nargout);


name  = { 'xr' 'xc' 'yr' 'yc' };

n = min(Nin,Nout);

ok = ones(1,n);

for ii = 1 : n

    z = varargin{ii};

    if strcmp(n{ii},'yr') & isequal(zi,-1)
       break 
    end

    s = size(z);

    ok(ii) = isequal(s,sz)

    if ~ok(ii)
        p = prod(s); v = ( p == max(s) );
        ok(ii) = ( v & any( v == sz ) );
        if ok(ii)
           if v == sz(1)
              z = z(:) * ones(1,sz(2));
           else
              z = ones(sz(1),1) * z(:)';
           end
        else
           ok(ii) = ( p == 1 );
           if ok(ii)
              z = z* ones(sz);
           end
        end
    end

    if ok(ii)
       varargout{ii} = z;
    end
end

if any(~ok)
  ii = find(~ok)
  m  = sprintf(1,'%s ',n{ii});
  m = sprintf('Invalid Size of %s.',m);
end
