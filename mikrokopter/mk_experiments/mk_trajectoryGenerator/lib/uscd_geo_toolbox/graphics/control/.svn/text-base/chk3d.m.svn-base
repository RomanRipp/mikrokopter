function [msg,x,y,varargout] = chk3d(x,y,varargin)

% CHK3D Check arguments to 3-D data routines.
%
%   [MSG,X,Y,Z1,Z2,...] = CHK3D(X,Y,Z1,Z2,...)
%

msg = '';

Nin  = nargin - 2;
Nout = nargout - 3;

if Nin < 0
   error('Inputs X ynd Y are missing.');
end


nn = 0;

if Nout > 0

   varargout = cell(1,Nout);

   if Nin > 0
      nn = min(Nin,Nout);
      varargout(1:nn) = varargin(1:nn);
   end

end


sx = size(x);  px = prod(sx);  vx = ( px == max(sx) );
sy = size(y);  py = prod(sy);  vy = ( py == max(sy) );

%****************************************************************
% X & Y only

if Nin == 0

   if isequal(sx,sy)
      return
   end

   if     vx & vy
      x = ones(py,1) * x(:)';
      y = y(:) * ones(1,px);
   elseif vx & ( px == sy(2) )
      x = ones(sy(1),1) * x(:)';
   elseif vy & ( py == sx(1) )
      y = y(:) * ones(1,sx(2));
   else
      msg = 'Size of X and Y must be agree.';
   end
   
   return

end

%****************************************************************
% Get Size of Matrices

sz = NaN * zeros(Nin,2);

for ii = 1 : Nin
 
    s = size(varargin{ii});

    if size(s,2) <= 3
       sz(ii,[1 2]) = s([1 2]);
    end

end

if any(isnan(sz(:,1)))
   msg = 'Following Inputs must be Matrices with max. 3 Dimensions.';
   return
end

ds = sz - sz(ones(Nin,1),:);

if ~all( ds(:) == 0 )
    pz = prod(sz,2);
    vz = ( pz == max(sz,[],2) );
    if ~all( vz & ( pz == pz(1) ) )
        msg = 'Following Inputs must have the same Size in X and Y.';
        return
    end
    for ii = 1 : nn
          varargout{ii} = varargout{ii}(:);
          if ( sz(1,2) == pz(1) )
             varargout{ii} = varargout{ii}';
          end
    end 
end
 
%****************************************************************
% Check Size of Matrices with X and Y
  
sz = sz(1,[1 2]);  pz = prod(sz);  vz = ( pz == max(sz) );

if isequal(sx,sy,sz)
   return
end

%---------------------------------------------------------------
% Check for Vectors

if vz

   %---------------------------------------
   % X and Y not Vectors
   if ~( vx & vy ) 
      msg = 'When Z is a vector, X and Y must also be vectors.';
      return
   end

   %---------------------------------------
   % Vectors of same Lenght
   if isequal(px,py,pz)
      x = x(:);
      y = y(:);
      for ii = 1 : nn
          varargout{ii} = varargout{ii}(:);
      end
      if ( sz(2) == pz )
         x = x';
         y = y';
         for ii = 1 : nn
             varargout{ii} = varargout{ii}';
         end
      end
      return
   end

   %---------------------------------------
   % X and Y build Z
   if isequal([px py],sz)
      if     ~( px == pz )
          x = x * ones(sz);
      elseif ~( py == pz )
          y = y * ones(sz);
      end
      return
   end

   msg = 'Vectors X and Y must match Size of Vector Z.';
   return

end

%---------------------------------------------------------------
% Check MatriceSize

if     vx & vy
   x = ones(py,1) * x(:)';
   y = y(:) * ones(1,px);
elseif vx & ( px == sy(2) )
   x = ones(sy(1),1) * x(:)';
elseif vy & ( py == sx(1) )
   y = y(:) * ones(1,sx(2));
elseif ~isequal(sx,sy)
   msg = 'Size of X and Y must be agree.';
   return
end

if ~isequal(size(x),size(y),sz)
    msg = 'Size of X and Y must match Size of Z.';
    return
end
