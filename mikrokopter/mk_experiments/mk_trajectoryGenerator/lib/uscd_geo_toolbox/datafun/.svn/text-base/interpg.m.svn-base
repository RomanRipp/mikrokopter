function x = interpg(x)

% INTERPG  Interpolate coordinates of polygon on edges of grid
%
% Z = INTERPG( X )
%
% X = [ N by M ]  N Points with M Coordinates in units [1] of Grid
%
% Z = [ K by M ]  with ( K >= N ) if no duplicate Points of X exist
%
% Z contains the unique Coordinates of X, with added Points
%   where the lines of the polygon X crossing edges of the grid.
%   at the added Points:  any( Z == round(Z) )
%

if nargin < 1
   error('Not enough InputArguments.')
end

if ~isnumeric(x)
    error('X must be numeric.')
end

sx = size(x);
ps = prod(sx);
ns = size(sx,2);

if ps <= 1
   return
end

prm = ( ps == max(sx) );
if prm
   dim = sum(cumprod(double(sx==1))) + 1;
   prm = ( dim > 1 );
   if prm
      perm = cat( 2 , ( dim : ns ) , ( 1 : dim-1) );
         x = permute( x , perm );
        sx = sx(perm);
   end
end

nx = sx(1);
mx = sx(2);

if ~( ps == nx*mx )
    error('X must be a 2-dimensional numeric Matrice');
end

n1 = nx - 1;

i1 = ( 1 : n1 );
i2 = i1 + 1;

dx = x(i2,:) - x(i1,:);

%-----------------------------------------
% Remove duplicate Points

jj = ( sum( abs(dx) < 1e2*eps , 2 ) == mx );

if any(jj)

   if nx-sum(jj) < 2
      if prm
         x = permute(x,perm);
      end
      return
   end

   jj = find(jj);

   x(jj+1,:) = [];

   nx = size(x,1);

   n1 = nx - 1;

   i1 = ( 1 : n1 );
   i2 = i1 + 1;

   dx = x(i2,:) - x(i1,:);

end
   
%-----------------------------------------
% Check for BorderCrossings

sg = sign(dx);

x0 =  ceil(sg.*x(i1,:));
x1 = floor(sg.*x(i2,:));

x0 = x0 + ( x0 == sg.*x(i1,:) );
x1 = x1 - ( x1 == sg.*x(i2,:) );

nn = x1 - x0 + 1;
nn = max(nn,0);

if ~any(nn(:))
    if prm
       x = ipermute(x,perm);
    end
    return
end

x0 = sg.*x0;
x1 = sg.*x1;

ns = sum(nn,2) + 1;

i0 = cumsum(cat(1,1,ns));

nz = sum(ns) + 1;

%-----------------------------------------
% Add Distance

x  = cat( 2 , x , cumsum(cat(1,0,sqrt(sum(dx.^2,2)))) );

%-----------------------------------------

x           = cat( 1 , x , zeros(nz-nx,mx+1) );
x(nz,:)     = x(nx,:);
x(i0(i1),:) = x(i1,:);

for ii = i1(find(ns(i1)>1))

    i1 = cumsum(cat(2,1,nn(ii,:))) + i0(ii) - 1;

    kk = ~( nn(ii,:) == 0 );

    for jj = find(kk)

        ix = ( x0(ii,jj) : sg(ii,jj) : x1(ii,jj) )' - x(i0(ii),jj);
 
        x(i1(jj)+(1:nn(ii,jj)),:) = ix * ( x(i0(ii+1),:)  - x(i0(ii),:)  ) / ...
                                         ( x(i0(ii+1),jj) - x(i0(ii),jj) ) + ...
                                     x( i0(ii)*ones(1,size(ix,1)) , : );

    end

    if sum(kk) > 1

       ind = i0(ii) + ( 1 : ns(ii)-1 );

       [h,si] = sort(x(ind,mx+1));

       x(ind,:) = x(ind(si),:);

    end

end


%-----------------------------------------
% Remove duplicate Points

x = x(:,1:mx);

dx = diff(x,1,1);

jj = ( sum( abs(dx) < 1e2*eps , 2 ) == mx );

if any(jj)

   jj = find(jj);

   x(jj+1,:) = [];

end

%-----------------------------------------
   
if prm
   x = ipermute(x,perm);
end
