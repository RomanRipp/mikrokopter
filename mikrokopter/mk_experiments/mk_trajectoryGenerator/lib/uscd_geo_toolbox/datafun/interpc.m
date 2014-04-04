function x = interpc(x,c);

% INTERPC  Interpolates Values of Matrice to Edges of specified Column
%
% Z = INTERPC( X , ColumnIndex )
%
% see also: INTERPG
%

if nargin < 1
   error('Not enough InputArguments.')
end

if ~isnumeric(x)
    error('X must be numeric.')
end

if nargin < 2
   c = [];
end

if isempty(x)
   return
end

if isempty(c)
   c = 1;
elseif ~( isfinite(c) & ( mod(c,1) == 0 ) & ...
          ( 1 <= c ) & ( c <= size(x,2) )        )
   error('ColumnIndex must be a finite Integer in Range of Columns of X.');
end
 
nx = size(x,1);

if nx <= 1
   return
end

n1 = nx - 1;

i1 = ( 1 : n1 );
i2 = i1 + 1;

dx = x(i2,c) - x(i1,c);

%-----------------------------------------
% Check for BorderCrossings

sg = sign(dx);

x0 =  ceil(sg.*x(i1,c));
x1 = floor(sg.*x(i2,c));

x0 = x0 + ( x0 == sg.*x(i1,c) );
x1 = x1 - ( x1 == sg.*x(i2,c) );

nn = x1 - x0 + 1;
nn = max(nn,0);

if ~any(nn)
    return
end

x0 = sg.*x0;
x1 = sg.*x1;

ns = nn + 1;

i0 = cumsum(cat(1,1,ns));

nz = sum(ns) + 1;

%-----------------------------------------

x           = cat( 1 , x , zeros(nz-nx,size(x,2)) );
x(nz,:)     = x(nx,:);
x(i0(i1),:) = x(i1,:);

for ii = i1(find(ns(i1)>1))

    ix = ( x0(ii) : sg(ii) : x1(ii) )' - x(i0(ii),c);
 
    x(i0(ii)+(1:nn(ii)),:) = ix * ( x(i0(ii+1),:) - x(i0(ii),:) ) / ...
                                  ( x(i0(ii+1),c) - x(i0(ii),c) ) + ...
                                 x( i0(ii)*ones(1,size(ix,1)) , : );

end

