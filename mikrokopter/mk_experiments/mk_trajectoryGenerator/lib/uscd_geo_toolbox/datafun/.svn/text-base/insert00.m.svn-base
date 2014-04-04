function [x,ind,ix] = insert(x,ind,y,dim)

% INSERT Inserts Values into Array
%
% [ Z , IY , IX ] = INSERT( X , Index , Y , DIM )
%
% inserts Y into X at the Index in the Dimension DIM.
%
% defaults: DIM = 1, for Vectors X the first non-singleton Dimension
%
% Index >  0  ==>  Y will inserted behind X(Index)
% Index <  0  ==>  Y will inserted before X(Index)
% Index == 0  ==>  Y will inserted before X(1)
%
% Make sure, that Index match the Size of X in DIM.
% 
% Y can be Empty, a Scalar, Vector or Matrice.
%
%  Make sure, that the Size of Y match the Size of X and 
%   the Length of Index in DIM.
%
% If Y is empty or not given, X(Index) will inserted (duplicated).
%
%--------------------------------------------------------------
% The Output-IndexVectors IY and IX refers to the inserted Y 
%   and the original X.
%
% IY = [ N by 2 ];  with ( N <= length(Index) )  
%
%   The 1. Column of IY is the Index in Z of the Inserted  Y
%   The 2. Column of IY is the Index in Y of the Inserted  Y
% 
% IX = [ size(X,DIM) by 1 ] is the Index in Z of the Original X
%
%--------------------------------------------------------------
% for Vectors X and Vectors or Scalars Y
%
%  Z( IY(:,1) ) == Y( IY(:,2) ) 
%
%  Z(IX)        == X
% 
%--------------------------------------------------------------
% for Matrice X and Y
%
%  Z( : , ... , IY(:,1) , ... , : ) == Y( : , ... , IY(:,2) , ... , : ) 
%
%  Z( : , ... , IX , ... , : )      == X
%


Nin = nargin;

Nout = nargout;

ix = [];

if Nin < 2
   error('Not enough Input Arguments.');
end

if Nin < 3
   y = [];
end

if Nin < 4
   dim = [];
end

%--------------------------------------------------------
% Check Size of X and Index

sx = size(x);

if isempty(dim)
   [ms,dim] = max(sx); 
   if ~( ms == prod( sx + ( sx == 0 ) ) )  % Not a Vector
      dim = 1;
   end
else
   sx = cat( 2 , sx , ones(1,dim-size(sx,2)) );
end

if ( all(sx) == 0 )
   sx      = ones( 1 , dim+(dim==1) );
   sx(dim) = 0;
   x       = zeros(sx);
end

ind = floor(ind(:));
 ni =  size(ind,1);

ind = abs(ind) - 1 * ( ind < 0 );  % Negative Inserts before !!!

if ~all( ( 0 <= ind ) & ( ind <= sx(dim) ) );
   error(sprintf('Index exceeds Matrix Dimension %.0f.',dim));
end

%--------------------------------------------------------
% Check Size of Y

sy = size(y);

not_y =  all( sy == 0 );   % ZERO-Size
one_y = prod( sy == 1 );   % Scalar

vec_y = ( max(sy) == prod( sy + ( sy == 0 ) ) ); % Scalar | Vector
vec_y = ( vec_y & ~one_y );                      % Vector

mat_y = ~( not_y | one_y );                      % Vector | Matrice

%--------------------------------------------------------
% Vector ==> reshape into Dimension DIM

if vec_y

   y = y(:);

   if dim > 1

      y = permute( y , cat( 2 , ( 2 : dim ) , 1 ) );

   end

   sy = size(y);

end

%--------------------------------------------------------
% Check Size of Y with Size of X and Length of Index at DIM

if mat_y

    sy = cat( 2 , sy , ones(1,dim-size(sy,2)) );

    s      = sx;
    s(dim) = ni;

    if ~isequal( s , sy )
       error(sprintf('Size of Y must match Size of X and Length of Index in Dimension %.0f.',dim));
    end

end


%--------------------------------------------------------
% Permute and Reshape X and Y to 1. Dimension at DIM

ns = size(sx,2);

perm = cat( 2 , dim , ( 1 : dim-1 ) , ( dim+1 : ns ) );

x = reshape( permute(x,perm) , sx(perm(1)) , prod(sx(perm(2:ns))) );

if mat_y
   y = reshape( permute(y,perm) , sy(perm(1)) , prod(sy(perm(2:ns))) );
end

%--------------------------------------------------------
% Remove Multiple Indize

[ind,si] = sort(ind,1);

    bad  = find( diff(ind,[],1) == 0 );

ind(bad) = [];

ni = size(ind,1);

if mat_y
     y        = y(si,:);
     y(bad,:) = [];
    si(bad)   = [];
else
    si        = ones(ni,1);
    if not_y
       si = NaN * si;
    end
end

%--------------------------------------------------------
% Build IndexVector to Insert

ii = ones( sx(dim)+ni , 1 );

ind = ind + ( 1 : ni )';

ii(ind) = 0;

if Nout > 2
   ix = find(ii);
end

if ~isempty(ii)

  ii    = cumsum(ii,1);

  ii(1) = ii(1) + ( ii(1) == 0 );   % Check for ZERO-Index

end

% Insert X (duplicate)
if ~isempty(x)
   x = x(ii,:);
end

% Insert Y
if ~not_y
   x(ind,:) = y;
end

if Nout > 1
   ind = cat( 2 , ind , si );
end

%--------------------------------------------------------
% Reshape and Permute back

sx(dim) = sx(dim) + ni;

x = reshape( x , sx(perm) );

perm = cat( 2 , ( 1 : dim-1 ) + 1 , 1 , ( dim+1 : ns ) );

x = permute( x , perm );
