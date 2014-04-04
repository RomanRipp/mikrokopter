function [x,si,varargout] = center(p,x0,x,dim,varargin);

% CENTER  Transforms Data into Intervall
%
%  XT = CENTER( Period , Offset , X )
%
% The Values of X will transformed to the Intervall:
%
%   [ Offset   Offset+Period ]
%
%
% [ XT , SortIndex ] = CENTER( Period , Offset , X , DIM )
%
%  Sorts the Transformated Coordinates XT along the Dimension DIM
%   and returns the IndexMatrix. 
%  Either DIM must defined or SortIndex must requested.
%   default: DIM = 2.
%
% [ XT , SortIndex , Y1T , Y2T , ... ] = CENTER( ... , Y1 , Y2 , ... )
%
%  Rearranges the Vectors or Matrices of Y# into the same Order.
%   Y# must correspond with the Size of X, i.e. the Length of
%  Dimension DIM in X and Y# must be the same, the Length of the other
%   Dimensions must be the same like in X, if this Value for X
%   is larger then 1 !
%     
%  example:
%
%   X = [ 1  by Nx ]  ==>  Y# = [ Ny by Nx ]  |  [ 1 by Nx ] 
%   X = [ Ny by Nx ]  ==>  Y# = [ Ny by Nx ] 
%   
%  If they not Correspond, the same Value will returned.
%
%


si = [];

Nin  = nargin;
Nout = nargout;

  is_sort = ( ( Nout >= 2 ) | ( Nin >= 4 ) );

     Nout =   Nout - 2;

varargout = cell(1,Nout);

if Nin < 2
 error('Not enough Input Arguments.');
end
 
if Nin < 3
  x = [];
end

if Nin < 4
  dim = 2;
end

if isempty(x)
  return
end

Nin = min( Nin-4 , Nout );

varargout(1:Nin) = varargin(1:Nin);

if isempty(x)
  return
end

if ~( isnumeric(p)  &  ( prod(size(p)) == 1 ) )
  error('Period must be single Number.');
end

if ~( isnumeric(x0)  &  ( prod(size(x0)) == 1 ) )
  error('Offset must be single Number.');
end

if ~isnumeric(x)
  error('X must be numeric.');
end


if p == 0
  return
end


%----------------------------
% Transform Offset to [ -p/2  p/2 )

%  x0 = x0 - p * floor( (x0+p/2) / p );


%----------------------------
% Remove Offset

  x = x - x0;


%-----------------------------
% Permute DIM --> 1. Dimension

  sx = size(x);

 nsx = size(sx,2);


%-------------------------------
if sx(dim) > 1
 
  %-----------------------------
  % Permute DIM --> 1. Dimension

  perm = cat( 2 , dim , (1:dim-1) , (dim+1:nsx) );

     x = reshape( permute(x,perm) , sx(dim) , ...
                                    prod(sx(perm(2:nsx))) );

  %-----------------------------
  % Check first and last Value
  %  along Dimension for upper PeriodLimit

  i01 = [ 1  sx(dim) ];

  % PeriodLimits
  % p will not matched from Transformation!!!

  il = ( mod(x,p) == 0 );

   y = x;  % Store original Vector
 
end


%---------------------------- 
% Transform to [ 0  p )

  x = x - p * floor( x / p );


%-------------------------------
if sx(dim) > 1

  %-----------------------------
  % Check for LimitValues at p

     ip = ( y(i01+[1 -1],:) < y(i01,:) ) ;

     x(i01,:) = ~il(i01,:) .* x(i01,:) + il(i01,:) .* ( p * ip );

  if  sx(dim) > 2
         
     i01 = ( 2 :  ( sx(dim) - 1 ) );
      
     ip = ( ( y(i01-1,:) < y(i01,:) )  |  ( y(i01+1,:) < y(i01,:) ) ) ;

     x(i01,:) = ~il(i01,:) .* x(i01,:) + il(i01,:) .* ( p * ip );

  end

  %-----------------------------
  % Permute back

     x = reshape( x , sx(perm) );

  perm = cat( 2 , (1:dim-1)+1 , 1 , (dim+1:nsx) );

     x = permute(x,perm);

end

%---------------------------- 
% Add Offset

x = x + x0;


%---------------------------- 

if ~is_sort  |  ( size(x,dim) == 1 )
  return
end


%-----------------------------
% Sort along 2. Dimension

[x,si] = sort(x,dim);

sok = ( sx == 1 );

for ii = 1 : Nin

  v = varargin{ii};

  varargout{ii} = v;

  if ~isempty(v)

    sv  = size(v);
    nsv = size(sv,2);

    ok = ( nsv >= nsx );
    if ok
       ok = all( ( sv(1:nsx) == sx ) | sok );
    end

    if ok
 
      % IndexMatrice for si
      % Valid for 2. Dimension

      ind = cell(nsv,1);

      for jj = 1 : nsv
         if size(x,jj) == 1
           ind{jj} = ones(1,sv(jj));
         else
           ind{jj} = ( 1 : sv(jj) );
         end
      end

      ind = si(ind{:});

      % Build Linear Index

      csv = cat( 2 , 1 , cumprod(sv) );

       iv = 1+0*ind;

      for jj = 1 : nsv

        if jj == dim
          iv = iv + ( ind - 1 ) * csv(jj);
        else
          iv = iv + ( cumsum(1+0*iv,jj) - 1 ) * csv(jj);
        end

      end    

      varargout{ii} = v(iv);

    end      
    % ok

  end
  % ~isempty(v)

end
