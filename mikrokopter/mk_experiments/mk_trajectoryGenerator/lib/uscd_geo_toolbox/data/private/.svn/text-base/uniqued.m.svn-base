function [x,ix,nn,iy] = uniqued(x,d);

% UNIQUED Returns a Matrice with no repetitions along a specified Dimension
%
%  [ Y , IX ] = UNIQUED( X , [DIM] )
%
%  If DIM is zero or empty the first nonsingleton Dimension of X is used.
%
%  IX is the corresponding IndexVector for Y in X at the Dimension DIM,
%    that Y == X( : , ... , IX , ... , : )
%
%  Use an NonZero imaginary part for DIM to sort the Matrix along DIM.
%
%  A positive imaginary part of DIM returns the sorted result,
%  a negative imaginary part returns the result in original order, 
%   where following equal Elements are removed.
%
%  [ Y , IX , NY , IY  ] = UNIQUED( X , [DIM] )
%
%  NY is the Number of repetitions of the Y-Elements in X,
%
%  IY the corresponding IndexVector for X in Y:
%    that X == Y( : , ... IY , ... , : );
%
% see also: UNIQUE, REMFOLLW
%

Nin = nargin;

isy = ( nargout >= 4 );

%--------------------------------------------------------
% Check Inputs

Nin = nargin;

if Nin < 1
   x = [];
elseif ~( isnumeric(x) | ischar(x) )
   error('X must be a Numeric or Character-Array.');
end

if Nin < 2
   d = [];
else
   ok = ( isnumeric(d) & ( prod(size(d)) == 1 ) );
   if ok
      rd = real(d);
      ok = ( ( mod(rd,1) == 0 ) & ( rd >= 0 ) );
   end
   if ~ok
       error('DIM must be a positive Integer.');
   end
end

%--------------------------------------------------------
% Check with DIM

ix = [];
iy = [];

nn = [];

if isempty(x)
   return
end

if isempty(d)
   d = 0;
end

id = imag(d);
d  = real(d);

srt = ~( id == 0 );  % True for Sort
str = ischar(x);     % True for Character

d = real(d);

s = size(x);

if d == 0
   d = sum(cumprod(double(s==1))) + 1;  % First with not ONE
elseif size(x,d) == 1
   ix = 1;
   iy = 1;
   nn = 1;
   return
end

%--------------------------------------------------------
% Permute that DIM --> 1. Dimension

m = size(s,2);

p = cat( 2 , d , ( 1 : d-1 ) , ( d+1 : m ) );

s1 = s(p(1));
s2 = prod(s(p(2:m)));

x = reshape(permute(x,p),s1,s2);

%--------------------------------------------------------
% Sort

ix = ( 1 : s1 )';

if srt

   for ii = s2 : -1 : 1
       [y,jj] = sort(x(ix,ii));
          ix  = ix(jj);
   end

   x = x(ix,:);

end

%--------------------------------------------------------
% Check X with Differences 

ii = ( 1 : (s1-1) );

jj = isnan(x); 

jj = ( ( diff(x,1,1) == 0 ) | ( jj(ii,:) & jj(ii+1,:) ) );
jj = ( sum( jj , 2 ) < s2 );


if ~all(jj)

    jj = find(cat(1,1,jj));  % True if Difference

    if isy
       iy     = zeros(s1,1);
       iy(jj) = 1;
       iy     = cumsum(iy,1);
       if srt
          if id < 0
             [h,ii] = sort(ix(jj));
                iy  = ii(iy);
          end
          iy(ix) = iy;
       end
    end

    x  = x(jj,:);

    ix = ix(jj);

    nn = diff(cat(1,jj,s1+1),1,1);

else

    nn = ones(s1,1);

end

%--------------------------------------------------------
% Sort back if negative Imaginary part of DIM

if id < 0
   [ix,jj] = sort(ix);
    x   =   x(jj,:);
    nn  =  nn(jj);
    if isy
       iy = jj(iy);
    end
end

%--------------------------------------------------------
% Reshape and Permute back

s(p(1)) = size(x,1);

x = reshape(x,s(p));

p = cat( 2 , ( 1 : d-1 ) + 1 , 1 , ( d+1 : m ) );

x = permute(x,p);

if d == 2
   ix = permute(ix,[2 1]);
   nn = permute(nn,[2 1]);
   if isy
      iy = permute(iy,[2 1]);
   end
end
