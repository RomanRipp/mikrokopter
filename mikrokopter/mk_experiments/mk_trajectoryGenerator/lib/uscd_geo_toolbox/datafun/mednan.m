function x = mednan(x,dim)

% MEDNAN   Median value, ignores NaN
%
%   For vectors, MEDNAN(X) is the median value of the elements in X.
%   For matrices, MEDNAN(X) is a row vector containing the median
%   value of each column.  For N-D arrays, MEDNAN(X) is the median
%   value of the elements along the first non-singleton dimension
%   of X.
%
%   MEDNAN(X,DIM) takes the median along the dimension DIM of X.
%
%   Example: If X = [0 1 2
%                    3 4 5]
%
%   then mednan(X,1) is [1.5 2.5 3.5] and median(X,2) is [1
%                                                         4]
%
% MEDNAN expects that SORT( X , 1 , 'ascend' ) returns ending NaN's.
%
%
%   See also MEDIAN, MEAN, STD, MIN, MAX, COV.

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.11 $  $Date: 1997/11/21 23:23:56 $

if nargin==1, 
  dim = min(find(size(x)~=1)); 
  if isempty(dim), dim = 1; end
end

if isempty(x), x = []; return, end

siz = [size(x) ones(1,dim-ndims(x))];
n = size(x,dim);
m = prod(siz)/n;

% Permute and reshape so that DIM becomes the row dimension of a 2-D array
perm = [dim:max(length(size(x)),dim) 1:dim-1];

siz(dim) = 1;

x = reshape(permute(x,perm),n,m);

%-----------------------------------------------------------

% Sort along first dimension

x = sort(x,1,'ascend');

% Check for NaN's

nn = n - sum(isnan(x),1);  % Number of NOT-NaN per Column

ok = ( nn == 0 );  % All NaN's

if all(ok)
   x = NaN * ones(siz);
   return
end

ok = ~ok;  % Not all NaN's

% Check that NaN's  sorted to End of X

if any( ok & isnan(x(1,:)) )
   warning('Leading NaN''s after SORT(ASCEND).');
   x = x(n:-1:1,:);
   if any( ok & isnan(x(1,:)) )
      error('NaN''s didn''t sorted at End or Begin.')
   end
end

%-----------------------------------------------------------
% Median of Columns

ok = find(ok);

nn = nn(ok);

i2 = mod(nn,2);    % [ Even | Odd ]

i1 = ( nn + i2 ) / 2 + n*(ok-1);

i2 = i1 + 1 - i2;  % ( i2 = i1 + 1  if  Even; i2 == i1 if Odd )

x(1,ok) = ( x(i1) + x(i2) ) / 2;  % Mean of X(i1) and X(i2)

%-----------------------------------------------------------

x = ipermute(reshape(x(1,:),siz(perm)),perm);
