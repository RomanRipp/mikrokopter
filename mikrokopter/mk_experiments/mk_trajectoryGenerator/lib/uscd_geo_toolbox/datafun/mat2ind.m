function [c,m] = mat2ind(c)

% MAT2IND  Converts Matrice to Indexed
%
% [ Y , Map ] = MAT2IND( X )
%
%   X   numeric Matrice with max. Dimensions
%
%   Y   Indexed Matrice, refering to Values of Map
%
%  Map  Unique Values of X: [ N  by  size(X,3) ] 
%
% Convert back:
%
%   s1 = size(Y,1); s2 = size(Y,2); s3 = size(Map,2);
%   X  = feval(class(Map),zeros(s1,s2,s3));
%   for ii = 1 : s3
%       X(:,:,ii) = reshape(Map(Y,ii),s1,s2);
%   end
%

if ~( isnumeric(c) & ( ndims(c) <= 3 ) )
    error('Input must be numeric with max. 3 Dimensions.');
elseif isempty(c)
    c = zeros(size(c));
    m = zeros(0,size(c,3));
    return
elseif ~all( isfinite(c(:)) | isnan(c(:)) )
    error('Input must have finite Values or NaN.');
end


cl = class(c);
if ~strcmp(cl,'double')
    c = double(c);
end

s1 = size(c,1);
s2 = size(c,2);
s3 = size(c,3);

n = s1 * s2;

c = reshape(c,n,s3);

si = ( 1 : n )';

for ii = s3 : -1 : 1
    [m,jj] = sort(c(si,ii));
       si  = si(jj);
end

c = c(si,:);

is_nan = find(isnan(c));
if ~isempty(is_nan)
    nn = ceil(max(c(:))) + 1;
    c(is_nan) = nn+1;
    is_nan = 1;
end

c(2:n,:) = c(2:n,:) - c(1:(n-1),:);

if is_nan    
    n = ~( sum( abs(c) < 1e3*eps , 2 ) == s3 );
else
    n = ~( sum(     c  < 1e3*eps , 2 ) == s3 );
end

n(1) = 1;

m = cumsum(c(find(n),:),1);

if is_nan
   m( find( m > nn ) ) = NaN;
end

n = cumsum(n);

c = reshape(n,s1,s2);

c(si) = c;

if ~strcmp(cl,'double')
    m = feval(cl,m);
end
