function x = remfollw(x,z,r)

% REMFOLLW  Remove Following equal Elements from Matrice
%
% X = REMFOLLW( X , Z , R )
%
%  Removes from X the following Values, matching the Elements of Z.
%
%  The single Values will replaced by R (optional)
%

Nin = nargin;

%------------------------------------------------------

if Nin < 1
   x = [];
elseif ~( isnumeric(x) | ischar(x) )
   error('X must be a Numeric or Character-Array.');
end

%------------------------------------------------------

if Nin < 2
   z = [];
elseif ~( isnumeric(z) | ischar(z) )
   error('Z must be a Numeric or Character-Array.');
end

%------------------------------------------------------

if Nin < 3
   r = [];
elseif ~( isnumeric(r) | ischar(r) )
   error('R must be a Numeric or Character-Array.');
elseif ~( isempty(r)  |  ( prod(size(r)) == 1 ) )
   error('R must be a single Element or EMPTY.');
end

%------------------------------------------------------

if isempty(x) | isempty(z)
   return
end

%*******************************************************

m = prod(size(z));

jj = zeros(size(x));

for ii = 1 : m
    
    jj = ( jj | ( x == z(ii) ) );

end

if ~any(jj)
    return
end

%------------------------------------------------------
% IND2GRP

jj = jj(:);

jj = find(jj);

n  = size(jj,1);

ii = cat( 1 , 1 , find( diff(jj,1,1) > 1 )+1 , n+1 );
ll = diff(ii,1,1);
ii = jj(ii(1:end-1));

jj = ( ll > 1 );

if ~any(jj)
    return
end

jj = find(jj);

ii = ii(jj);   % StartIndex
ll = ll(jj);   % Length

%------------------------------------------------------
% GRP2IND

ii = ii + 1;
ll = ll - 1;

m = size(ll,1);

jj = ones(sum(ll),1);
kk = cumsum( cat(1,1,ll) , 1 );

jj(kk(1:m)) = ii;

if m > 1
   jj(kk(2:m)) = jj(kk(2:m))-(ii(1:m-1,1)+(ll(1:m-1)-1));
end

jj = cumsum(jj,1);

%------------------------------------------------------
% Replace and Remove

if ~isempty(r)
    x(ii-1) = r;
end

x(jj) = [];

