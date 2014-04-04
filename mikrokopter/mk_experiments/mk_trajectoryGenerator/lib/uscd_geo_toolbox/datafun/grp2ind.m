function [ii,nn] = grp2ind(i0,l,s);

% GRP2IND  Built IndexVector from StartIndex and Length
%
% [Index,GroupNumber] = GRP2IND( StartIndex , GroupLength )
%
% where LENGTH(GroupNumber) == MAX(Index)
%
% GRP2IND( ... , LowSampleStep  ) LowSamples in Groups
%
% See also: IND2GRP, GRPMEAN, CUMSUM, CUMPROD
%

ii = [];
nn = [];

if isempty(i0);
   return
end

if nargin < 3
   s = 1;
end

si = size(i0);

if ( sum( si > 1 ) > 1 )
   error('StartIndex must be a Vector.');
end

i0 = i0(:);
l  =  l(:);

if ~isequal(size(i0,1),size(l,1))
   error('Size of StartIndex and GroupLenght must be the same.');
end

jj = ( l == 0 );
if all(jj)
   return
end

if any(jj)
      jj  = find(jj);
   i0(jj) = [];
    l(jj) = [];
end

n = size(l,1);

l = ceil( l / s );

ii = s * ones(sum(l),1);
jj = cumsum( cat(1,1,l) , 1 );

ii(jj(1:n)) = i0;

if n > 1
   ii(jj(2:n)) = ii(jj(2:n))-(i0(1:n-1,1)+s*(l(1:n-1)-1));
end

ii = cumsum(ii,1);

if ( nargout == 2 ) & ~isempty(ii)
   nn = zeros(max(ii),1);
   kk = zeros(size(ii));
   kk(jj(1:n)) = 1;
   kk = cumsum(kk);
   nn(ii) = kk;
end

if n == 1
   return
end

jj = find( si > 1 );

if jj == 1
   return
end
  
perm = cat( 2 , (1:jj-1)+1 , 1 , (jj+1:size(si,2)) );

ii = permute(ii,perm);

if ( nargout == 2 )
   nn = permute(nn,perm);
end