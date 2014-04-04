function x = sawind(n,k)

% SAWIND  Create SAW-Shaped IndexVector
%
% IND = SAWIND( N , K )
%
% Creates an IndexVector [ 1 .. N ] with
%  a branch after K Values.
%
% Example to create a Colormap with high contrasts:
%
% z = peaks(360);
%
% nc = 256; cmap = hsv(nc); si = sawind(nc,6);
%
% figure('colormap',cmap)
% image(z,'cdatamapping','scaled'), colorbar
%
% figure('colormap',cmap(si,:))
% image(z,'cdatamapping','scaled'), colorbar
%

ok = ( isnumeric(n) & ( prod(size(n)) == 1 ) );
if ok
   ok = ( ( n >= 0 ) & ( mod(n,1) == 0 ) );
end
if ~ok
    error('N must be a positive Integer.');
end

ok = ( isnumeric(k) & ( prod(size(k)) == 1 ) );
if ok
   ok = ( ( k >  0 ) & ( mod(k,1) == 0 ) );
end
if ~ok
    error('K must be a positive nonzero Integer.');
end

if n == 0
   x = [];
   return
end

m = ceil(n/k);

x = m * ones(k,m);

x(1,:) = ( 1 : m );

x = cumsum(x,1);

x = x(:);

m = m*k;          % Elements in X

if m == n
   return
end

% Remove the Elements, exeeding N
%
% Last Row, last (m-n) Columns
%

ii = ( 1 : (m-n) ) - 1;  % [ 0 .. m-n-1 ]

x(m-k*ii) = [];
