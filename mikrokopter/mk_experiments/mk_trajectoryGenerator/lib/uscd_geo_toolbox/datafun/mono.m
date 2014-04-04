function ii = mono(x);

% MONO monotonize vector
%
% Index = MONO( V )
%
%  returns Index of the ascending elements of V
%   i.e.  all( diff(V(Index)) > 0 ) == 1
%  

x = x(:);

ok = ones( size(x) );

jj = 0;

while ~isempty(jj) & ( sum(ok) > 1 )  

   ii = find(ok);

   jj = find( diff(x(ii)) <= 0 ) + 1;

   ok(ii(jj)) = 0;

end

ii = find(ok);