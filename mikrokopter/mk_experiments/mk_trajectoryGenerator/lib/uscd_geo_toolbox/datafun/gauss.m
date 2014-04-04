function y=gauss(x,s,m)

% GAUSS  Gaussian function
%
%  Y = GAUSS( X , S , M )
%
%  Y = EXP(-(X-M).^2./S.^2)./(sqrt(2*pi).*S);
%
%  sum( Y(X=(-inf..inf)) * dX ) = 1/sqrt(2)
%

Nin = nargin;

if Nin < 2, s = 1; end
if Nin < 3, m = 0; end

x = ((x-m).^2) ./ (s.^2);

s = sqrt(2*pi) * s;

y = exp(-x) ./ s;
