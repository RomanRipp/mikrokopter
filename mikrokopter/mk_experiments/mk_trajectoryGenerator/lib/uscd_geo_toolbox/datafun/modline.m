function y = modline(mode,x,c,A);

% MODLINE   returns modified Line with fixed derivatives at Start and End
%
%  Y = MODLINE( Mode , X , C , Amplitude )
%
%  Mode = 'cos' | 'cubic'
%
%  X = N-Element Vector, should be monotonic
%
%  C = 2-Element Vector with Derivatives at X([1 N])
%
% Boundary Values:
%
%  (dy/dx)(X==X(1)) = c(1);   % Derivative at X(1)
%  (dy/dx)(X==X(N)) = c(2);   % Derivative at X(N)
%
% 'cubic':  y(X==X(1)) = 0;
%           y(X==X(N)) = 0;
%
% 'cos':    y(X==X(1)) = +A;
%           y(X==X(N)) = -A;
%           y(X==X_2 ) = 0;  X_2 = (X(N)-X(1))/2
%     (dy/dx)(X==X_2 ) = -A*pi/(X(N)-X(1));
%           

Nin = nargin;

if ( Nin == 0 )
   mode = '';
end

if Nin < 3
   c = [];
end

if Nin < 4
   A = [];
end

%---------------------------------------

mode = strcmp(lower(mode),'cos');

if isempty(A)
   A = 1;
else
   A = A(1);
end

if isempty(c)
   c = [ 0  0 ];
else
   c = cat( 2  , c(:)' , [ 0  0 ] );
   c = c([1 2]);
end

%---------------------------------------

if Nin < 1
   d = [ 1  pi ];
   d = d(1+mode);
else

  si = size(x);
  ps = prod(si);

  if isempty(x)
     y = zeros(si);
     return
  end

  if ps == 1
     error('X must have Length >= 2.')
  end

  if ~any( ps == si )
     warning('X should be a Vector');
  end

  d = x([1 end]);
  x = ( x - d(1) );
  d = d(2) - d(1);

end

p = { [ 0  c(1)  -(2*c(1)+c(2))/d      (c(1)+c(2))/d^2 ] 
      [ A  c(1)  -(4*c(1)+c(2))/d  (5*c(1)+3*c(2))/d^2  -2*(c(1)+c(2))/d^3 ] };

p = p{ 1+mode };
 
if Nin < 1
   y = p;
   return
end

y = p(1);
for ii = 2 : size(p,2)
    y = y + p(ii) * x.^(ii-1);
end

if mode
   y = y .* cos(pi/d*x); 
end
