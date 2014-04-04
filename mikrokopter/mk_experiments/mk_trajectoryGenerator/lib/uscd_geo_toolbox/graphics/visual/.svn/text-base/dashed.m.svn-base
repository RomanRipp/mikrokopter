function x = dashed(n,d)

% DASHED  returns a dashed Line between 0 and 1
%
% X = DASHED( N + V*i , D + M*i );
%
%----------------------------------------------------
%
% N:  Number of Dashes between 0 and 1
% V:  Deviation, normalized to Length + Distance
%
% D:  Distance between, normalized to Length of single Dash
%
% M:   -1:  Gab_Number = Dash_Number - 1
%          Dash_Length = 1 / ( N + ( N - 1 ) * D )
%
%       0:  Gab_Number = Dash_Number
%          Dash_Length = 1 / ( N + N * D )
%
%       1:  Gab_Number = Dash_Number + 1
%          Dash_Length = 1 / ( N + ( N + 1 ) * D )
%
% defaults: N = 10; V =  0;
%           D =  1; M = -1;
%
%----------------------------------------------------
%
% X:  [ Start ; End ] = [ 2 by N ] - Matrice
%
%----------------------------------------------------
% Examples:
%----------------------------------------------------
%
% m = [ -1  0  1 ];
%
% figure, axis([ 0  1  -1.5 1.5]), hold on, box on
% set( gca , 'ytick' , m , 'tickdir' , 'out' );
% ylabel('Mode');
%
% for w = m 
%     x = dashed(6,1+w*i);
%     plot(x,w+0*x,'linewidth',3);
% end
%
%----------------------------------------------------
%
% v = ( 0 : 0.1 : 1 );
%
% figure, axis([ 0  1  -0.1  1.1]), hold on, box on
% set( gca , 'ytick' , v , 'tickdir' , 'out' );
% ylabel('Deviation');
%
% for w = v 
%     x = dashed(10+w*i,1-i);
%     plot(x,w+0*x,'linewidth',3);
% end
%
%----------------------------------------------------
%
% figure('units','pixels','position',[100 100 800 200]);
% axis([ 0  1  -1  1 ]),  hold on, box on
% set( gca , 'ytick' , [] , 'tickdir' , 'out' , ...
%            'position' , [ 0.02 0.15 0.96 0.75 ] );
%
% h = line('erasemode','xor','color','k','linewidth',10);
%
% while 1
%    for v = 0.01 : 0.01 : 1 
%        x = dashed(6+v*i,0.5);
%        x = x([1 2 2],:);
%        x(3,:) = NaN;
%        set(h,'xdata',x(:),'ydata',zeros(prod(size(x)),1));
%        pause(0.01)
%    end
% end
%

if nargin < 1
   n = [];
end

if nargin < 2
   d = [];
end

if isempty(n)
   n = 10;
end

if isempty(d)
   d = 1 - i;
end


if ~( isnumeric(n) & ( prod(size(n)) == 1 ) & ...
      isnumeric(d) & ( prod(size(d)) == 1 ) )
    error('Inputs must be single Numerics.');
end

m = imag(n);    % Offset, relative to Length+istance
n = real(n);    % Number of Dash's

if ~( ( mod(n,1) == 0 ) & ( n > 0 ) )
    error('Number must be a positive Integer.');
end

s = sign(imag(d));  % -1 | 0 | 1
d =  abs(real(d));  % Distance, relative to Length

% Length of Dash
l = 1 / ( n + (n+s) * d );

% Distance between StartPoints
d = l * ( 1 + d );

m = m - floor(m);      % Offset: [ 0 .. 1 )

nn = n + 2 * sign(m);  % Add 2 if Offset

x = zeros(2,nn);

x(1,:) = ( 0 : (nn-1) ) * d;
x(1,:) = x(1,:) + ( d - l ) * ( s == 1 );

x(2,:) = x(1,:) + l;

if m == 0 
   return
end

x = x + ( m - 1 ) * d;

x = min(max(x,0),1);

x = x( : , find( x(2,:) > x(1,:)+1e-10 ) );
