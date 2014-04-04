function c = rgb(n)

% RGB  Red-Yellow-Green-Cyan-Blue Colormap
%
% ColorMap = RGB( N )
%

if nargin < 1
  fig = get(0,'currentfigure');
  if isempty(fig)
    n = size( get(0,'defaultfigurecolormap') , 1 );
  else
    n = size( get(fig,'colormap') , 1 );
  end
end


n0 = 64;

x = linspace( 0 , pi/2 , n0+1 );
y = sin(x) + 1;
y = cat( 2 , ( 2 - y(end:-1:2) ) , y );
y = cat( 2 , y , y(end-1:-1:1) );
y = y/2;

y = y(:);

ny = size(y,1);  %  n0 == 4*n0+1  !!!

n2 = ( ny - 1 ) / 2;

cc = [ 0  0  1   -1
       0  1  1    0
       0  1  0    1
       1  1  0    2
       1  0  0    3   ];

no = n2/2;  % Offset for Red and Blue
% Overlapp: n2
% Red and Blue n2+no !!!
n0 = 4*n2 + 1 + 2*no;
c0 = zeros( n0 , 3 );

ic = ( 1 : ny );

for ii = 1 : size(cc,1);

  jj = cc(ii,4) * n2 + no + ic;

  ok = find( ( 1 <= jj )  &  ( jj <= n0 ) );
 
  for kk = 1 : 3
    c0(jj(ok),kk) = c0(jj(ok),kk) + y(ok) * cc(ii,kk);
  end

end

c0 = c0 + ( 1 - c0 ) .* ( c0 > 1 );

c0 = c0 .^ 0.5;

c = interp1( ( 1 : n0 )' , c0 , linspace( 1 , n0 , n )' );