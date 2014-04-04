function c = brc(n);


% BRC  Blue->Cyan | Yellow->Red - ColorMap
%
% ColorMap = BRC( ColorNumber )
% 


%**************************************************
% Huge-Base
%
% Magenta  5/6
% Blue     4/6   2/3
% Cyan     3/6  
% Green    2/6   1/3
% Yellow   1/6
% Red      0/6   0/3
%

h = [ 0  1/3  2/3  ];             %    Red -->     Green -->    Blue

h = h + 1/120 * [ -4  -4  12 ];   % MagRed --> YellGreen --> MagBlue

%**************************************************
% Scaling

m = 6/10;
v = 4/5 + 0.3*i;                  % Amplitude  + i*Potenz


%**************************************************
% Get Number of Colors

if nargin == 0

   fig = get( 0 , 'currentfigure' );

   if isempty(fig)
      c = get( 0 , 'defaultfigurecolormap' );
   else
      c = get( fig , 'colormap' );
   end
  
   n = size(c,1);

end

%----------------------------------------------

n = 2 * round(n/2);

n2 = n / 2;

c = ones( n , 3 );


%**************************************************
% Huge

h = h - floor(h);
 
%----------------------------------------------
% Blue --> Yellow :  1 .. n2
 
i1    = h([3 2]);
i1(1) = i1(1) + ( i1(1) < i1(2) );

i1(2) = i1(1) - m * ( i1(1) - i1(2) );


%----------------------------------------------
% Yellow --> Red : n2+1 .. n

i2    = h([1 2]);
i2(1) = i2(1) - ( i2(1) > i2(2) );

i2(2) = i2(1) + m * ( i2(2) - i2(1) );


c(:,1) = cat( 2 , linspace(i1(1),i1(2),n2) , ...
                  linspace(i2(2),i2(1),n2)       )';

%**************************************************
% Value

a = real(v);  % Amplitude
p = imag(v);  % Potenz


if 0

  % Use COSIN, p < 1 

  c(:,3) = cos( 2*pi * linspace(-n2,n2,n)' / n );
  c(:,3) = a * sign(c(:,3)) .* ( abs(c(:,3)) .^ p );  % [ -a .. a ]
 
  % [ -a  a ] --> [ a  1 ]
  c(:,3) = ( ( 1 - a ) / (2*a) ) * ( c(:,3) + a )  + a;    %  [ a .. 1 ]

else

  % Use Potenz, p >= 1


  c3 = linspace(0,n2,n2)';
  c3 = cat( 1 , -c3(n2:-1:1) , c3 );

  c(:,3) = ( 1 - abs( c3 / n2 ) ) .^ p;

  % [ 0  1 ] --> [ a  1 ]

  c(:,3) = ( 1 - a ) *  c(:,3)  + a;    %  [ a .. 1 ]


end

%**************************************************
% HSV --> RGB

c = hsv2rgb( c - (ceil(c)-1) );
