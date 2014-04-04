function c = mhsv(n)

% MHSV  modified Hue-Saturation-Value color map.
%
% ColorMap = MHSV( ColorNumber )
%

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

if n == 0
   c = zeros(0,3);
   return
end


%*********************************************
% Definitions

% Basic Colors

c = [ ...

 0.9  0   1     % m
 0    0   1     % b
 0    1   1     % c
 0    1   0     % g
 1    1   0     % y
 1    0   0     % r

];

% Distances between Colors

d = [ ...

 0.5  % m - m/b
 0.4  %     m/b - b
 0.5  %           b - b/c
 0.4  %               b/c - c
 0.2  %                     c - g/c
 0.2  % g/c - g 
 0.3  %       g - g/y
 0.5  %           g/y - y
 0.5  %                 y - y/r
 0.5  %                     y/r - r

];

%*********************************************
% Scaling of HUE

d = cumsum(cat(1,0,d));
c = rgb2hsv(c);
c = c(:,1);

nc = size(c,1);
nd = size(d,1);

xc = ( 0 : (nc-1) )' / (nc-1);
xd = ( 0 : (nd-1) )' / (nd-1);

c = interp1(xc,c,xd);

d = d / d(nd);     % 0 .. 1

%*********************************************
% Interpolate HUE

x = (0:(n-1))' / (n-1);

c = interp1(d,c,x);

c = hsv2rgb(cat(2,c,ones(n,2))) .^ 0.8;

return

%***************************************************
% RGB-Plot

figure('colormap',c);

hold on, box on, xlim([0.5 n+0.5])

plot(c(:,1),'r');
plot(c(:,2),'g');
plot(c(:,3),'b');

h = rgb2hsv(c);

plot(h(:,1),'k-');

colorbar('horiz');