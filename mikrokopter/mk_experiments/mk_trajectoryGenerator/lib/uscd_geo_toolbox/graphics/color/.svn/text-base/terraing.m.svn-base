function cmap = terraing(N,p);

% TERRAING   TerrainColormap  [ R  G  B ], GrayScaled LandColors
%            
%           DarkCyan --> LightCyan / DarkGray --> LightGray
%
% Palette = TERRAING(N)  returns N Colors (default: 128) 
%          
%  Note:  N is automaticly set to an EVEN Number
%
% See also: TERRA, TERRAIN
%

if nargin < 1
   f = get(0,'currentfigure');
   if isempty(f)
      N = size(get(0,'defaultfigurecolormap'),1);
   else
      N = size(get(f,'colormap'),1);
   end
end

if nargin < 2
   p = 2/3 + 0.2*i;
end
     
  N = 2*ceil(N/2);

cmap = zeros(N,3);

c0 = [ 0.40  0.60  0.75
       0.85  0.90  0.90 ];

c1 = [ 0.60  0.60  0.60
       0.85  0.85  0.85  ];

br = imag(p);

c0 = cat( 3 ,  ( ( 1 - br ) * c0 + br ) , c0*br );
c1 = cat( 3 ,  ( ( 1 - br ) * c1 + br ) , c1*br );
 
ii = 1 + ( br < 0 );

c0 = c0(:,:,ii);
c1 = c1(:,:,ii);
  
ind = ( 1 : N/2 );

x0  = [ 0   1 ]';
x1  = linspace(x0(1),x0(2),N/2)';

x1 = x1 .^ real(p);

cmap(ind,:)     = interp1(x0,c0,x1);

cmap(ind+N/2,:) = interp1(x0,c1,x1);



