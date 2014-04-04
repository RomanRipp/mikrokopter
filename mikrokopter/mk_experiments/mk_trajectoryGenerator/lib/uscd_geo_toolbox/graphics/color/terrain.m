function cmap = terrain(n);

% TERRAIN   TerrainColormap  [ R  G  B ]
%            
%           DarkCyan --> LightCyan / Green --> Brown
%
% Palette = TERRAIN(N)  returns N Colors
%          
%  Note:  N is automaticly set to an EVEN Number
%
% see also: TERRA, TERRAING
%

cmap = [ ...
   0.     10.  30.  90.
  10.     20.  40.  95.
  49.999  80. 100. 100.
  50.001  20.  90.   0.
  55.     90.  80.   0.
  90.     70.  15.   0.
 100.     65.  30.  20. ];

if nargin == 0

   fig = get( 0 , 'currentfigure' );

   if isempty(fig)
      c = get( 0 , 'defaultfigurecolormap' );
   else
      c = get( fig , 'colormap' );
   end
  
   n = size(c,1);

else

  ok = ( isnumeric(n)  &  ( prod(size(n)) == 1 ) );
  if ok
     ok = ( ( n >= 0 ) & ( mod(n,1) == 0 ) );
  end

  if ~ok
      error('ColorNumber must be a positive Integer.');
  elseif ( n == 0 )
      cmap = zeros(0,3);
      return
  end

end

n = 2*ceil(n/2);

c = linspace(cmap(1,1),cmap(end,1),n)';

cmap = interp1( cmap(:,1) , cmap(:,[2 3 4]) , c ) / 100;


