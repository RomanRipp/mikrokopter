function c = mhsv(n)

% MJET  modified NCSA fluid jet image color map.
%
% ColorMap = MJET( ColorNumber )
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

ini = [ 0     0  0  1 
        1/3   0  1  1
        2/3   1  1  0
        1     1  0  0  ];

n1 = n-1;
if n1 == 1
   c  = 0.5;
else
   c = ( 0 : n1 )' / n1;
end

c = interp1(ini(:,1),ini(:,[2 3 4]),c);
 
