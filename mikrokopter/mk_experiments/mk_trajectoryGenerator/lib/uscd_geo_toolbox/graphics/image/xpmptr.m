function [Msg,c] = xpmptr(file,fig);

% XPMPTR  Set FigurePointer to XPM-Pointer
%
% [Msg,Data] = XPMPTR( XPM_FileName )
%
%    returns the PointerData of an XPM-File.
%
%  XPMPTR( XPM_FileName , FigureHandle )
%
%    sets the Pointer of the specified Figure(s) customized
%
%  The XPM-File has to define a [ 16 by 16 ]-Pointer-Image:
%   - All Colors different from White will interpred as Black.
%   - The Color None will set to NaN.
%   - The last HexColor gives the PointerHotSpot,
%       the RedValue is the Row, the BlueValue the Column. 
%
%  required M-File:  READXPM
%
%  example for an Pointer-XPM-File:
%
%/* XPM */
%/* Pointer MagnifyGlass */
%static char*glass[]={
%/* width hight num_colors chars_per_pixel */
%"16 16 4 1",
%/* colors */
%". c None",
%"# c #000000",
%"o c #ffffff",
%/* Pointer HotSpot */
%/* Row  00  Column   24-bit-ColorNotation*/
%"p c #060006", 
%/* Row   0  Column */
%/* "p c #606",        8-bit-ColorNotation*/        
%/* pixels */
%"....###.........",
%"..##o.o##.......",
%".#.o...o.#......",
%".#o..#..o#......",
%"#o...#...o#.....",
%"#..#####..#.....",
%"#o...#...o#.....",
%".#o..#..o#......",
%".#.o...o.#o.....",
%"..##o.o####o....",
%"....###.o###o...",
%".........o###o..",
%"..........o###o.",
%"...........o###.",
%"............o#o.",
%"................"};
%


Nin  = nargin;
Nout = nargout;

nl = char(10);

Msg = '';

c = struct( 'cdata'   , { default  } , ...
            'hotspot' , { [ 1  1 ] }     );


if Nin < 1

   Msg = 'Not enough InputArguments.';

   return

end


%-------------------------------------------------------------------------
% Read XPM-File

[Msg,d,m] = readxpm(file);

if ~isempty(Msg)
   return
end


%-------------------------------------------------------------------------
% Check ImageSize

if ~isequal( size(d) , [16 16] )

   Msg = 'Image must have Size [ 16 by 16 ].';

   return

end

%-------------------------------------------------------------------------
% Get CData

nm = size(m,1);

w = find( sum(m,2) == 3 );    % White Color

c.cdata = d./d;

if ~isempty(w)
   for ww = w(:)'
      c.cdata = c.cdata + ( d == ww );   %   White --> 2
   end
end

%-------------------------------------------------------------------------
% Get HotSpot from last Color:  [ Red  Blue ]

p = round( [ 16 ; 16 ; 255 ] * m(nm,[1 3]) );

p(1,:) = 1;

ok = find( sum( ( ( 1 <= p ) & ( p <= 16 ) ) , 2 ) == 2 );

ok = max(ok);

if ok == 1
   warning(['PointerHotSpot in "'  file  '"  out of Range.']);
end

c.hotspot = p(ok,:);


%-------------------------------------------------------------
% Set FigurePointers

if ( Nin == 1 ) 

    if ( Nout == 2 )
       return

    end


   fig = get(0,'currentfigure');

   if isempty(fig)
 
      fig = figure;

   end


end

if isempty(fig)
   return
end


if ~isnumeric(fig)
   Msg = '2. Input must be FigureHandle(s).'
   return
end


ok = find( ishandle(fig) );

if ~isempty(fig)

   set( fig(ok) , 'Pointer'             , 'custom'  , ...
                  'PointerShapeCData'   , c.cdata   , ...
                  'PointerShapeHotSpot' , c.hotspot        );

end




%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function c = default

% DefaultPointerSpapeCData


c = [ ...
'ooo.............'
'o##oo...........'
'o####oo.........'
'.o#####oo.......'
'.o#######oo.....'
'..o########oo...'
'..o##########o..'
'...o#####oooo...'
'...o#####o......'
'....o##oo#o.....'
'....o##o.o#o....'
'.....o#o..o#o...'
'.....o#o...o#o..'
'......o.....o#o.'
'.............o#o'
'..............oo'  ];

c = double(c);

c( find( c == double('.') ) ) = NaN;
c( find( c == double('#') ) ) =  1 ;
c( find( c == double('o') ) ) =  2 ;

