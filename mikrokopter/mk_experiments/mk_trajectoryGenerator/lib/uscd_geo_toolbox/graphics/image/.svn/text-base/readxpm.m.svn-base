function [Msg,c,cmap,ccm] = readxpm(file);

% READXPM  Reads ImageData from XPM-Format
%
% [ Msg , C , Colormap , CharacterColorMap ] = READXPM(FileName)
%
%  returns the Indexed Colormatrice C and the corresponding ColorMap.
%
%
% [ Msg , C ] = READXPM(FileName)
%
%   returns the 3-dimensional True-ColorMatrice C, contains
%    the ColorValues  ( 0 <= C <= 1 );
%
% READXPM supports multiple Frames, concatenated in the XPM-File.
%
%  The Frames will added to the 3. Dimension of C, in case of indexed Colors,
%                        or the 4. Dimension of C, in case of RGB-True-Colors.
%
%
% Note:  READXPM  requires Hexadecimal ColorDefinitions,
%                  like:    '#ffffff'
%
%                 or X11-ColorNames,
%                  like:    'white'
%
%  The X11-ColorNames are defined by the function RGB at the end of READXPM
%
%  see also: WRTXPM, READPPM, WRTPPM, IMREAD, IMWRITE
%
%----------------------------------------------------------------------------
%  example for an XPM-File:
%
%/* XPM */
%/* MagnifyGlass */
%static char *glass[]={
%/* width hight Ncolor chars_per_pixel Nframe */
%"16 16 4 1 1",
%/* colors */
%". c None",
%"# c #000000",
%"o c #ffffff",
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
%"................"
%};
%


Msg    = '';
c      = [];
cmap   = zeros(0,3);
ccm    = [];

nl  = char(10);

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
Msg0 = sprintf('%s: ',upper(fcn));

Nout = nargout;


if nargin < 1
  Msg = [ Msg0  'Input FileName is undefined.' ];
  return
end

if ~ischar(file) | isempty(file) | ( prod(size(file)) ~= size(file,2) )
  Msg = [ Msg0  'FileName must be a string.' ];
  return
end


% File starts with:

% /* XPM */
% static char *name[] = {
% /* width hight Ncolor chars_per_pixel Nframe */
% "    70    70   7      1    1",
% /* colors */
% ". c #c0c0c0",
% 

% ColorIdentifer
cid = 'c ';

%********************************************************************
% Open & Read File

fid = fopen(file,'r');

if fid == -1
 file1 = cat( 2 , file(1:(end-strcmp(file(end),'.'))) , '.xpm' );
 fid = fopen(file1,'r');
 if fid == -1
   Msg = [ Msg0  'Can''t open file: '  file ];
   return
 end
 file = file1;
end

%-------------------------------------------------
% Get First Line and Check Header

bb = fgetl(fid);

if ~ischar(bb)
 if  bb == -1
   Msg = [ Msg0  'Empty file: '  file ];
   return
 end
end
 
header = '/* XPM */';
nh     = size(header,2);
nb     = size(bb,2);
n      = min(nh,nb);

if ~strcmp(bb(1:n),header(1:n))
 Msg = [ Msg0  'Invalid Header, '''  header  ''' expected.' ];
 return
end
  
%-------------------------------------------------
% Read Following Data

bb = fread(fid,'char');

fclose(fid);

%********************************************************************

bb = bb(:)';

%-----------------------------------------------------------------
% Remove CR

bb( find( bb == 13 ) ) = [];

nb = size(bb,2);

%-----------------------------------------------------------------
% Remove Blanks, leading NewLine

jj = 1;

while ~isempty(jj)

   jj = find( ( bb(1:end-1) == 32 ) & ( bb(2:end) == 10 ) );

   bb(jj) = [];

end

nb = size(bb,2);

%-----------------------------------------------------------------
% Get Start-EndIndex of DataLines: " .. "

i01 = find( bb == double('"') );

ni = prod(size(i01));

ok = zeros(1,ni);

%-------------------------------------------------------------
% Check Previous Character of Start-EndIndex: NewLine

is1 = ( i01(1) == 1 );

if is1
   ok(1) = ok(1) + 1;
end

ind = ( (1+is1) : ni );

ok(ind) = ( bb(i01(ind)-1) == 10  );  % Previous is NewLine

%-------------------------------------------------------------
% Check Following Characters of Start-EndIndex: NewLine or ',\n'

% Following is NewLine

ind = ( 1 : (ni-1) );

ok(ind) = ok(ind) + ( bb(i01(ind)+1) == 10 ); 


% Following is ',\n'

is2 = ( i01(ni-1) == nb-1 );

ind = ( 1 : (ni-1-is2) );

ok(ind) = ok(ind) + ( ( bb(i01(ind)+1) == double(',') ) & ( bb(i01(ind)+2) == 10 ) );  

% The last should be allways ok

ok(ni) = 1;  

%-------------------------------------------------------------
% Only ONE ok is allowed

ok = ( ok == 1 );

i01 = i01( find( ok == 1 ) );

%-------------------------------------------------------------

if ( mod(i01,2) ~= 0 )  |  isempty(i01)
 Msg = [ Msg0  'Syntax Error.' ];
 return
end


i01 = reshape(i01,2,length(i01)/2)';

i01(:,1) = i01(:,1)+1;  % Start  of '"'
i01(:,2) = i01(:,2)-1;  % End    of '"'

%********************************************************************

%-----------------------------------------------------------------
% Read Initialisation from 1. Line
% /* width height num_colors chars_per_pixel */

msg = '';
try
  eval(['ini = ['  char(bb(i01(1,1):i01(1,2)))  '];'] );
catch
  msg = lasterr;
end

ok = isempty(msg);
if ok 
   ini = ini(:)';
   ok = ( isnumeric(ini) & ( size(ini,2) >= 4 ) & ~isempty(ini) );
   if ok
      ok = all( isfinite(ini) & ( ini > 0 ) & ( mod(ini,1) == 0 ) );
   end
end


if ~ok
  Msg = [ Msg0 'Invalid Initalisation for [ Width Hight NumColors CharPerPixel ]. ' ];
  return
end 

si  = ini([2 1]);  % [ Height  Width ]
nc  = ini(3);      % ColorNumber
cpp = ini(4);      % CharPerPixel

%********************************************************************
% Check Color- and PixelNumber

npr = size(i01,1) - 1 - nc;   % Number pf Pixel-Rows 

if npr < si(1)
  Msg = [ Msg0 'Invalid Number of Colors and Pixel. ' ];
  return
end 

nf = floor( npr / si(1) );  % Number of Frames

if nf*si(1) < npr

   ww = warnstat;

   warning('on');

   fprintf(nl)
   warning(['Incomplete Frames in XPM-File: '  file  ]);
   fprintf(nl)

   warning(ww);

   i01 = i01( 1 : 1+nc+nf*si(1) , : );
end

npr = nf * si(1);

ind = 1 + nc + ( 1 : npr );

if ~all( i01(ind,2)-i01(ind,1)+1 == cpp*si(2) )
  Msg = [ Msg0 'Invalid Number of Pixel per Row. ' ];
  return
end 


p = 10 .^ ( 3 * (0:cpp-1) );  % Potenz for add multiple PixelCharacter


%-----------------------------------------------------------------
% Get Colors

if ~all( i01(2:nc+1,2)-i01(2:nc+1,1)+1 >= cpp+1 )
  Msg = [ Msg0 'Error reading Colormap.' ];
  return
end

cmap = NaN*zeros(nc,3);
cini = NaN*zeros(nc,2);   % [ IniValue  CmapIndex ]

crgb = rgb;  % { ColorName    ColorHex     ColorDec        }
             %   'white'      '#FFFFFF'   [ 255 255 255 ] 

crgb(:,1) = lower(crgb(:,1));

% CharacterMap

ccm = bb( ones(nc,1)*(1:cpp) + i01((1:nc)+1,ones(1,cpp)) - 1 );

if isequal(size(ccm),[cpp nc])
   ccm = permute(ccm,[2 1]);
end

cini(:,1) = sum( ccm.*p(ones(nc,1),:) , 2);

i01( 1 + (1:nc) , 1 ) = i01( 1 + (1:nc) , 1 ) + cpp;

for ii = 1 : nc

  ind = ( i01(ii+1,1) : i01(ii+1,2) );
 
  % Search for HexColor
  jj = find( bb(ind) == double('#') );

  if isempty(jj)
 
     cc = cat(2,32,bb(ind),32);    % ColorString, Blanks added
    
     % Remove Single Characters

      ch = ~( cc == 32 );  % NonBlanks
     dch = diff(ch);
     is1 = ( ( dch(1:end-1) ==  1 )  &  ... 
             ( dch(2:end-0) == -1 )         );
     cc(find(is1)+1) = [];
        
     % Remove Blanks from left and right
     jj = find( ( cc == 32 )  |   ( cc == 9 ) );
     if ~isempty(jj) 
        jj1 = find( jj == ( 1 : length(jj) ) );
        jj2 = find( jj == ( size(cc,2)-length(jj)+1 : size(cc,2) ) );
        cc(jj([jj1 jj2])) = []; 
     end

     % Check ColorName  or 'none'
     kk  = find( strcmp( lower(char(cc)) , crgb(:,1) ) );

     if isempty(kk)  &  isempty( findstr( lower(char(cc)) , 'none' ) )

       Msg = [ Msg nl(1:(end*(~isempty(Msg))))  ...
               'Color "'  char(bb(ind)) '" not supported.' ];

     elseif ~isempty(kk)

        kk = kk(1);

        cc = crgb{kk,2}(2:end);
        s  = size(cc,2)/3;

        cmax = hex2dec(char(double('f')*ones(1,s)));
        for jj = 1 : 3
          cmap(ii,jj) = hex2dec( char( cc( (jj-1)*s+1 : jj*s ) )  ) / cmax;       
        end

        cini(ii,2) = ii;

     end


  else

     cc = bb(ind(jj)+1:ind(end));
     
     % Remove Blanks from left and right
     jj = find( ( cc == 32 )  |   ( cc == 9 ) );
     if ~isempty(jj) 
        jj1 = find( jj == ( 1 : length(jj) ) );
        jj2 = find( jj == ( size(cc,2)-length(jj)+1 : size(cc,2) ) );
        cc(jj([jj1 jj2])) = []; 
     end

     % Step for single Color
     s = size(cc,2) / 3;

     if ~any( s == [ 1  2  4 ] )
       Msg = [ Msg0 'Invalid Color. ' ];
       return
     end

     cmax = hex2dec(char(double('f')*ones(1,s)));
     for jj = 1 : 3
       cmap(ii,jj) = hex2dec( char( cc( (jj-1)*s+1 : jj*s ) )  ) / cmax;       
     end

     cini(ii,2) = ii;

  end

end

if ~isempty(Msg)
  Msg = [ Msg0 'Error reading Colormap.'  nl  Msg  ];
  return
end

jj = find(~isnan(cini(:,2)));
cini = cini(jj,:);
cmap = cmap(jj,:);
ccm  =  ccm(jj,:);

cini(:,2) = cumsum(ones(size(cini,1),1),1);


%-----------------------------------------------------------------
% Get Pixel

ind = 1 + nc + ( 1 : npr );

ii      = ones(npr,cpp*si(2));
ii(:,1) = i01(ind,1);

ii = cumsum(ii,2);

c = bb(ii);

c = permute(c,[2 1]);            % [ cpp*si(2)    by npr ]
c = reshape(c,cpp,si(2),npr);    % [ cpp by si(2) by npr ]
c = permute(c,[3 2 1]);          % [ npr by si(2) by cpp ]

p = permute(p(:),[2 3 1]);          % [ 1 by 1 by cpp ]

c = sum( c .* p(ones(1,npr),ones(1,si(2)),:) , 3 );

c0 = c;
c  = NaN*c;

for ii = 1 : size(cini,1)

  c(find(c0==cini(ii,1))) = cini(ii,2);

end

%-----------------------------------------------
% RGB-Image

if Nout < 3

 cmap = cat( 1, cmap , NaN*ones(1,3) );

 n = size(cmap,1);

 c(find(isnan(c))) = n;

 c = cmap(cat(3,c,c+n,c+2*n));

end

%-----------------------------------------------
% Multiple Frames

if nf > 1

   sc = size(c);

   nd = size(sc,2);

   perm = cat( 2 , ( 2 : nd ) , 1 );
   
   sc   = cat( 2 , sc(2:nd) , si(1) , nf );

   c = reshape( permute( c , perm ) , sc );

   perm = cat( 2 , nd , ( 1 : nd-1 ) , nd+1 );

   c = permute( c , perm );

end

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ww = warnstat

% WARNSTAT  Returns global WarningStatus
%
%  WARNSTAT returns the Status of WARNING
%
% Matlab R<13   WARNING
% Matlab R>12   WARNING for Identifier ALL
%

ww = warning;

if isstruct(ww)   % New Matlab R>12 Syntax
   try
      id = strcmp({ww.identifier},'all');
      if any(id)
         id = find(id);
         ww = ww(id(1)).state;
      else
         ww = '';
      end
   catch
      ww = '';
   end
elseif ~chkstr(ww)
   ww = '';
end

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = chkstr(str,opt)


% CHKSTR  Checks Input for String
%
%  ok = chkstr(str,Option)
%
%  Option ~= 0 ==> only true for nonempty Strings
%
%   default:   Option == 0  (true for empty Strings)
%

 
if nargin < 2
   opt = 0;
end

ok = ( strcmp( class(str) , 'char' )      & ...
       ( prod(size(str)) == size(str,2) ) & ...
       ( isequal(opt,0) | ~isempty(str) )         );



%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  c = rgb

% RGB  returns X-Windows RGB-Colors 
%
% RGB is a 3-Column CellArray: 
%             { ColorName   HexString   DecColors }
%  example:     'white'      '#FFFFFF'   [ 255 255 255 ] 

 c = { ... 

  'snow'                    '#FFFAFA'    [ 255 250 250 ]
  'ghost white'             '#F8F8FF'    [ 248 248 255 ]
  'GhostWhite'              '#F8F8FF'    [ 248 248 255 ]
  'white smoke'             '#F5F5F5'    [ 245 245 245 ]
  'WhiteSmoke'              '#F5F5F5'    [ 245 245 245 ]
  'gainsboro'               '#DCDCDC'    [ 220 220 220 ]
  'floral white'            '#FFFAF0'    [ 255 250 240 ]
  'FloralWhite'             '#FFFAF0'    [ 255 250 240 ]
  'old lace'                '#FDF5E6'    [ 253 245 230 ]
  'OldLace'                 '#FDF5E6'    [ 253 245 230 ]
  'linen'                   '#FAF0E6'    [ 250 240 230 ]
  'antique white'           '#FAEBD7'    [ 250 235 215 ]
  'AntiqueWhite'            '#FAEBD7'    [ 250 235 215 ]
  'papaya whip'             '#FFEFD5'    [ 255 239 213 ]
  'PapayaWhip'              '#FFEFD5'    [ 255 239 213 ]
  'blanched almond'         '#FFEBCD'    [ 255 235 205 ]
  'BlanchedAlmond'          '#FFEBCD'    [ 255 235 205 ]
  'bisque'                  '#FFE4C4'    [ 255 228 196 ]
  'peach puff'              '#FFDAB9'    [ 255 218 185 ]
  'PeachPuff'               '#FFDAB9'    [ 255 218 185 ]
  'navajo white'            '#FFDEAD'    [ 255 222 173 ]
  'NavajoWhite'             '#FFDEAD'    [ 255 222 173 ]
  'moccasin'                '#FFE4B5'    [ 255 228 181 ]
  'cornsilk'                '#FFF8DC'    [ 255 248 220 ]
  'ivory'                   '#FFFFF0'    [ 255 255 240 ]
  'lemon chiffon'           '#FFFACD'    [ 255 250 205 ]
  'LemonChiffon'            '#FFFACD'    [ 255 250 205 ]
  'seashell'                '#FFF5EE'    [ 255 245 238 ]
  'honeydew'                '#F0FFF0'    [ 240 255 240 ]
  'mint cream'              '#F5FFFA'    [ 245 255 250 ]
  'MintCream'               '#F5FFFA'    [ 245 255 250 ]
  'azure'                   '#F0FFFF'    [ 240 255 255 ]
  'alice blue'              '#F0F8FF'    [ 240 248 255 ]
  'AliceBlue'               '#F0F8FF'    [ 240 248 255 ]
  'lavender'                '#E6E6FA'    [ 230 230 250 ]
  'lavender blush'          '#FFF0F5'    [ 255 240 245 ]
  'LavenderBlush'           '#FFF0F5'    [ 255 240 245 ]
  'misty rose'              '#FFE4E1'    [ 255 228 225 ]
  'MistyRose'               '#FFE4E1'    [ 255 228 225 ]
  'white'                   '#FFFFFF'    [ 255 255 255 ]
  'black'                   '#000000'    [   0   0   0 ]
  'dark slate gray'         '#2F4F4F'    [  47  79  79 ]
  'DarkSlateGray'           '#2F4F4F'    [  47  79  79 ]
  'dark slate grey'         '#2F4F4F'    [  47  79  79 ]
  'DarkSlateGrey'           '#2F4F4F'    [  47  79  79 ]
  'dim gray'                '#696969'    [ 105 105 105 ]
  'DimGray'                 '#696969'    [ 105 105 105 ]
  'dim grey'                '#696969'    [ 105 105 105 ]
  'DimGrey'                 '#696969'    [ 105 105 105 ]
  'slate gray'              '#708090'    [ 112 128 144 ]
  'SlateGray'               '#708090'    [ 112 128 144 ]
  'slate grey'              '#708090'    [ 112 128 144 ]
  'SlateGrey'               '#708090'    [ 112 128 144 ]
  'light slate gray'        '#778899'    [ 119 136 153 ]
  'LightSlateGray'          '#778899'    [ 119 136 153 ]
  'light slate grey'        '#778899'    [ 119 136 153 ]
  'LightSlateGrey'          '#778899'    [ 119 136 153 ]
  'gray'                    '#BEBEBE'    [ 190 190 190 ]
  'grey'                    '#BEBEBE'    [ 190 190 190 ]
  'light grey'              '#D3D3D3'    [ 211 211 211 ]
  'LightGrey'               '#D3D3D3'    [ 211 211 211 ]
  'light gray'              '#D3D3D3'    [ 211 211 211 ]
  'LightGray'               '#D3D3D3'    [ 211 211 211 ]
  'midnight blue'           '#191970'    [  25  25 112 ]
  'MidnightBlue'            '#191970'    [  25  25 112 ]
  'navy'                    '#000080'    [   0   0 128 ]
  'navy blue'               '#000080'    [   0   0 128 ]
  'NavyBlue'                '#000080'    [   0   0 128 ]
  'cornflower blue'         '#6495ED'    [ 100 149 237 ]
  'CornflowerBlue'          '#6495ED'    [ 100 149 237 ]
  'dark slate blue'         '#483D8B'    [  72  61 139 ]
  'DarkSlateBlue'           '#483D8B'    [  72  61 139 ]
  'slate blue'              '#6A5ACD'    [ 106  90 205 ]
  'SlateBlue'               '#6A5ACD'    [ 106  90 205 ]
  'medium slate blue'       '#7B68EE'    [ 123 104 238 ]
  'MediumSlateBlue'         '#7B68EE'    [ 123 104 238 ]
  'light slate blue'        '#8470FF'    [ 132 112 255 ]
  'LightSlateBlue'          '#8470FF'    [ 132 112 255 ]
  'medium blue'             '#0000CD'    [   0   0 205 ]
  'MediumBlue'              '#0000CD'    [   0   0 205 ]
  'royal blue'              '#4169E1'    [  65 105 225 ]
  'RoyalBlue'               '#4169E1'    [  65 105 225 ]
  'blue'                    '#0000FF'    [   0   0 255 ]
  'dodger blue'             '#1E90FF'    [  30 144 255 ]
  'DodgerBlue'              '#1E90FF'    [  30 144 255 ]
  'deep sky blue'           '#00BFFF'    [   0 191 255 ]
  'DeepSkyBlue'             '#00BFFF'    [   0 191 255 ]
  'sky blue'                '#87CEEB'    [ 135 206 235 ]
  'SkyBlue'                 '#87CEEB'    [ 135 206 235 ]
  'light sky blue'          '#87CEFA'    [ 135 206 250 ]
  'LightSkyBlue'            '#87CEFA'    [ 135 206 250 ]
  'steel blue'              '#4682B4'    [  70 130 180 ]
  'SteelBlue'               '#4682B4'    [  70 130 180 ]
  'light steel blue'        '#B0C4DE'    [ 176 196 222 ]
  'LightSteelBlue'          '#B0C4DE'    [ 176 196 222 ]
  'light blue'              '#ADD8E6'    [ 173 216 230 ]
  'LightBlue'               '#ADD8E6'    [ 173 216 230 ]
  'powder blue'             '#B0E0E6'    [ 176 224 230 ]
  'PowderBlue'              '#B0E0E6'    [ 176 224 230 ]
  'pale turquoise'          '#AFEEEE'    [ 175 238 238 ]
  'PaleTurquoise'           '#AFEEEE'    [ 175 238 238 ]
  'dark turquoise'          '#00CED1'    [   0 206 209 ]
  'DarkTurquoise'           '#00CED1'    [   0 206 209 ]
  'medium turquoise'        '#48D1CC'    [  72 209 204 ]
  'MediumTurquoise'         '#48D1CC'    [  72 209 204 ]
  'turquoise'               '#40E0D0'    [  64 224 208 ]
  'cyan'                    '#00FFFF'    [   0 255 255 ]
  'light cyan'              '#E0FFFF'    [ 224 255 255 ]
  'LightCyan'               '#E0FFFF'    [ 224 255 255 ]
  'cadet blue'              '#5F9EA0'    [  95 158 160 ]
  'CadetBlue'               '#5F9EA0'    [  95 158 160 ]
  'medium aquamarine'       '#66CDAA'    [ 102 205 170 ]
  'MediumAquamarine'        '#66CDAA'    [ 102 205 170 ]
  'aquamarine'              '#7FFFD4'    [ 127 255 212 ]
  'dark green'              '#006400'    [   0 100   0 ]
  'DarkGreen'               '#006400'    [   0 100   0 ]
  'dark olive green'        '#556B2F'    [  85 107  47 ]
  'DarkOliveGreen'          '#556B2F'    [  85 107  47 ]
  'dark sea green'          '#8FBC8F'    [ 143 188 143 ]
  'DarkSeaGreen'            '#8FBC8F'    [ 143 188 143 ]
  'sea green'               '#2E8B57'    [  46 139  87 ]
  'SeaGreen'                '#2E8B57'    [  46 139  87 ]
  'medium sea green'        '#3CB371'    [  60 179 113 ]
  'MediumSeaGreen'          '#3CB371'    [  60 179 113 ]
  'light sea green'         '#20B2AA'    [  32 178 170 ]
  'LightSeaGreen'           '#20B2AA'    [  32 178 170 ]
  'pale green'              '#98FB98'    [ 152 251 152 ]
  'PaleGreen'               '#98FB98'    [ 152 251 152 ]
  'spring green'            '#00FF7F'    [   0 255 127 ]
  'SpringGreen'             '#00FF7F'    [   0 255 127 ]
  'lawn green'              '#7CFC00'    [ 124 252   0 ]
  'LawnGreen'               '#7CFC00'    [ 124 252   0 ]
  'green'                   '#00FF00'    [   0 255   0 ]
  'chartreuse'              '#7FFF00'    [ 127 255   0 ]
  'medium spring green'     '#00FA9A'    [   0 250 154 ]
  'MediumSpringGreen'       '#00FA9A'    [   0 250 154 ]
  'green yellow'            '#ADFF2F'    [ 173 255  47 ]
  'GreenYellow'             '#ADFF2F'    [ 173 255  47 ]
  'lime green'              '#32CD32'    [  50 205  50 ]
  'LimeGreen'               '#32CD32'    [  50 205  50 ]
  'yellow green'            '#9ACD32'    [ 154 205  50 ]
  'YellowGreen'             '#9ACD32'    [ 154 205  50 ]
  'forest green'            '#228B22'    [  34 139  34 ]
  'ForestGreen'             '#228B22'    [  34 139  34 ]
  'olive drab'              '#6B8E23'    [ 107 142  35 ]
  'OliveDrab'               '#6B8E23'    [ 107 142  35 ]
  'dark khaki'              '#BDB76B'    [ 189 183 107 ]
  'DarkKhaki'               '#BDB76B'    [ 189 183 107 ]
  'khaki'                   '#F0E68C'    [ 240 230 140 ]
  'pale goldenrod'          '#EEE8AA'    [ 238 232 170 ]
  'PaleGoldenrod'           '#EEE8AA'    [ 238 232 170 ]
  'light goldenrod yellow'  '#FAFAD2'    [ 250 250 210 ]
  'LightGoldenrodYellow'    '#FAFAD2'    [ 250 250 210 ]
  'light yellow'            '#FFFFE0'    [ 255 255 224 ]
  'LightYellow'             '#FFFFE0'    [ 255 255 224 ]
  'yellow'                  '#FFFF00'    [ 255 255   0 ]
  'gold'                    '#FFD700'    [ 255 215   0 ]
  'light goldenrod'         '#EEDD82'    [ 238 221 130 ]
  'LightGoldenrod'          '#EEDD82'    [ 238 221 130 ]
  'goldenrod'               '#DAA520'    [ 218 165  32 ]
  'dark goldenrod'          '#B8860B'    [ 184 134  11 ]
  'DarkGoldenrod'           '#B8860B'    [ 184 134  11 ]
  'rosy brown'              '#BC8F8F'    [ 188 143 143 ]
  'RosyBrown'               '#BC8F8F'    [ 188 143 143 ]
  'indian red'              '#CD5C5C'    [ 205  92  92 ]
  'IndianRed'               '#CD5C5C'    [ 205  92  92 ]
  'saddle brown'            '#8B4513'    [ 139  69  19 ]
  'SaddleBrown'             '#8B4513'    [ 139  69  19 ]
  'sienna'                  '#A0522D'    [ 160  82  45 ]
  'peru'                    '#CD853F'    [ 205 133  63 ]
  'burlywood'               '#DEB887'    [ 222 184 135 ]
  'beige'                   '#F5F5DC'    [ 245 245 220 ]
  'wheat'                   '#F5DEB3'    [ 245 222 179 ]
  'sandy brown'             '#F4A460'    [ 244 164  96 ]
  'SandyBrown'              '#F4A460'    [ 244 164  96 ]
  'tan'                     '#D2B48C'    [ 210 180 140 ]
  'chocolate'               '#D2691E'    [ 210 105  30 ]
  'firebrick'               '#B22222'    [ 178  34  34 ]
  'brown'                   '#A52A2A'    [ 165  42  42 ]
  'dark salmon'             '#E9967A'    [ 233 150 122 ]
  'DarkSalmon'              '#E9967A'    [ 233 150 122 ]
  'salmon'                  '#FA8072'    [ 250 128 114 ]
  'light salmon'            '#FFA07A'    [ 255 160 122 ]
  'LightSalmon'             '#FFA07A'    [ 255 160 122 ]
  'orange'                  '#FFA500'    [ 255 165   0 ]
  'dark orange'             '#FF8C00'    [ 255 140   0 ]
  'DarkOrange'              '#FF8C00'    [ 255 140   0 ]
  'coral'                   '#FF7F50'    [ 255 127  80 ]
  'light coral'             '#F08080'    [ 240 128 128 ]
  'LightCoral'              '#F08080'    [ 240 128 128 ]
  'tomato'                  '#FF6347'    [ 255  99  71 ]
  'orange red'              '#FF4500'    [ 255  69   0 ]
  'OrangeRed'               '#FF4500'    [ 255  69   0 ]
  'red'                     '#FF0000'    [ 255   0   0 ]
  'hot pink'                '#FF69B4'    [ 255 105 180 ]
  'HotPink'                 '#FF69B4'    [ 255 105 180 ]
  'deep pink'               '#FF1493'    [ 255  20 147 ]
  'DeepPink'                '#FF1493'    [ 255  20 147 ]
  'pink'                    '#FFC0CB'    [ 255 192 203 ]
  'light pink'              '#FFB6C1'    [ 255 182 193 ]
  'LightPink'               '#FFB6C1'    [ 255 182 193 ]
  'pale violet red'         '#DB7093'    [ 219 112 147 ]
  'PaleVioletRed'           '#DB7093'    [ 219 112 147 ]
  'maroon'                  '#B03060'    [ 176  48  96 ]
  'medium violet red'       '#C71585'    [ 199  21 133 ]
  'MediumVioletRed'         '#C71585'    [ 199  21 133 ]
  'violet red'              '#D02090'    [ 208  32 144 ]
  'VioletRed'               '#D02090'    [ 208  32 144 ]
  'magenta'                 '#FF00FF'    [ 255   0 255 ]
  'violet'                  '#EE82EE'    [ 238 130 238 ]
  'plum'                    '#DDA0DD'    [ 221 160 221 ]
  'orchid'                  '#DA70D6'    [ 218 112 214 ]
  'medium orchid'           '#BA55D3'    [ 186  85 211 ]
  'MediumOrchid'            '#BA55D3'    [ 186  85 211 ]
  'dark orchid'             '#9932CC'    [ 153  50 204 ]
  'DarkOrchid'              '#9932CC'    [ 153  50 204 ]
  'dark violet'             '#9400D3'    [ 148   0 211 ]
  'DarkViolet'              '#9400D3'    [ 148   0 211 ]
  'blue violet'             '#8A2BE2'    [ 138  43 226 ]
  'BlueViolet'              '#8A2BE2'    [ 138  43 226 ]
  'purple'                  '#A020F0'    [ 160  32 240 ]
  'medium purple'           '#9370DB'    [ 147 112 219 ]
  'MediumPurple'            '#9370DB'    [ 147 112 219 ]
  'thistle'                 '#D8BFD8'    [ 216 191 216 ]
  'snow1'                   '#FFFAFA'    [ 255 250 250 ]
  'snow2'                   '#EEE9E9'    [ 238 233 233 ]
  'snow3'                   '#CDC9C9'    [ 205 201 201 ]
  'snow4'                   '#8B8989'    [ 139 137 137 ]
  'seashell1'               '#FFF5EE'    [ 255 245 238 ]
  'seashell2'               '#EEE5DE'    [ 238 229 222 ]
  'seashell3'               '#CDC5BF'    [ 205 197 191 ]
  'seashell4'               '#8B8682'    [ 139 134 130 ]
  'AntiqueWhite1'           '#FFEFDB'    [ 255 239 219 ]
  'AntiqueWhite2'           '#EEDFCC'    [ 238 223 204 ]
  'AntiqueWhite3'           '#CDC0B0'    [ 205 192 176 ]
  'AntiqueWhite4'           '#8B8378'    [ 139 131 120 ]
  'bisque1'                 '#FFE4C4'    [ 255 228 196 ]
  'bisque2'                 '#EED5B7'    [ 238 213 183 ]
  'bisque3'                 '#CDB79E'    [ 205 183 158 ]
  'bisque4'                 '#8B7D6B'    [ 139 125 107 ]
  'PeachPuff1'              '#FFDAB9'    [ 255 218 185 ]
  'PeachPuff2'              '#EECBAD'    [ 238 203 173 ]
  'PeachPuff3'              '#CDAF95'    [ 205 175 149 ]
  'PeachPuff4'              '#8B7765'    [ 139 119 101 ]
  'NavajoWhite1'            '#FFDEAD'    [ 255 222 173 ]
  'NavajoWhite2'            '#EECFA1'    [ 238 207 161 ]
  'NavajoWhite3'            '#CDB38B'    [ 205 179 139 ]
  'NavajoWhite4'            '#8B795E'    [ 139 121  94 ]
  'LemonChiffon1'           '#FFFACD'    [ 255 250 205 ]
  'LemonChiffon2'           '#EEE9BF'    [ 238 233 191 ]
  'LemonChiffon3'           '#CDC9A5'    [ 205 201 165 ]
  'LemonChiffon4'           '#8B8970'    [ 139 137 112 ]
  'cornsilk1'               '#FFF8DC'    [ 255 248 220 ]
  'cornsilk2'               '#EEE8CD'    [ 238 232 205 ]
  'cornsilk3'               '#CDC8B1'    [ 205 200 177 ]
  'cornsilk4'               '#8B8878'    [ 139 136 120 ]
  'ivory1'                  '#FFFFF0'    [ 255 255 240 ]
  'ivory2'                  '#EEEEE0'    [ 238 238 224 ]
  'ivory3'                  '#CDCDC1'    [ 205 205 193 ]
  'ivory4'                  '#8B8B83'    [ 139 139 131 ]
  'honeydew1'               '#F0FFF0'    [ 240 255 240 ]
  'honeydew2'               '#E0EEE0'    [ 224 238 224 ]
  'honeydew3'               '#C1CDC1'    [ 193 205 193 ]
  'honeydew4'               '#838B83'    [ 131 139 131 ]
  'LavenderBlush1'          '#FFF0F5'    [ 255 240 245 ]
  'LavenderBlush2'          '#EEE0E5'    [ 238 224 229 ]
  'LavenderBlush3'          '#CDC1C5'    [ 205 193 197 ]
  'LavenderBlush4'          '#8B8386'    [ 139 131 134 ]
  'MistyRose1'              '#FFE4E1'    [ 255 228 225 ]
  'MistyRose2'              '#EED5D2'    [ 238 213 210 ]
  'MistyRose3'              '#CDB7B5'    [ 205 183 181 ]
  'MistyRose4'              '#8B7D7B'    [ 139 125 123 ]
  'azure1'                  '#F0FFFF'    [ 240 255 255 ]
  'azure2'                  '#E0EEEE'    [ 224 238 238 ]
  'azure3'                  '#C1CDCD'    [ 193 205 205 ]
  'azure4'                  '#838B8B'    [ 131 139 139 ]
  'SlateBlue1'              '#836FFF'    [ 131 111 255 ]
  'SlateBlue2'              '#7A67EE'    [ 122 103 238 ]
  'SlateBlue3'              '#6959CD'    [ 105  89 205 ]
  'SlateBlue4'              '#473C8B'    [  71  60 139 ]
  'RoyalBlue1'              '#4876FF'    [  72 118 255 ]
  'RoyalBlue2'              '#436EEE'    [  67 110 238 ]
  'RoyalBlue3'              '#3A5FCD'    [  58  95 205 ]
  'RoyalBlue4'              '#27408B'    [  39  64 139 ]
  'blue1'                   '#0000FF'    [   0   0 255 ]
  'blue2'                   '#0000EE'    [   0   0 238 ]
  'blue3'                   '#0000CD'    [   0   0 205 ]
  'blue4'                   '#00008B'    [   0   0 139 ]
  'DodgerBlue1'             '#1E90FF'    [  30 144 255 ]
  'DodgerBlue2'             '#1C86EE'    [  28 134 238 ]
  'DodgerBlue3'             '#1874CD'    [  24 116 205 ]
  'DodgerBlue4'             '#104E8B'    [  16  78 139 ]
  'SteelBlue1'              '#63B8FF'    [  99 184 255 ]
  'SteelBlue2'              '#5CACEE'    [  92 172 238 ]
  'SteelBlue3'              '#4F94CD'    [  79 148 205 ]
  'SteelBlue4'              '#36648B'    [  54 100 139 ]
  'DeepSkyBlue1'            '#00BFFF'    [   0 191 255 ]
  'DeepSkyBlue2'            '#00B2EE'    [   0 178 238 ]
  'DeepSkyBlue3'            '#009ACD'    [   0 154 205 ]
  'DeepSkyBlue4'            '#00688B'    [   0 104 139 ]
  'SkyBlue1'                '#87CEFF'    [ 135 206 255 ]
  'SkyBlue2'                '#7EC0EE'    [ 126 192 238 ]
  'SkyBlue3'                '#6CA6CD'    [ 108 166 205 ]
  'SkyBlue4'                '#4A708B'    [  74 112 139 ]
  'LightSkyBlue1'           '#B0E2FF'    [ 176 226 255 ]
  'LightSkyBlue2'           '#A4D3EE'    [ 164 211 238 ]
  'LightSkyBlue3'           '#8DB6CD'    [ 141 182 205 ]
  'LightSkyBlue4'           '#607B8B'    [  96 123 139 ]
  'SlateGray1'              '#C6E2FF'    [ 198 226 255 ]
  'SlateGray2'              '#B9D3EE'    [ 185 211 238 ]
  'SlateGray3'              '#9FB6CD'    [ 159 182 205 ]
  'SlateGray4'              '#6C7B8B'    [ 108 123 139 ]
  'LightSteelBlue1'         '#CAE1FF'    [ 202 225 255 ]
  'LightSteelBlue2'         '#BCD2EE'    [ 188 210 238 ]
  'LightSteelBlue3'         '#A2B5CD'    [ 162 181 205 ]
  'LightSteelBlue4'         '#6E7B8B'    [ 110 123 139 ]
  'LightBlue1'              '#BFEFFF'    [ 191 239 255 ]
  'LightBlue2'              '#B2DFEE'    [ 178 223 238 ]
  'LightBlue3'              '#9AC0CD'    [ 154 192 205 ]
  'LightBlue4'              '#68838B'    [ 104 131 139 ]
  'LightCyan1'              '#E0FFFF'    [ 224 255 255 ]
  'LightCyan2'              '#D1EEEE'    [ 209 238 238 ]
  'LightCyan3'              '#B4CDCD'    [ 180 205 205 ]
  'LightCyan4'              '#7A8B8B'    [ 122 139 139 ]
  'PaleTurquoise1'          '#BBFFFF'    [ 187 255 255 ]
  'PaleTurquoise2'          '#AEEEEE'    [ 174 238 238 ]
  'PaleTurquoise3'          '#96CDCD'    [ 150 205 205 ]
  'PaleTurquoise4'          '#668B8B'    [ 102 139 139 ]
  'CadetBlue1'              '#98F5FF'    [ 152 245 255 ]
  'CadetBlue2'              '#8EE5EE'    [ 142 229 238 ]
  'CadetBlue3'              '#7AC5CD'    [ 122 197 205 ]
  'CadetBlue4'              '#53868B'    [  83 134 139 ]
  'turquoise1'              '#00F5FF'    [   0 245 255 ]
  'turquoise2'              '#00E5EE'    [   0 229 238 ]
  'turquoise3'              '#00C5CD'    [   0 197 205 ]
  'turquoise4'              '#00868B'    [   0 134 139 ]
  'cyan1'                   '#00FFFF'    [   0 255 255 ]
  'cyan2'                   '#00EEEE'    [   0 238 238 ]
  'cyan3'                   '#00CDCD'    [   0 205 205 ]
  'cyan4'                   '#008B8B'    [   0 139 139 ]
  'DarkSlateGray1'          '#97FFFF'    [ 151 255 255 ]
  'DarkSlateGray2'          '#8DEEEE'    [ 141 238 238 ]
  'DarkSlateGray3'          '#79CDCD'    [ 121 205 205 ]
  'DarkSlateGray4'          '#528B8B'    [  82 139 139 ]
  'aquamarine1'             '#7FFFD4'    [ 127 255 212 ]
  'aquamarine2'             '#76EEC6'    [ 118 238 198 ]
  'aquamarine3'             '#66CDAA'    [ 102 205 170 ]
  'aquamarine4'             '#458B74'    [  69 139 116 ]
  'DarkSeaGreen1'           '#C1FFC1'    [ 193 255 193 ]
  'DarkSeaGreen2'           '#B4EEB4'    [ 180 238 180 ]
  'DarkSeaGreen3'           '#9BCD9B'    [ 155 205 155 ]
  'DarkSeaGreen4'           '#698B69'    [ 105 139 105 ]
  'SeaGreen1'               '#54FF9F'    [  84 255 159 ]
  'SeaGreen2'               '#4EEE94'    [  78 238 148 ]
  'SeaGreen3'               '#43CD80'    [  67 205 128 ]
  'SeaGreen4'               '#2E8B57'    [  46 139  87 ]
  'PaleGreen1'              '#9AFF9A'    [ 154 255 154 ]
  'PaleGreen2'              '#90EE90'    [ 144 238 144 ]
  'PaleGreen3'              '#7CCD7C'    [ 124 205 124 ]
  'PaleGreen4'              '#548B54'    [  84 139  84 ]
  'SpringGreen1'            '#00FF7F'    [   0 255 127 ]
  'SpringGreen2'            '#00EE76'    [   0 238 118 ]
  'SpringGreen3'            '#00CD66'    [   0 205 102 ]
  'SpringGreen4'            '#008B45'    [   0 139  69 ]
  'green1'                  '#00FF00'    [   0 255   0 ]
  'green2'                  '#00EE00'    [   0 238   0 ]
  'green3'                  '#00CD00'    [   0 205   0 ]
  'green4'                  '#008B00'    [   0 139   0 ]
  'chartreuse1'             '#7FFF00'    [ 127 255   0 ]
  'chartreuse2'             '#76EE00'    [ 118 238   0 ]
  'chartreuse3'             '#66CD00'    [ 102 205   0 ]
  'chartreuse4'             '#458B00'    [  69 139   0 ]
  'OliveDrab1'              '#C0FF3E'    [ 192 255  62 ]
  'OliveDrab2'              '#B3EE3A'    [ 179 238  58 ]
  'OliveDrab3'              '#9ACD32'    [ 154 205  50 ]
  'OliveDrab4'              '#698B22'    [ 105 139  34 ]
  'DarkOliveGreen1'         '#CAFF70'    [ 202 255 112 ]
  'DarkOliveGreen2'         '#BCEE68'    [ 188 238 104 ]
  'DarkOliveGreen3'         '#A2CD5A'    [ 162 205  90 ]
  'DarkOliveGreen4'         '#6E8B3D'    [ 110 139  61 ]
  'khaki1'                  '#FFF68F'    [ 255 246 143 ]
  'khaki2'                  '#EEE685'    [ 238 230 133 ]
  'khaki3'                  '#CDC673'    [ 205 198 115 ]
  'khaki4'                  '#8B864E'    [ 139 134  78 ]
  'LightGoldenrod1'         '#FFEC8B'    [ 255 236 139 ]
  'LightGoldenrod2'         '#EEDC82'    [ 238 220 130 ]
  'LightGoldenrod3'         '#CDBE70'    [ 205 190 112 ]
  'LightGoldenrod4'         '#8B814C'    [ 139 129  76 ]
  'LightYellow1'            '#FFFFE0'    [ 255 255 224 ]
  'LightYellow2'            '#EEEED1'    [ 238 238 209 ]
  'LightYellow3'            '#CDCDB4'    [ 205 205 180 ]
  'LightYellow4'            '#8B8B7A'    [ 139 139 122 ]
  'yellow1'                 '#FFFF00'    [ 255 255   0 ]
  'yellow2'                 '#EEEE00'    [ 238 238   0 ]
  'yellow3'                 '#CDCD00'    [ 205 205   0 ]
  'yellow4'                 '#8B8B00'    [ 139 139   0 ]
  'gold1'                   '#FFD700'    [ 255 215   0 ]
  'gold2'                   '#EEC900'    [ 238 201   0 ]
  'gold3'                   '#CDAD00'    [ 205 173   0 ]
  'gold4'                   '#8B7500'    [ 139 117   0 ]
  'goldenrod1'              '#FFC125'    [ 255 193  37 ]
  'goldenrod2'              '#EEB422'    [ 238 180  34 ]
  'goldenrod3'              '#CD9B1D'    [ 205 155  29 ]
  'goldenrod4'              '#8B6914'    [ 139 105  20 ]
  'DarkGoldenrod1'          '#FFB90F'    [ 255 185  15 ]
  'DarkGoldenrod2'          '#EEAD0E'    [ 238 173  14 ]
  'DarkGoldenrod3'          '#CD950C'    [ 205 149  12 ]
  'DarkGoldenrod4'          '#8B6508'    [ 139 101   8 ]
  'RosyBrown1'              '#FFC1C1'    [ 255 193 193 ]
  'RosyBrown2'              '#EEB4B4'    [ 238 180 180 ]
  'RosyBrown3'              '#CD9B9B'    [ 205 155 155 ]
  'RosyBrown4'              '#8B6969'    [ 139 105 105 ]
  'IndianRed1'              '#FF6A6A'    [ 255 106 106 ]
  'IndianRed2'              '#EE6363'    [ 238  99  99 ]
  'IndianRed3'              '#CD5555'    [ 205  85  85 ]
  'IndianRed4'              '#8B3A3A'    [ 139  58  58 ]
  'sienna1'                 '#FF8247'    [ 255 130  71 ]
  'sienna2'                 '#EE7942'    [ 238 121  66 ]
  'sienna3'                 '#CD6839'    [ 205 104  57 ]
  'sienna4'                 '#8B4726'    [ 139  71  38 ]
  'burlywood1'              '#FFD39B'    [ 255 211 155 ]
  'burlywood2'              '#EEC591'    [ 238 197 145 ]
  'burlywood3'              '#CDAA7D'    [ 205 170 125 ]
  'burlywood4'              '#8B7355'    [ 139 115  85 ]
  'wheat1'                  '#FFE7BA'    [ 255 231 186 ]
  'wheat2'                  '#EED8AE'    [ 238 216 174 ]
  'wheat3'                  '#CDBA96'    [ 205 186 150 ]
  'wheat4'                  '#8B7E66'    [ 139 126 102 ]
  'tan1'                    '#FFA54F'    [ 255 165  79 ]
  'tan2'                    '#EE9A49'    [ 238 154  73 ]
  'tan3'                    '#CD853F'    [ 205 133  63 ]
  'tan4'                    '#8B5A2B'    [ 139  90  43 ]
  'chocolate1'              '#FF7F24'    [ 255 127  36 ]
  'chocolate2'              '#EE7621'    [ 238 118  33 ]
  'chocolate3'              '#CD661D'    [ 205 102  29 ]
  'chocolate4'              '#8B4513'    [ 139  69  19 ]
  'firebrick1'              '#FF3030'    [ 255  48  48 ]
  'firebrick2'              '#EE2C2C'    [ 238  44  44 ]
  'firebrick3'              '#CD2626'    [ 205  38  38 ]
  'firebrick4'              '#8B1A1A'    [ 139  26  26 ]
  'brown1'                  '#FF4040'    [ 255  64  64 ]
  'brown2'                  '#EE3B3B'    [ 238  59  59 ]
  'brown3'                  '#CD3333'    [ 205  51  51 ]
  'brown4'                  '#8B2323'    [ 139  35  35 ]
  'salmon1'                 '#FF8C69'    [ 255 140 105 ]
  'salmon2'                 '#EE8262'    [ 238 130  98 ]
  'salmon3'                 '#CD7054'    [ 205 112  84 ]
  'salmon4'                 '#8B4C39'    [ 139  76  57 ]
  'LightSalmon1'            '#FFA07A'    [ 255 160 122 ]
  'LightSalmon2'            '#EE9572'    [ 238 149 114 ]
  'LightSalmon3'            '#CD8162'    [ 205 129  98 ]
  'LightSalmon4'            '#8B5742'    [ 139  87  66 ]
  'orange1'                 '#FFA500'    [ 255 165   0 ]
  'orange2'                 '#EE9A00'    [ 238 154   0 ]
  'orange3'                 '#CD8500'    [ 205 133   0 ]
  'orange4'                 '#8B5A00'    [ 139  90   0 ]
  'DarkOrange1'             '#FF7F00'    [ 255 127   0 ]
  'DarkOrange2'             '#EE7600'    [ 238 118   0 ]
  'DarkOrange3'             '#CD6600'    [ 205 102   0 ]
  'DarkOrange4'             '#8B4500'    [ 139  69   0 ]
  'coral1'                  '#FF7256'    [ 255 114  86 ]
  'coral2'                  '#EE6A50'    [ 238 106  80 ]
  'coral3'                  '#CD5B45'    [ 205  91  69 ]
  'coral4'                  '#8B3E2F'    [ 139  62  47 ]
  'tomato1'                 '#FF6347'    [ 255  99  71 ]
  'tomato2'                 '#EE5C42'    [ 238  92  66 ]
  'tomato3'                 '#CD4F39'    [ 205  79  57 ]
  'tomato4'                 '#8B3626'    [ 139  54  38 ]
  'OrangeRed1'              '#FF4500'    [ 255  69   0 ]
  'OrangeRed2'              '#EE4000'    [ 238  64   0 ]
  'OrangeRed3'              '#CD3700'    [ 205  55   0 ]
  'OrangeRed4'              '#8B2500'    [ 139  37   0 ]
  'red1'                    '#FF0000'    [ 255   0   0 ]
  'red2'                    '#EE0000'    [ 238   0   0 ]
  'red3'                    '#CD0000'    [ 205   0   0 ]
  'red4'                    '#8B0000'    [ 139   0   0 ]
  'DeepPink1'               '#FF1493'    [ 255  20 147 ]
  'DeepPink2'               '#EE1289'    [ 238  18 137 ]
  'DeepPink3'               '#CD1076'    [ 205  16 118 ]
  'DeepPink4'               '#8B0A50'    [ 139  10  80 ]
  'HotPink1'                '#FF6EB4'    [ 255 110 180 ]
  'HotPink2'                '#EE6AA7'    [ 238 106 167 ]
  'HotPink3'                '#CD6090'    [ 205  96 144 ]
  'HotPink4'                '#8B3A62'    [ 139  58  98 ]
  'pink1'                   '#FFB5C5'    [ 255 181 197 ]
  'pink2'                   '#EEA9B8'    [ 238 169 184 ]
  'pink3'                   '#CD919E'    [ 205 145 158 ]
  'pink4'                   '#8B636C'    [ 139  99 108 ]
  'LightPink1'              '#FFAEB9'    [ 255 174 185 ]
  'LightPink2'              '#EEA2AD'    [ 238 162 173 ]
  'LightPink3'              '#CD8C95'    [ 205 140 149 ]
  'LightPink4'              '#8B5F65'    [ 139  95 101 ]
  'PaleVioletRed1'          '#FF82AB'    [ 255 130 171 ]
  'PaleVioletRed2'          '#EE799F'    [ 238 121 159 ]
  'PaleVioletRed3'          '#CD6889'    [ 205 104 137 ]
  'PaleVioletRed4'          '#8B475D'    [ 139  71  93 ]
  'maroon1'                 '#FF34B3'    [ 255  52 179 ]
  'maroon2'                 '#EE30A7'    [ 238  48 167 ]
  'maroon3'                 '#CD2990'    [ 205  41 144 ]
  'maroon4'                 '#8B1C62'    [ 139  28  98 ]
  'VioletRed1'              '#FF3E96'    [ 255  62 150 ]
  'VioletRed2'              '#EE3A8C'    [ 238  58 140 ]
  'VioletRed3'              '#CD3278'    [ 205  50 120 ]
  'VioletRed4'              '#8B2252'    [ 139  34  82 ]
  'magenta1'                '#FF00FF'    [ 255   0 255 ]
  'magenta2'                '#EE00EE'    [ 238   0 238 ]
  'magenta3'                '#CD00CD'    [ 205   0 205 ]
  'magenta4'                '#8B008B'    [ 139   0 139 ]
  'orchid1'                 '#FF83FA'    [ 255 131 250 ]
  'orchid2'                 '#EE7AE9'    [ 238 122 233 ]
  'orchid3'                 '#CD69C9'    [ 205 105 201 ]
  'orchid4'                 '#8B4789'    [ 139  71 137 ]
  'plum1'                   '#FFBBFF'    [ 255 187 255 ]
  'plum2'                   '#EEAEEE'    [ 238 174 238 ]
  'plum3'                   '#CD96CD'    [ 205 150 205 ]
  'plum4'                   '#8B668B'    [ 139 102 139 ]
  'MediumOrchid1'           '#E066FF'    [ 224 102 255 ]
  'MediumOrchid2'           '#D15FEE'    [ 209  95 238 ]
  'MediumOrchid3'           '#B452CD'    [ 180  82 205 ]
  'MediumOrchid4'           '#7A378B'    [ 122  55 139 ]
  'DarkOrchid1'             '#BF3EFF'    [ 191  62 255 ]
  'DarkOrchid2'             '#B23AEE'    [ 178  58 238 ]
  'DarkOrchid3'             '#9A32CD'    [ 154  50 205 ]
  'DarkOrchid4'             '#68228B'    [ 104  34 139 ]
  'purple1'                 '#9B30FF'    [ 155  48 255 ]
  'purple2'                 '#912CEE'    [ 145  44 238 ]
  'purple3'                 '#7D26CD'    [ 125  38 205 ]
  'purple4'                 '#551A8B'    [  85  26 139 ]
  'MediumPurple1'           '#AB82FF'    [ 171 130 255 ]
  'MediumPurple2'           '#9F79EE'    [ 159 121 238 ]
  'MediumPurple3'           '#8968CD'    [ 137 104 205 ]
  'MediumPurple4'           '#5D478B'    [  93  71 139 ]
  'thistle1'                '#FFE1FF'    [ 255 225 255 ]
  'thistle2'                '#EED2EE'    [ 238 210 238 ]
  'thistle3'                '#CDB5CD'    [ 205 181 205 ]
  'thistle4'                '#8B7B8B'    [ 139 123 139 ]
  'gray0'                   '#000000'    [   0   0   0 ]
  'grey0'                   '#000000'    [   0   0   0 ]
  'gray1'                   '#030303'    [   3   3   3 ]
  'grey1'                   '#030303'    [   3   3   3 ]
  'gray2'                   '#050505'    [   5   5   5 ]
  'grey2'                   '#050505'    [   5   5   5 ]
  'gray3'                   '#080808'    [   8   8   8 ]
  'grey3'                   '#080808'    [   8   8   8 ]
  'gray4'                   '#0A0A0A'    [  10  10  10 ]
  'grey4'                   '#0A0A0A'    [  10  10  10 ]
  'gray5'                   '#0D0D0D'    [  13  13  13 ]
  'grey5'                   '#0D0D0D'    [  13  13  13 ]
  'gray6'                   '#0F0F0F'    [  15  15  15 ]
  'grey6'                   '#0F0F0F'    [  15  15  15 ]
  'gray7'                   '#121212'    [  18  18  18 ]
  'grey7'                   '#121212'    [  18  18  18 ]
  'gray8'                   '#141414'    [  20  20  20 ]
  'grey8'                   '#141414'    [  20  20  20 ]
  'gray9'                   '#171717'    [  23  23  23 ]
  'grey9'                   '#171717'    [  23  23  23 ]
  'gray10'                  '#1A1A1A'    [  26  26  26 ]
  'grey10'                  '#1A1A1A'    [  26  26  26 ]
  'gray11'                  '#1C1C1C'    [  28  28  28 ]
  'grey11'                  '#1C1C1C'    [  28  28  28 ]
  'gray12'                  '#1F1F1F'    [  31  31  31 ]
  'grey12'                  '#1F1F1F'    [  31  31  31 ]
  'gray13'                  '#212121'    [  33  33  33 ]
  'grey13'                  '#212121'    [  33  33  33 ]
  'gray14'                  '#242424'    [  36  36  36 ]
  'grey14'                  '#242424'    [  36  36  36 ]
  'gray15'                  '#262626'    [  38  38  38 ]
  'grey15'                  '#262626'    [  38  38  38 ]
  'gray16'                  '#292929'    [  41  41  41 ]
  'grey16'                  '#292929'    [  41  41  41 ]
  'gray17'                  '#2B2B2B'    [  43  43  43 ]
  'grey17'                  '#2B2B2B'    [  43  43  43 ]
  'gray18'                  '#2E2E2E'    [  46  46  46 ]
  'grey18'                  '#2E2E2E'    [  46  46  46 ]
  'gray19'                  '#303030'    [  48  48  48 ]
  'grey19'                  '#303030'    [  48  48  48 ]
  'gray20'                  '#333333'    [  51  51  51 ]
  'grey20'                  '#333333'    [  51  51  51 ]
  'gray21'                  '#363636'    [  54  54  54 ]
  'grey21'                  '#363636'    [  54  54  54 ]
  'gray22'                  '#383838'    [  56  56  56 ]
  'grey22'                  '#383838'    [  56  56  56 ]
  'gray23'                  '#3B3B3B'    [  59  59  59 ]
  'grey23'                  '#3B3B3B'    [  59  59  59 ]
  'gray24'                  '#3D3D3D'    [  61  61  61 ]
  'grey24'                  '#3D3D3D'    [  61  61  61 ]
  'gray25'                  '#404040'    [  64  64  64 ]
  'grey25'                  '#404040'    [  64  64  64 ]
  'gray26'                  '#424242'    [  66  66  66 ]
  'grey26'                  '#424242'    [  66  66  66 ]
  'gray27'                  '#454545'    [  69  69  69 ]
  'grey27'                  '#454545'    [  69  69  69 ]
  'gray28'                  '#474747'    [  71  71  71 ]
  'grey28'                  '#474747'    [  71  71  71 ]
  'gray29'                  '#4A4A4A'    [  74  74  74 ]
  'grey29'                  '#4A4A4A'    [  74  74  74 ]
  'gray30'                  '#4D4D4D'    [  77  77  77 ]
  'grey30'                  '#4D4D4D'    [  77  77  77 ]
  'gray31'                  '#4F4F4F'    [  79  79  79 ]
  'grey31'                  '#4F4F4F'    [  79  79  79 ]
  'gray32'                  '#525252'    [  82  82  82 ]
  'grey32'                  '#525252'    [  82  82  82 ]
  'gray33'                  '#545454'    [  84  84  84 ]
  'grey33'                  '#545454'    [  84  84  84 ]
  'gray34'                  '#575757'    [  87  87  87 ]
  'grey34'                  '#575757'    [  87  87  87 ]
  'gray35'                  '#595959'    [  89  89  89 ]
  'grey35'                  '#595959'    [  89  89  89 ]
  'gray36'                  '#5C5C5C'    [  92  92  92 ]
  'grey36'                  '#5C5C5C'    [  92  92  92 ]
  'gray37'                  '#5E5E5E'    [  94  94  94 ]
  'grey37'                  '#5E5E5E'    [  94  94  94 ]
  'gray38'                  '#616161'    [  97  97  97 ]
  'grey38'                  '#616161'    [  97  97  97 ]
  'gray39'                  '#636363'    [  99  99  99 ]
  'grey39'                  '#636363'    [  99  99  99 ]
  'gray40'                  '#666666'    [ 102 102 102 ]
  'grey40'                  '#666666'    [ 102 102 102 ]
  'gray41'                  '#696969'    [ 105 105 105 ]
  'grey41'                  '#696969'    [ 105 105 105 ]
  'gray42'                  '#6B6B6B'    [ 107 107 107 ]
  'grey42'                  '#6B6B6B'    [ 107 107 107 ]
  'gray43'                  '#6E6E6E'    [ 110 110 110 ]
  'grey43'                  '#6E6E6E'    [ 110 110 110 ]
  'gray44'                  '#707070'    [ 112 112 112 ]
  'grey44'                  '#707070'    [ 112 112 112 ]
  'gray45'                  '#737373'    [ 115 115 115 ]
  'grey45'                  '#737373'    [ 115 115 115 ]
  'gray46'                  '#757575'    [ 117 117 117 ]
  'grey46'                  '#757575'    [ 117 117 117 ]
  'gray47'                  '#787878'    [ 120 120 120 ]
  'grey47'                  '#787878'    [ 120 120 120 ]
  'gray48'                  '#7A7A7A'    [ 122 122 122 ]
  'grey48'                  '#7A7A7A'    [ 122 122 122 ]
  'gray49'                  '#7D7D7D'    [ 125 125 125 ]
  'grey49'                  '#7D7D7D'    [ 125 125 125 ]
  'gray50'                  '#7F7F7F'    [ 127 127 127 ]
  'grey50'                  '#7F7F7F'    [ 127 127 127 ]
  'gray51'                  '#828282'    [ 130 130 130 ]
  'grey51'                  '#828282'    [ 130 130 130 ]
  'gray52'                  '#858585'    [ 133 133 133 ]
  'grey52'                  '#858585'    [ 133 133 133 ]
  'gray53'                  '#878787'    [ 135 135 135 ]
  'grey53'                  '#878787'    [ 135 135 135 ]
  'gray54'                  '#8A8A8A'    [ 138 138 138 ]
  'grey54'                  '#8A8A8A'    [ 138 138 138 ]
  'gray55'                  '#8C8C8C'    [ 140 140 140 ]
  'grey55'                  '#8C8C8C'    [ 140 140 140 ]
  'gray56'                  '#8F8F8F'    [ 143 143 143 ]
  'grey56'                  '#8F8F8F'    [ 143 143 143 ]
  'gray57'                  '#919191'    [ 145 145 145 ]
  'grey57'                  '#919191'    [ 145 145 145 ]
  'gray58'                  '#949494'    [ 148 148 148 ]
  'grey58'                  '#949494'    [ 148 148 148 ]
  'gray59'                  '#969696'    [ 150 150 150 ]
  'grey59'                  '#969696'    [ 150 150 150 ]
  'gray60'                  '#999999'    [ 153 153 153 ]
  'grey60'                  '#999999'    [ 153 153 153 ]
  'gray61'                  '#9C9C9C'    [ 156 156 156 ]
  'grey61'                  '#9C9C9C'    [ 156 156 156 ]
  'gray62'                  '#9E9E9E'    [ 158 158 158 ]
  'grey62'                  '#9E9E9E'    [ 158 158 158 ]
  'gray63'                  '#A1A1A1'    [ 161 161 161 ]
  'grey63'                  '#A1A1A1'    [ 161 161 161 ]
  'gray64'                  '#A3A3A3'    [ 163 163 163 ]
  'grey64'                  '#A3A3A3'    [ 163 163 163 ]
  'gray65'                  '#A6A6A6'    [ 166 166 166 ]
  'grey65'                  '#A6A6A6'    [ 166 166 166 ]
  'gray66'                  '#A8A8A8'    [ 168 168 168 ]
  'grey66'                  '#A8A8A8'    [ 168 168 168 ]
  'gray67'                  '#ABABAB'    [ 171 171 171 ]
  'grey67'                  '#ABABAB'    [ 171 171 171 ]
  'gray68'                  '#ADADAD'    [ 173 173 173 ]
  'grey68'                  '#ADADAD'    [ 173 173 173 ]
  'gray69'                  '#B0B0B0'    [ 176 176 176 ]
  'grey69'                  '#B0B0B0'    [ 176 176 176 ]
  'gray70'                  '#B3B3B3'    [ 179 179 179 ]
  'grey70'                  '#B3B3B3'    [ 179 179 179 ]
  'gray71'                  '#B5B5B5'    [ 181 181 181 ]
  'grey71'                  '#B5B5B5'    [ 181 181 181 ]
  'gray72'                  '#B8B8B8'    [ 184 184 184 ]
  'grey72'                  '#B8B8B8'    [ 184 184 184 ]
  'gray73'                  '#BABABA'    [ 186 186 186 ]
  'grey73'                  '#BABABA'    [ 186 186 186 ]
  'gray74'                  '#BDBDBD'    [ 189 189 189 ]
  'grey74'                  '#BDBDBD'    [ 189 189 189 ]
  'gray75'                  '#BFBFBF'    [ 191 191 191 ]
  'grey75'                  '#BFBFBF'    [ 191 191 191 ]
  'gray76'                  '#C2C2C2'    [ 194 194 194 ]
  'grey76'                  '#C2C2C2'    [ 194 194 194 ]
  'gray77'                  '#C4C4C4'    [ 196 196 196 ]
  'grey77'                  '#C4C4C4'    [ 196 196 196 ]
  'gray78'                  '#C7C7C7'    [ 199 199 199 ]
  'grey78'                  '#C7C7C7'    [ 199 199 199 ]
  'gray79'                  '#C9C9C9'    [ 201 201 201 ]
  'grey79'                  '#C9C9C9'    [ 201 201 201 ]
  'gray80'                  '#CCCCCC'    [ 204 204 204 ]
  'grey80'                  '#CCCCCC'    [ 204 204 204 ]
  'gray81'                  '#CFCFCF'    [ 207 207 207 ]
  'grey81'                  '#CFCFCF'    [ 207 207 207 ]
  'gray82'                  '#D1D1D1'    [ 209 209 209 ]
  'grey82'                  '#D1D1D1'    [ 209 209 209 ]
  'gray83'                  '#D4D4D4'    [ 212 212 212 ]
  'grey83'                  '#D4D4D4'    [ 212 212 212 ]
  'gray84'                  '#D6D6D6'    [ 214 214 214 ]
  'grey84'                  '#D6D6D6'    [ 214 214 214 ]
  'gray85'                  '#D9D9D9'    [ 217 217 217 ]
  'grey85'                  '#D9D9D9'    [ 217 217 217 ]
  'gray86'                  '#DBDBDB'    [ 219 219 219 ]
  'grey86'                  '#DBDBDB'    [ 219 219 219 ]
  'gray87'                  '#DEDEDE'    [ 222 222 222 ]
  'grey87'                  '#DEDEDE'    [ 222 222 222 ]
  'gray88'                  '#E0E0E0'    [ 224 224 224 ]
  'grey88'                  '#E0E0E0'    [ 224 224 224 ]
  'gray89'                  '#E3E3E3'    [ 227 227 227 ]
  'grey89'                  '#E3E3E3'    [ 227 227 227 ]
  'gray90'                  '#E5E5E5'    [ 229 229 229 ]
  'grey90'                  '#E5E5E5'    [ 229 229 229 ]
  'gray91'                  '#E8E8E8'    [ 232 232 232 ]
  'grey91'                  '#E8E8E8'    [ 232 232 232 ]
  'gray92'                  '#EBEBEB'    [ 235 235 235 ]
  'grey92'                  '#EBEBEB'    [ 235 235 235 ]
  'gray93'                  '#EDEDED'    [ 237 237 237 ]
  'grey93'                  '#EDEDED'    [ 237 237 237 ]
  'gray94'                  '#F0F0F0'    [ 240 240 240 ]
  'grey94'                  '#F0F0F0'    [ 240 240 240 ]
  'gray95'                  '#F2F2F2'    [ 242 242 242 ]
  'grey95'                  '#F2F2F2'    [ 242 242 242 ]
  'gray96'                  '#F5F5F5'    [ 245 245 245 ]
  'grey96'                  '#F5F5F5'    [ 245 245 245 ]
  'gray97'                  '#F7F7F7'    [ 247 247 247 ]
  'grey97'                  '#F7F7F7'    [ 247 247 247 ]
  'gray98'                  '#FAFAFA'    [ 250 250 250 ]
  'grey98'                  '#FAFAFA'    [ 250 250 250 ]
  'gray99'                  '#FCFCFC'    [ 252 252 252 ]
  'grey99'                  '#FCFCFC'    [ 252 252 252 ]
  'gray100'                 '#FFFFFF'    [ 255 255 255 ]
  'grey100'                 '#FFFFFF'    [ 255 255 255 ]
  'dark grey'               '#A9A9A9'    [ 169 169 169 ]
  'DarkGrey'                '#A9A9A9'    [ 169 169 169 ]
  'dark gray'               '#A9A9A9'    [ 169 169 169 ]
  'DarkGray'                '#A9A9A9'    [ 169 169 169 ]
  'dark blue'               '#00008B'    [   0   0 139 ]
  'DarkBlue'                '#00008B'    [   0   0 139 ]
  'dark cyan'               '#008B8B'    [   0 139 139 ]
  'DarkCyan'                '#008B8B'    [   0 139 139 ]
  'dark magenta'            '#8B008B'    [ 139   0 139 ]
  'DarkMagenta'             '#8B008B'    [ 139   0 139 ]
  'dark red'                '#8B0000'    [ 139   0   0 ]
  'DarkRed'                 '#8B0000'    [ 139   0   0 ]
  'light green'             '#90EE90'    [ 144 238 144 ]
  'LightGreen'              '#90EE90'    [ 144 238 144 ]

  };

