function [Msg,img] = str2img( str , c , s );

% STR2IMG  Convert String to Image
%
% [ Msg , IMG ] = STR2IMG( STR , C , SIZE );
%
%  STR  = String to convert
%  C    = [ N by 2 ] CellArray: { Character  ImageData }, ImageValues = 0 | 1
%  SIZE = [ N by 2 ] Matrice: [ Width Height ] of CharacterImage
% 
%  C and SIZE are outputs from CHAR2IMG, all Characters of STR should be 
%   contained in 1. Column of C.
%
% see also: CHAR2IMG
%
%-----------------------------------------------------------------------
% Example:
%
%  FontIn = { 'FontName'   , 'helvetica' , ...
%             'FontSize'   , 12          , ...
%             'FontWeight' , 'bold'            };
%
%  [Msg,c,s] = char2img( ' 0123456789.:%' , FontIn{:} );
%
%  cl  = clock;  % [ YY MM DD hh mm ss ]
%
%  pc  = floor( 100 * ( cl(4) + cl(5)/60 + cl(6)/3600 ) / 24 );  % PercentDay
%
%  cl  = [ floor(cl(4:6)) pc ];  % [ hh mm ss PercenDay ]
%
%  str = sprintf( '%2.0f:%2.2d:%2.2d   %3.0f%%' , cl );
%
%  [Msg,img] = str2img( str , c , s );
%  
%  fig = figure('units','pixels','position',[200 200 400 300],'color',[1 1 1]);
%
%  offs = 10;
%    ww = size(img,2) + 2*offs;
%    hh = size(img,1) + 2*offs;
%   
%  uicontrol( 'parent'   , fig                 , ...
%             'units'    , 'pixels'            , ...
%             'position' , [ 10  10  ww  hh ]  , ...
%             'style'    , 'pushbutton'        , ...
%             'string'   , str                 , ...
%              FontIn{:}                       , ...
%             'horizontalalignment' , 'center' , ...
%             'backgroundcolor'     , [ 0.5  1.0  1.0 ] , ...
%             'foregroundcolor'     , [ 0.0  0.0  0.0 ] , ...
%             'tooltipstring'       , 'Text'               ); 
%
%  cdata        = zeros( hh , ww , 3 );
%  cdata(:,:,1) = 1;
%
%  i1 = ( 1 : size(img,1) ) + offs;
%  i2 = ( 1 : size(img,2) ) + offs;
% 
%  cdata(i1,i2,:) = cdata(i1,i2,:) .* img(:,:,[1 1 1]);
%
%  uicontrol( 'parent'   , fig                 , ...
%             'units'    , 'pixels'            , ...
%             'position' , [ 10  10+hh+10  ww  hh ]  , ...
%             'style'    , 'pushbutton'        , ...
%             'cdata'    , cdata               , ...
%             'tooltipstring' , 'Image'                  );
%

Nin = nargin;

Msg = '';
img = zeros(0,0);

if Nin < 2
  Msg = 'Not enough Input Arguments.';
  return
end


%------------------------------------------
% Check String

if iscellstr(str)
    str = char(str);
end

if isempty(str)
   return
end

m = size(str,2);

if ~( ischar(str)  &  ( prod(size(str)) == m )  )

  Msg = '1. Input must be a String.';

  return

end

%------------------------------------------
% Check Characters

if isempty(c)
   return
end

n = size(c,1);

ok = ( iscell(c)  &  ( size(c,2) == 2 )  &  ...
       ( prod(size(c)) == 2*n )    );

if ok
   ok = iscellstr(c(:,1));
end

if ~ok

   Msg = 'Input C must be 2-Column CellArry: { Character ImageData }.';

   return

end


%------------------------------------------
% Check Size

if Nin < 3   

   s = [];

else

  ok = ( isnumeric(s)  &  isequal( size(s) , [n 2] ) );
  if ok
     ok = ( all(isfinite(s(:)))  &  all( s(:) >= 0 ) );  
  end

  if ~ok

     s = [];

  end

end

if isempty(s)

   s = zeros(n,2);

   for ii = 1 : n

     s(ii,:) = [ size(c{ii,2},2) size(c{ii,2},1) ];

   end

end

%*************************************************


ind = zeros(1,m);

for ii = 1 : m

  ind(ii) = sum( cumsum( strcmp(str(ii),c(:,1)) ) == 0 ) + 1;

end

ind( find( ind > n ) ) = [];

if all(  ( s(:,2) - s(1,2) ) == 0  )

   img = cat( 2 , c{ind,2} );

   return

end


h = max(s(:,2));

img = ones( h , sum(s(ind,1)) );

m   = size(ind,2);

i0  = cumsum([ 0 ; s(ind,1) ]) + 1;

for ii = 1 : m

    img( h-s(ind(ii),2)+1 : h , ( i0(ii) : i0(ii+1)-1 ) ) = c{ind(ii),2};

end
