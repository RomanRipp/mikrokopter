function [Msg,c,s] = char2img(c,varargin);


% CHAR2IMG Converts Single Character to Image
%
% [Msg,IMG,SIZE] = CHAR2IMG( C , FontProperty , FontPropertyValue , ... );
%
%  C    = N-Element Chararcter- or CellStringArray, with single Characters
%
%  IMG  = [ N by 2 ] CellArray: { Character  ImageData }, ImageValues = 0 | 1
%  SIZE = [ N by 2 ] Matrice: [ Width Height ] of CharacterImage
%
% see also: STR2IMG
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
%  fig = figure('position',[200 200 400 300],'color',[1 1 1]);
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

Msg = '';

s = zeros(0,2);

Nin = nargin;

if Nin == 0

  Msg = 'Input S is undefined.';
  return

end


is_win = strcmp( upper(computer) , 'PCWIN' );

%------------------------------------------
% Check Characters

if isempty(c)
   return
end

if ischar(c)

   c     = c(:);
  
     ii  = find( double(c) == 32 );

   c(ii) = 'H';

   c     = cellstr(c(:));

   c(ii) = { char(32) };

end

ok = iscellstr(c);

if ok

   c = c(:);

   n = size(c,1);

   try

     d  = c;
     d( find(strcmp(d,char(32))) ) = { 'H' };
 
     ok = isequal( size(cat(1,d{:})) , [ n  1 ] );

   catch

     ok = 0;

   end

end

if ~ok

   Msg = 'Input C must be a Char- or CellStringArray, containing single Characters.';

   return

end

%------------------------------------------
% Check additional Properties

VarIn = varargin;

if ~isempty(VarIn)

  VarIn = VarIn(:);

    Nin = size(VarIn,1);

     ok = ( mod(Nin,2) == 0 );

  if ok

       Nin = Nin / 2;

     VarIn = reshape(VarIn,2,Nin);

        ok = iscellstr(VarIn(1,:));

     if ok

        chk = 'font';
        m   = size(chk,2);

        p = double( lower( char( VarIn(1,:)' ) ) );
        ok = ( size(p,2) > m );
        if ok
           ok = ( ones(Nin,1)*double(chk) == p(:,1:m) );
           ok = ( sum(ok(:)) == Nin*m );
        end

     end

  end

  if ~ok

     Msg = 'Additional Inputs must be Font-Property-Value-Pairs.';
  
     return

  end

end


%*********************************************************

cmap = [ 0  0  0 
         1  1  1  ];

fig = figure('units'      , 'pixels'          , ...
             'position'   , [ 100 100 10 10 ] , ...
             'menubar'    , 'none'            , ...
             'createfcn'  , ''                , ...
             'color'      , [ 1  1  1 ]       , ...
             'colormap'   , cmap              , ...
             'nextplot'   , 'add'             , ...
            'windowstyle' , 'modal'                 );

axe = axes( 'parent'     , fig              , ...
            'units'      , 'pixels'         , ...
            'position'   , [ 1  1  10  10 ] , ...
            'color'      , 'none'           , ...
            'visible'    , 'off'            , ...
            'nextplot'   , 'add'                   );


%--------------------------------------------------------

try

  ht = text( 'parent'    , axe         , ...
             'units'     , 'pixels'    , ...
             'position'  , [ 1 1 0 ]   , ...
             'color'     , [ 0  0  0 ] , ...
             'string'    , ''          , ...
             'clipping'  , 'off'       , ...
             'visible'   , 'on'        , ...
             'horizontalalignment' , 'left'   , ...
             'verticalalignment'   , 'middle' , ...
              VarIn{:}                              );        

catch

 Msg = lasterr;

end

if ~isempty(Msg)

   Msg = [ 'Invalid FontPropertyInputs.' char(10) Msg ];

   delete(axe)
   delete(fig)

   return

end


%--------------------------------------------------------

c = c(:,[1 1]);

for ii = 1 : n

    ext = zeros(2,4);

    for jj = [ 1  2 ];

        set( ht , 'string' , char( double(c{ii,1})*ones(jj,jj) ) );

        ext(jj,:) = get(ht,'extent');

    end

    set( ht , 'string' , c{ii,1} );

    wh = ext(2,[3 4]) - ext(1,[3 4]);

    pos = cat( 2 , 1 , floor(wh(2)/2)+1 , 0 );
    
    pos([1 2]) = pos([1 2]) + 1*is_win;
    
    set(ht,'position',pos); 

    set( fig , 'position' , [ 100 100 wh ] );


    c{ii,2} = getframe(fig);   

    sc = size(c{ii,2}.cdata);
    
    i1 = ( 1 : wh(2) ) + sc(1) - wh(2);
    i2 = ( 1 : wh(1) );
    
    c{ii,2}.cdata = c{ii,2}.cdata(i1,i2,:);
    
    if isempty(c{ii,2}.colormap)

       c{ii,2} = 1 - ( double(c{ii,2}.cdata(:,:,1)) == 0 );

    else

       c{ii,2} = ind2rgb( c{ii,2}.cdata , c{ii,2}.colormap );

       c{ii,2} = 1 - ( c{ii,2}(:,:,1) == 0 );
        
    end

    s(ii,:) = [ size(c{ii,2},2)  size(c{ii,2},1) ];
    
end

%--------------------------------------------------------

delete(ht)
delete(axe)
delete(fig)
