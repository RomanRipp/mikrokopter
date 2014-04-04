function [ok,str,txt] = comment(varargin)

% COMMENT  Opens DialogWindow to enter or modify a Commentary
%
%  [Ok,CommentArray,CommentString] = COMMENT( FigureName , Titel , Text , ...
%                                              Comment   , [ NRow NCol ] );
%
%  OK             1  if  Window finished by Ok
%                 0   otherwise (Cancel or Closed)
%
%  CommentArray   CharArray or CellStringArray of Comment
%
%  CommentString  String of CommentArray, 
%                 Rows are  separated by NewLine-Character ( char(10) )
%


nl = char(10);

VarIn    = cell(5,1);
VarIn(:) = { '' };

n = min( size(VarIn,1) , prod(size(varargin)) );

VarIn(1:n) = varargin(1:n);

name = VarIn{1};
if isempty(name)
  name = 'Please Enter a Comment';
end

si = VarIn{5};

VarIn = VarIn([2 3 4]);

VarIn = cat(1,VarIn(:),{'Cancel'});


if ~isempty(si)

  ok = ( isnumeric( si )  &  ( prod(size(si)) == 2 ) );
  if ok
     ok = all( isfinite(si)  &  ( si >= 1 )  & ( mod(si,1) == 0 ) );
  end

  if ~ok
     error([ 'Size [ NRow NCol ] must be a 2-Element numeric Vector ' nl ...
             'with Elements larger ZERO.' ]);
  end 

end


FontUnits  = 'points';

[FontSize,ssi,ppi,is_win,fn]   = bestfont;


FontIni = { 'helvetica'   'bold'    'Titel'
            'helvetica'   'normal'  'Text'
             fn           'normal'  'Edit'
            'helvetica'   'normal'  'Control'  };

val = 0;

UIColor  = [ 0.9  0.95  1.0 ];

UIBorderWidth = 3;



     % StringFormat for Handle
     clear eps
     Hform = sprintf( '%%.%.0fg' , ceil(abs(log(eps)/log(10)))+1 );



fig = figure('color'       , [ 1  1  1 ], ...
             'numbertitle' , 'off' , ...
             'menubar'     , 'none', ...
             'toolbar'     , 'none', ...
             'name'        ,  name , ...
             'resize'      , 'off' , ... 
             'visible'     , 'off' , ...
             'createfcn'   , ''    , ...
             'handlevisibility' , 'callback' , ...
             'WindowStyle'      , 'modal'           );

%--------------------------------------------------

ResumeStatus = 'inactive';

ResumeCB     = [ 'set('  sprintf(Hform,fig) ',''waitstatus'','''  ResumeStatus  ''');' ];


%***************************************************
% Right Positioning of UIControls

  % DummyAxes

   axe = axes( 'parent'   , fig          , ...
               'units'    , 'normalized' , ...
               'position' , [ 0 0 1 1 ]  , ...
               'visible'  , 'off'              );
  % DummyText

  ht = text('parent'     , axe        , ...
            'units'      ,'pixels'    , ...
            'position'   , [ 1 1 0 ]  , ...
            'string'     , ''         , ...
            'FontUnits'  , FontUnits  , ...
            'FontSize'   , FontSize   , ...
            'FontAngle'  , get(0,'DefaultUIControlFontAngle')  , ...
            'visible'    , 'off'      , ...
            'interpreter' , 'none'                              );

  NV = size(VarIn,1);

  ext = zeros(NV,2);

  Msg = '';  
  for ii = 1 : NV

    if ~( ( iscellstr(VarIn{ii})  | ischar(VarIn{ii}) )  & ...
          ( prod(size(VarIn{ii})) == size(VarIn{ii},1)*size(VarIn{ii},2) )   )

       Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
               'Input '   FontIni{ii,3}  ...
               ' must be a CharArray or CellStringArray.' ];

    else
 
      str = rmblank( char(VarIn{ii}) , 2 );

      if strcmp( FontIni{ii,3} , 'Edit' )  

         if isempty(si)

            si = [ 6  60 ]; 

            if ~isempty(str)
              inl = cat( 2 , ( double(str) == 10 ) , ones( size(str,1) , 1 ) );
              inl = inl';
              inl = cat( 1 , 1 , inl(:) );
              si  = max( si , [ sum(inl,1)-1  max(diff(find(inl))-1) ] );
            end 
  
         end

         str = char( 'H' * ones(si(1),si(2)) );
 
      end

      set( ht , 'FontName'   , FontIni{ii,1}   , ...
                'FontWeight' , FontIni{ii,2}   , ...
                'string'     , str         );

      e         = get( ht , 'extent' );

      ext(ii,:) = e([3 4]);
 
    end

  end

  delete(ht)
  delete(axe)

  if ~isempty(Msg)
    delete(fig); 
    error(Msg)
    return
  end


 is_titel = ( ext(1,2) > 0 );
 is_text  = ( ext(2,2) > 0 );

 ext(:,1) = ext(:,1) + 4/3 * ext(4,2);
 
 ext      = ext + 2*UIBorderWidth * ( ext > 0 );

       ww = ext(:,1);
 
 UIHeight = 1.2 * ext(4,2);
 
 UIWidth          = ww;
 UIWidth(4)       = 3*ww(4) + 2*UIHeight/2 ;
 UIWidth([1 2 4]) = UIWidth([1 2 4]) + UIHeight;

 ext(3,1) = max(UIWidth);

 xoff = 5;

FigHeight = UIHeight + sum(ext(:,2));
FigWidth  = max(UIWidth);

      ld = ( FigWidth - 3*ww(4) ) / 4;  % Horizontal Distance between lower Controls

FigWidth = FigWidth + 2*xoff;


FigPos = [ NaN  NaN   FigWidth  FigHeight ];  

FigPos([3 4]) = FigPos([3 4]) + ( 5/6*ssi([3 4]) - FigPos([3 4]) ) .* ...
                                ( 5/6*ssi([3 4]) < FigPos([3 4]) );

FigPos([1 2]) = floor( ( ssi([3 4]) - FigPos([3 4]) ) / 2 );


      ww = ww * FigPos(3) / FigWidth ;
      ld = ld * FigPos(3) / FigWidth;

UIHeight = UIHeight * FigPos(4) / FigHeight;

% FontSize = FontSize * min(FigPos([3 4])./[FigWidth FigHeight]);

ext(3,2) = ext(3,2) + FigPos(4) - FigHeight;


set( fig , 'units'       , 'pixels' , ...
           'position'    , FigPos         );
             
%******************************************************************

if ( ext(1,2) > 0 )

  HTitel = uicontrol('Parent',fig, ...
       	'Units','pixels', ...
	'BackgroundColor',[ 1  1  1 ], ...
	'Position'   , ...
  [ (FigPos(3)-ww(1))/2  1.5*UIHeight+ext(2,2)+ext(3,2)+1  ww(1)  ext(1,2) ], ...
	'String'     , VarIn{1}, ...
        'FontUnits'  , FontUnits  , ...
        'FontSize'   , FontSize   , ...
        'FontName'   , FontIni{1,1} , ...
        'FontWeight' , FontIni{1,2} , ...
        'horizontalalignment' , 'center' ,  ...
	'Style','text', ...
        'min'  , 0 , ...
        'max'  , 1  );
end

if ( ext(2,2) > 0 )

  HText = uicontrol('Parent',fig, ...
       	'Units','pixels', ...
	'BackgroundColor',[ 1  1  1 ], ...
	'Position'   , ...
  [ (FigPos(3)-ww(2))/2  1.5*UIHeight+ext(3,2)+1  ww(2)  ext(2,2) ], ...
	'String'     , VarIn{2}     , ...
        'FontUnits'  , FontUnits    , ...
        'FontSize'   , FontSize     , ...
        'FontName'   , FontIni{2,1} , ...
        'FontWeight' , FontIni{2,2} , ...
        'horizontalalignment' , 'center' ,  ...
	'Style','text', ...
        'min'  , 0 , ...
        'max'  , 1  );
end

HEdit = uicontrol('Parent',fig, ...
	'Units','pixels', ...
	'BackgroundColor',UIColor, ...
	'Position'   ,[ 1+xoff  1.5*UIHeight+1 ext(3,1)  ext(3,2) ], ...
	'String'     , VarIn{3}, ...
        'FontUnits'  , FontUnits  , ...
        'FontSize'   , FontSize   , ...
        'FontName'   , FontIni{3,1} , ...
        'FontWeight' , FontIni{3,2} , ...
        'horizontalalignment' , 'left' ,  ...
	'Style','edit', ...
        'min'  , 0    , ...
        'max'  , 2    , ...
        'userdata' , VarIn{3}  );


%------------------------------------------------------------------

HOk = uicontrol('Parent',fig, ...
	'Units','pixels', ...
	'BackgroundColor',UIColor, ...
	'Position',[ ld+xoff 1/4*UIHeight  ww(4)  UIHeight], ...
	'String'     , 'Ok'         , ...
        'FontUnits'  , FontUnits    , ...
        'FontSize'   , FontSize     , ...
        'FontName'   , FontIni{4,1} , ...
        'FontWeight' , FontIni{4,2} , ...
        'horizontalalignment' , 'center'    , ...
	'Style'               , 'pushbutton', ...
        'userdata'            , 0                  );

CB = [ 'set('  sprintf(Hform,HOk) ',''userdata'',1);'  ResumeCB ];

set( HOk , 'callback' , CB );


%------------------------------------------------------------------

CB = [ 'set('  sprintf(Hform,HEdit) ',''string'',' ...
       'get('  sprintf(Hform,HEdit) ',''userdata''));'    ];

HReset = uicontrol('Parent',fig, ...
	'Units','pixels', ...
	'BackgroundColor',UIColor, ...
	'Position',[ 2*ld+ww(4)+xoff  1/4*UIHeight  ww(4)  UIHeight], ...
	'String'     , 'Reset', ...
        'FontUnits'  , FontUnits    , ...
        'FontSize'   , FontSize     , ...
        'FontName'   , FontIni{4,1} , ...
        'FontWeight' , FontIni{4,2} , ...
        'horizontalalignment' , 'center'     , ...
	'Style'               , 'pushbutton' , ...
        'callback'            , CB                   );

%------------------------------------------------------------------

HCancel = uicontrol('Parent',fig, ...
	'Units','pixels', ...
	'BackgroundColor',UIColor, ...
	'Position',[ FigWidth-ld-ww(4)-xoff  1/4*UIHeight  ww(4)  UIHeight], ...
	'String'     , 'Cancel'     , ...
        'FontUnits'  , FontUnits    , ...
        'FontSize'   , FontSize     , ...
        'FontName'   , FontIni{4,1} , ...
        'FontWeight' , FontIni{4,2} , ...
        'horizontalalignment' , 'center'     , ...
	'Style'               , 'pushbutton' , ...
        'callback'            , ResumeCB           );


%--------------------------------------------------

set( fig , 'visible'         , 'on'      , ...
           'waitstatus'      , 'waiting' , ...
           'CloseRequestFcn' , ResumeCB           );

waitfor( fig , 'waitstatus' , ResumeStatus );

%--------------------------------------------------

ok  = get( HOk   , 'userdata' );
str = get( HEdit , 'string'   );

if ischar(str)  &  ( size(str,1) == 1 )

   txt = str;

else
  
   txt        = cellstr(str);
   txt( : ,2) = { nl };
   txt(end,2) = { '' };
 
   txt = txt';

   txt = cat(2,txt{:});

end

delete(fig);



%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = rmblank(str,dim,cc)

% RMBLANK  Remove Blanks, NewLines at Begin and End of CharacterArrays
%
% String = RMBLANK( CharArray )
%
% CharArray  2-dimensional CharacterArray
%
% further Options:
%
% String = RMBLANK( CharArray , DIM , CHAR )
%
%  
%  DIM  specifies Dimension to work, 
%       default: 2
%
%  CHAR specifies BlankCharacters to remove
%       default:  [ 32  13  10  9 ];  % [ Space CR LF TAB ]
%

  
msg = '';
 nl = char(10);

Nin = nargin;

if Nin < 1
  error('Not enough Input Arguments.')
else
  str0 = str;
  if ischar(str)
    str = double(str);
  end
  ok = isnumeric(str);
  if ok
    ok = all( ( mod(str(:),1) == 0 )  & ...
              ( str(:) >= 0 ) & isfinite(str(:))  );
  end
  if ~ok
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CharArray must be a String or ASCII-Codes.'];
  end
  if size(str,1)*size(str,2) ~= prod(size(str))
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CharArray must be 2-dimensional.'];
  end     
end

if Nin < 2
  dim = 2;
else
  if ~isnumeric(dim)
    msg = [ msg nl(1:(end*(~isempty(msg)))) ...
            'Input DIM must be numeric.' ];
  elseif ~isempty(dim)
    dim = dim(:)';
    if ~all( ( dim == 1 ) |  ( dim == 2 ) )
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
             'Values for Input DIM must Integers larger ZERO.' ];
    end
  end 
end

if Nin < 3
  cc = [ 32  13  10  9 ];  % [ Space CR LF TAB ]
else
  if ischar(cc)
    cc = double(cc);
  end
  ok = isnumeric(cc);
  if ok & ~isempty(cc)
    cc = cc(:)';
    ok = all( ( mod(cc,1) == 0 )  & ...
              ( cc >= 0 ) & isfinite(cc)  );
  end
  if ~ok
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CHAR must be a String or ASCII-Codes.'];
  end
end

if ~isempty(msg)
  error(msg)
end


if isempty(str)
 str = '';
 return
end

if isempty(dim) | isempty(cc)
  str = str0;
  return
end

     jj  = find(str == 0 );
 str(jj) = cc(1);

  blank  = 0*str;

  for ii = cc
    blank = ( blank | ( str == ii ) );
  end

  si = size(str);

  for d = dim

    bad = ( sum(blank,3-d) == si(3-d) );
    jj  = find( bad );
    if ~isempty(jj) 

         p  = [ 3-d  d ];
        str = permute(str,p);

         jj = jj(:)';
         nb = size(jj,2);
        jj1 = find( jj ==   ( 1 : nb ) );       % Blank at Begin
        jj2 = find( jj == ( ( 1 : nb ) + ...    % Blank at End
                            ( si(d) - nb ) ) );
        str(:,jj([jj1 jj2])) = [];

        str = permute(str,p);

    end
    
  end

  str = char(str);


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  str = strhcat(str,del,n,nl)

% STRHCAT  Concatenates Strings into ONE
%
% STRHCAT( StringArray , Delimiter )
%   Forms one long String from the Strings in the
%   StringArray, delimited with the delimiter.
%   The EndDelimiter will be removed.
%
% STRHCAT( StringArray , Delimiter , N , NewLine )
%   Build a  NewLine after each N-th String.
%   default: N = 10;  NewLine = char(10);
%
% Example:  
%         >> strhcat({'apples' 'pies' 'prunes'},', ')
%    
%         ans =
%
%         apples, pies, prunes
%
%         >> strhcat({'apples';'pies';'prunes'},', ',2)
%    
%         ans =
%
%         apples, pies
%         prunes
%



Nin = nargin;

if Nin < 4
 nl = char(10);
end
if Nin < 3
 n = 10;
end
if Nin < 2
 del = char((32*ones(1,3)));
end


if isempty(str)
 str = '';
 return
end


if ischar(str)
  str = cellstr(str);
end

str = str(:);

str( : , 2 ) = { del };

str(n:n:size(str,1),2) = {nl};


str(    size(str,1),2) = {''};


str = str';

str = cat(2,str{:});


