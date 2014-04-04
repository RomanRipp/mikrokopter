function val = list_sel(varargin)

% LIST_SEL  Opens DialogWindow to select from a List and returns Value
%
%  Value = LIST_SEL( List , FigureName , Titel , Text );
%
%  The List of Strings, which will shown in a PopUpMenu 
%   to select one of them. This Value will returned.
% 


if nargin < 1
 error('No enough Input Arguments.');
end

nl = char(10);

VarIn    = cell(4,1);
VarIn(:) = { '' };

n = min( size(VarIn,1) , prod(size(varargin)) );

VarIn(1:n) = varargin(1:n);

name = VarIn{2};
if isempty(name)
  name = 'Please Select';
end

VarIn = VarIn([3 4 1]);

VarIn = cat(1,VarIn(:),{'Cancel'});


FontUnits  = 'points';
FontSize   = 12;

FontIni = { 'helvetica'   'bold'    'Titel'
            'helvetica'   'normal'  'Text'
            'courier'     'normal'  'List'
            'helvetica'   'normal'  'Control'  };

val = 0;

UIColor  = [ 0.8  0.9  1.0 ];

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

    if ~( iscellstr(VarIn{ii})  |  ischar(VarIn{ii}) )

       Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
               'Input '   FontIni{ii,3}  ...
               ' must be a CharArray or CellStringArray.' ];

    else
 
      set( ht , 'FontName'   , FontIni{ii,1}   , ...
                'FontWeight' , FontIni{ii,2}           );
            
      str = rmblank( char(VarIn{ii}) , 2 );

      if strcmp( FontIni{ii,3} , 'List' )

        for jj = 1 : size(str,1)
          set( ht , 'string' , str(jj,:) );
          e         = get( ht , 'extent' );
          e(3)      = e(3) + 2*2/3*e(4);
          ext(ii,:) = e([3 4]) + ( ext(ii,:) - e([3 4]) ) .* ( ext(ii,:) > e([3 4]) );
        end
  
      else

        set( ht , 'string' , str       );
        e         = get( ht , 'extent' );
        e(3)      = e(3) + 2*2/3*e(4);
        ext(ii,:) = e([3 4]);
 
      end

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


 ext      = ext + 2*UIBorderWidth * ( ext > 0 );

 ext(3,1) = ext(3,1) + ext(3,2);  % PopUp !!!

       ww = ext(:,1);
 
 UIHeight = max(ext([ 3  4 ],2));
 
 UIWidth  = max( cat( 1 , ww([1 2 3]) , (2*ww(4)+UIHeight/2) ) );

FigHeight = 1.5 * UIHeight + sum(ext(:,2));
FigWidth  = UIWidth + UIHeight;

       ld = ( FigWidth - 2*ww(4) ) / 3;  % Horizontal Distance between lower Controls


uni = get(0,'units');
      set(0,'units','pixels');
ssi = get(0,'screensize');
      set(0,'units',uni);


FigPos = [ NaN  NaN   FigWidth  FigHeight ];  

FigPos([3 4]) = FigPos([3 4]) + ( 2/3*ssi([3 4]) - FigPos([3 4]) ) .* ...
                                ( 2/3*ssi([3 4]) < FigPos([3 4]) );

FigPos([1 2]) = floor( ( ssi([3 4]) - FigPos([3 4]) ) / 2 );


      ww = ww * FigPos(3) / FigWidth ;
      ld = ld * FigPos(3) / FigWidth;

UIHeight = UIHeight * FigPos(4) / FigHeight;
FontSize = FontSize * min(FigPos([3 4])./[FigWidth FigHeight]);


set( fig , 'units'       , 'pixels' , ...
           'position'    , FigPos         );
             
%******************************************************************

if ( ext(1,2) > 0 )

  HTitel = uicontrol('Parent',fig, ...
       	'Units','pixels', ...
	'BackgroundColor',[ 1  1  1 ], ...
	'Position'   ,[ (FigPos(3)-ww(1))/2  3*UIHeight+ext(2,2)  ww(1)  ext(1,2) ], ...
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
	'Position'   ,[ (FigPos(3)-ww(2))/2  3*UIHeight  ww(2)  ext(2,2) ], ...
	'String'     , VarIn{2}, ...
        'FontUnits'  , FontUnits  , ...
        'FontSize'   , FontSize   , ...
        'FontName'   , FontIni{2,1} , ...
        'FontWeight' , FontIni{2,2} , ...
        'horizontalalignment' , 'center' ,  ...
	'Style','text', ...
        'min'  , 0 , ...
        'max'  , 1  );
end

HList = uicontrol('Parent',fig, ...
	'Units','pixels', ...
	'BackgroundColor',UIColor, ...
	'Position'   ,[ (FigPos(3)-ww(3))/2  2*UIHeight  ww(3)  UIHeight ], ...
	'String'     , VarIn{3}, ...
        'FontUnits'  , FontUnits  , ...
        'FontSize'   , FontSize   , ...
        'FontName'   , FontIni{3,1} , ...
        'FontWeight' , FontIni{3,2} , ...
        'horizontalalignment' , 'left' ,  ...
	'Style','popupmenu', ...
        'min'  , 0 , ...
        'max'  , 1 , ...
        'userdata' , 0  );


CB = [ 'set('  sprintf(Hform,HList) ',''userdata'',1);'  ResumeCB  ];


HOk = uicontrol('Parent',fig, ...
	'Units','pixels', ...
	'BackgroundColor',UIColor, ...
	'Position',[ ld 0.5*UIHeight  ww(4)  UIHeight], ...
	'String','Ok', ...
        'FontUnits'  , FontUnits  , ...
        'FontSize'   , FontSize   , ...
        'FontName'   , FontIni{4,1} , ...
        'FontWeight' , FontIni{4,2} , ...
        'horizontalalignment' , 'center' ,  ...
	'Style','pushbutton', ...
        'callback' , CB );

HCancel = uicontrol('Parent',fig, ...
	'Units','pixels', ...
	'BackgroundColor',UIColor, ...
	'Position',[ FigWidth-ld-ww(4)  0.5*UIHeight  ww(4)  UIHeight], ...
	'String','Cancel', ...
        'FontUnits'  , FontUnits  , ...
        'FontSize'   , FontSize   , ...
        'FontName'   , FontIni{4,1} , ...
        'FontWeight' , FontIni{4,2} , ...
        'horizontalalignment' , 'center' ,  ...
	'Style'      , 'pushbutton', ...
        'callback'   , ResumeCB            );


%--------------------------------------------------

set( fig , 'visible'         , 'on'      , ...
           'waitstatus'      , 'waiting' , ...
           'CloseRequestFcn' , ResumeCB           );

waitfor( fig , 'waitstatus' , ResumeStatus );


%--------------------------------------------------

if get(HList,'userdata')
  val=get(HList,'value');
end

% set(fig,'visible','off');

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
