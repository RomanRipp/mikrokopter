function [ok,val,str,txt] = propdlg(name,ttl,txt,frm,str,siz)

% PROPDLG  Opens DialogWindow to enter or modify a Property
%
%  [OK,Value,String,Text] = PROPDLG( FigureName , Titel  , Comment , ...
%                                    Format , String/Value , [ NRow NCol ] );
%
%
%  OK             1  if  Window finished by Ok
%                 0   otherwise (Cancel or Closed)
%
%  Text          String of Editor
%
% [String,Value] = VAL2STR(Text,Format)
%
% See also: VAL2STR (required), STR2VAL
%

Nin = nargin;

if Nin < 1, name = ''; end
if Nin < 2, ttl  = ''; end
if Nin < 3, txt  = ''; end
if Nin < 4, frm  = ''; end
if Nin < 5, str  = ''; end
if Nin < 6, siz  = []; end

%****************************************************************
% Check Inputs

msg = cell(0,1);

if isempty(name)
   name = 'Property-Editor';
elseif ~chkstr(name)
   msg = cat(1,msg,{'Title must be a String.'});
end


[ok,ttl] = chkcstr(ttl,0);
if ~( ok | isempty(ttl) ) 
    msg = cat(1,msg,{'Title must be a CharacterArray or CellArray of Strings.'});
elseif isempty(ttl)
    ttl = '';
end

[ok,txt] = chkcstr(txt,0);
if ~( ok | isempty(txt) ) 
    msg = cat(1,msg,{'Text must be a CharacterArray or CellArray of Strings.'});
elseif isempty(txt)
    txt = '';
end

if isempty(frm)
   frm = 'char';
elseif ~chkstr(frm) 
   msg = cat(1,msg,{'Format must be a String.'});
end

if chkstr(frm)
   try
      [val,str,frm,m] = str2val(str,frm);
   catch
       m = lasterr;
   end
   if ~isempty(m)
       msg = cat( 1 , msg , {sprintf('Error using STR2VAL.\n%s',m)} );
   end
else
   val = str;
end

if ~isempty(siz)
    ok = ( isnumeric( siz )  &  ( prod(size(siz)) == 2 ) );
    if ok
       ok = all( isfinite(siz)  &  ( siz >= 1 )  & ( mod(siz,1) == 0 ) );
    end

    if ~ok
        msg = cat(1,msg,{'Size must be a 2-Element Integer larger ZERO.'})
    end 
elseif ischar(str)
    siz = size(str);
end

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Invalid Inputs.\n%s',msg));
end

%****************************************************************

FontUnits  = 'points';

[FontSize,ssi,ppi,is_win,fn]   = bestfont;

FontIni = { 'helvetica'   'bold'    'Titel'
            'helvetica'   'normal'  'Text'
             fn           'normal'  'Edit'
            'helvetica'   'normal'  'Control'  };

Label = { ttl  txt  str  'Cancel' };
  
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

  nl = prod(size(Label));

  ext = zeros(nl,2);

  siz = siz + [ 1  32 ] .* ( siz == 0 );

  for ii = 1 : nl

      if strcmp( FontIni{ii,3} , 'Edit' )  
         lab = char( 'H' * ones(siz) );
      else
         lab = rmblank( char(Label{ii}) , 2 );
      end


      set( ht , 'FontName'   , FontIni{ii,1}   , ...
                'FontWeight' , FontIni{ii,2}   , ...
                'string'     , lab         );

      e         = get( ht , 'extent' );

      ext(ii,:) = e([3 4]);

  end

  delete(ht)
  delete(axe)


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
	'String'     , Label{1}, ...
        'FontUnits'  , FontUnits  , ...
        'FontSize'   , FontSize   , ...
        'FontName'   , FontIni{1,1} , ...
        'FontWeight' , FontIni{1,2} , ...
        'horizontalalignment' , 'center' ,  ...
        'tag' , 'Title' , ...
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
	'String'     , Label{2}     , ...
        'FontUnits'  , FontUnits    , ...
        'FontSize'   , FontSize     , ...
        'FontName'   , FontIni{2,1} , ...
        'FontWeight' , FontIni{2,2} , ...
        'horizontalalignment' , 'center' ,  ...
        'tag' , 'Text' , ...
	'Style','text', ...
        'min'  , 0 , ...
        'max'  , 1  );
end

CB = 'set(gcbo,''string'',val2str(get(gcbo,''string''),get(gcbo,''UserData'')))';

ht = { 'center'  'left' };

mx = 1 + ( siz(1) > 1 ) ;

ht = ht{mx};

HEdit = uicontrol('Parent',fig, ...
	'Units','pixels', ...
	'BackgroundColor',UIColor, ...
	'Position'   ,[ 1+xoff  1.5*UIHeight+1 ext(3,1)  ext(3,2) ], ...
	'String'     , Label{3}, ...
        'FontUnits'  , FontUnits  , ...
        'FontSize'   , FontSize   , ...
        'FontName'   , FontIni{3,1} , ...
        'FontWeight' , FontIni{3,2} , ...
        'horizontalalignment' , ht , ...
        'tag'  , 'Edit' , ...
	'Style','edit', ...
        'min'  , 0    , ...
        'max'  , mx , ...
        'CallBack' , CB , ...
        'userdata' , frm  );


%------------------------------------------------------------------

ResumeStatus = 'inactive';

ResumeCB     = [ 'set(get(gcbo,''parent''),''waitstatus'','''  ResumeStatus  ''');' ];

CB = [ 'set(gcbo,''userdata'',1);'  ResumeCB ];

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
        'tag'                 , 'Ok' , ...
	'Style'               , 'pushbutton', ...
        'CallBack'            , CB  , ...
        'userdata'            , 0                  );

%------------------------------------------------------------------

CB = 'set(findobj(get(gcbo,''parent''),''tag'',''Edit''),''string'',get(gcbo,''userdata''))';

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
        'tag'        , 'Reset' , ...
	'Style'      , 'pushbutton' , ...
        'UserData'   , str , ... 
        'callback'   , CB                   );

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
        'tag'        , 'Cancel' , ...
	'Style'      , 'pushbutton' , ...
        'callback'   , ResumeCB           );


%--------------------------------------------------

set( fig , 'visible'         , 'on'      , ...
           'waitstatus'      , 'waiting' , ...
           'CloseRequestFcn' , ResumeCB           );

waitfor( fig , 'waitstatus' , ResumeStatus );

%--------------------------------------------------

ok  = get( HOk   , 'userdata' );
txt = get( HEdit , 'string'   );

[str,val] = val2str(txt,frm);


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


