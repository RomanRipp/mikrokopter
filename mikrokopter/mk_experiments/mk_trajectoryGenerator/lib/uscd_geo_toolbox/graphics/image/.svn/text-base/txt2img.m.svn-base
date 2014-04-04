function  [Msg,c] = txt2img(txt,varargin);


% TXT2IMG  Creates ImageCData or MultiLines Text
%
%  [ Msg, CDATA ] = TXT2IMG( TEXT )
%
%  TXT2IMG opens a Figure, creates a TextUIControl with the
%   string TEXT. TEXT could be a CharArray or CellStringArray.
%  The Number of Rows of TEXT is the Number of Rows of the Text
%   in the UIControl. CDATA contains the RGB-ColorValues of the 
%   TextUIControl, for using in the CData-Property for PushButtons.
%   CData = [ PixelHigh by PixelWidth by 3 ] double
%      3. Dimension == [ R G B ],   0 <= CData <= 1 
%
%  Msg contains ErrorMessages, if the Inputs are invalid.
%  Msg is empty, if no error detected.
%
%  The Figure, which has been opend, will be closed automaticly.
%
%  Note:  TXT2IMG does temporaly changes in the Root 
%          DefaultUIControlProperties. The originaly Properties
%          will set to default automaticly. If TXT2IMG terminates
%          before, you have to take care if the Root Properties are correct.
%
%
%  Additional Inputs:  TXT2IMG( Text , Property1 , PropertyValue1 , ... );
% 
%  Additional Inputs are the TextUIControlProperties and their Values,
%   like FontSize, FontWeight, FontName, BackGroundColor ... .
%  A Special Input is an OverSize in Width and High, which is sometimes nessesary
%   for large and bold Fonts. The Units of this OverSize are Pixels:
%     TXT2IMG( Text , ... , 'OverSize' , [ Width  High ] , ... )
%  
%  If you use this CData for a PushButton, give an OverSize in Width
%   and High for the PushButton, and take care on the BackGroundColor.
%
%  Note:  The Size of the Text may be various on different machines,
%          depends on the ScreenPixelsPerInch of the Root.
%         If you use the created CData on an other Machine, the Size 
%          may be differs from the FontSize in which the CData was created.
%         
%--------------------------------------------------------------------------
%  Example:
%--------------------------------------------------------------------------
%
% txt = { 'Why?'
%         'The computer did it.' };
%
% % Create 3 ImageCData
%            
% [Msg1,c1] = txt2img(txt,'fontsize',12);
%
% [Msg2,c2] = txt2img(txt,'fontsize'  ,14    , ...
%                            'fontweight','bold', ...
%                            'backgroundcolor',[1 0 0]);
%
% [Msg3,c3] = txt2img(txt,'fontsize'  ,14    , ...
%                            'fontweight','bold', ...
%                            'backgroundcolor',[1 0 0], ...
%                            'oversize'   , [ 20  20 ]);
%
% if ~isempty(Msg1) | ~isempty(Msg2) | ~isempty(Msg3)
% % An Error occured
%   fprintf([ Msg1 char(10) Msg2 char(10) Msg3 char(10) ]);
%   return
% end
%
% figure('units'   , 'pixels' , ...
%        'position', [ 100 200 200 200 ]);
%
% % Create UIControls with the Size of the CData and 20 Pixel more
%
%    uicontrol('style'    , 'pushbutton',...
%              'units'    , 'pixels' , ...
%              'position' , [ 10 140 size(c1,2)+20 size(c1,1)+20] , ...
%              'cdata'    ,  c1 , ...
%              'callback' , 'disp(''1: ok'')'  )             
%
%    uicontrol('style'    , 'pushbutton',...
%              'units'    , 'pixels' , ...
%              'position' , [ 10 80 size(c2,2)+20 size(c2,1)+20] , ...
%              'cdata'    ,  c2  , ...
%              'callback' , 'disp(''2: good'')'    )             
%
%    % c3 has allready OverSize of 20 !!!
%    % uicontrol didn't need an Oversize
%    uicontrol('style'    , 'pushbutton',...
%              'units'    , 'pixels' , ...
%              'position' , [ 10 20 size(c3,2) size(c3,1)] , ...
%              'cdata'    ,  c3  , ...
%              'callback' , 'disp(''3: better'')'    )             
%

c = [];

nl = char(10);

Msg  = '';
Msg0 = ' TXT2IMG: ';

OverSize = [ 0  0 ];


if nargin < 1
  Msg = [ Msg0  'Input Text is undefined.' ];
  return
end


if ischar(txt)
 txt = cellstr(txt);
end

if ~iscellstr(txt)
  Msg = 'Text must be a  CharArray or CellStringArray.';
end

txt = txt(:);

VarArg = varargin;
VarArg = VarArg(:);

if ( mod(size(VarArg,1),2) ~= 0 )

  Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
          'Additional Inputs must contain UIControl-Property-Value-Pairs.' ];

end

VarArg = reshape(VarArg,2,size(VarArg,1)/2)';


if ~iscellstr(VarArg(:,1))
 Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
         ' Properties must be Strings.' ];
end

if ~isempty(Msg)
  Msg = [ Msg0  Msg ];
  return
end

% Look for Oversize

ii = find(strcmp('OVERSIZE',upper(VarArg(:,1))));
if ~isempty(ii)

  OverSize = VarArg{ii,2};

  VarArg(ii,:) = [];

  OverSize = OverSize(:)';
  if ~isnumeric(OverSize) | isempty(OverSize) | size(OverSize,2) > 2
    Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
             ' OverSize must be a 2-Element Numeric: [ Width High ]' ];
  elseif ~all(isfinite(OverSize)) | any( OverSize < 0 )
    Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
             ' OverSize must contain Values >= 0 ' ];
  end  
 
  if size(OverSize,2) == 1
    OverSize = [ OverSize  0 ];
  end   

end


if ~isempty(Msg)
  Msg = [ Msg0  Msg ];
  return
end


VarArg(:,3) = { [] };  % DefaultValues


% Set DefaultProperties 

for ii = 1 : size(VarArg,1)

  msg = '';

  eval( ' VarArg{ii,3} = get(0,[''DefaultUIControl'' VarArg{ii,1}]);'  ,  ...
        ' msg = lasterr;' );

  if isempty(msg)
     set(0,['DefaultUIControl' VarArg{ii,1}],VarArg{ii,2});
  else
    Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
            'Invalid UIControlProperty: ' VarArg{ii,1} ];
    break
  end

end

if ~isempty(Msg)

  % Set DefaultProperties back
  for jj = 1 : ii-1
    set(0,['DefaultUIControl' VarArg{ii,1}],VarArg{ii,3});
  end

  Msg = [ Msg0  Msg ];
  return

end


nrow = size(txt,1);
ncol = size(char(txt),2);


fig = figure('menubar'  , 'none', ...
             'units'    , 'pixels' , ...
             'color'    ,   get(0,'DefaultUIControlBackGroundColor') , ... 
             'colormap' , [ get(0,'DefaultUIControlForeGroundColor') ; ...
                            get(0,'DefaultUIControlBackGroundColor')      ] , ...
             'visible'  , 'on' , ...
             'createfcn', ''     );

figpos = get(fig,'position');


% Create Axes and Text to Determine TextWidth

axe = axes('parent' , fig , ...
           'units','pixels' , ...
           'position',[1 1  figpos([3 4]) ] , ...
           'xlim'    , [ 0  1 ] , ...
           'ylim'    , [ 0  1 ] , ...
           'xtick'   , [] , ...
           'ytick'   , [] , ...
           'visible' , 'off'  );


ext = zeros(nrow,4);  % TextExtension

ht = text('parent' , axe , ...
          'units'   , 'pixels' , ...
          'position',[ 1 1 0 ] , ...
          'string'  , txt{1,:} , ... 
          'fontunits'  , get(0,'DefaultUIControlFontUnits') , ... 
          'fontsize'   , get(0,'DefaultUIControlFontSize')  , ... 
          'fontname'   , get(0,'DefaultUIControlFontName')  , ... 
          'fontangle'  , get(0,'DefaultUIControlFontAngle') , ... 
          'fontweight' , get(0,'DefaultUIControlFontWeight') , ...
          'visible'    , 'off' , ...
          'interpreter' , 'none'      );

ext(1,:) = get(ht,'extent');

for ii = 2 : nrow
  set(ht,'string',txt{ii});
  ext(ii,:) = get(ht,'extent');
end

delete(axe)

ext = max(ext(:,3)) + 4;   % TextWidth


% Create TextUIControl

h = uicontrol('parent'   , fig  , ...
              'visible'  , 'on' , ...
              'selected' , 'off' , ...
              'units'    , 'characters', ...
              'position' , [0 0 ncol nrow ] , ...
              'style'    , 'text' , ...
              'string'   , txt  );


%  [ characters ]  --- get(0,'defaultuicontrolfont...') --> [ pixels ]

set(h,'units','pixels');


% Set DefaultProperties back
for ii = 1 : size(VarArg,1)
  set(0,['DefaultUIControl' VarArg{ii,1}],VarArg{ii,3});
end


figpos = get(fig,'position');
   pos = get(h,'position');

   pos(3) = ext;

set(fig,'position',[ figpos([1 2])  pos([3 4])+OverSize ] );

set(h,'position',[ 1  1+OverSize(2)/2  pos([3 4])+[ OverSize(1) 0 ] ]);


figure(fig);

c = getframe(fig);

m = c.colormap;
c = c.cdata;

if isempty(m)  &  size(c,3) == 3

  c = double(c)/255;

else

  c = double(c);

  n = size(m,1);

  c(:,:,2) = c(:,:,1)+1*n;
  c(:,:,3) = c(:,:,1)+2*n;

  c = m(c);

end

delete(fig)

