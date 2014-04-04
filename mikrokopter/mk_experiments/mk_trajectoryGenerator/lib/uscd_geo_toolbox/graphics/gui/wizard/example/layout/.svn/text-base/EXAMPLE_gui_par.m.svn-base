function Config = gui_par

% GUI_PAR  GUI-Structure for Konifguration-Parameter-Section
%
% called by: GUI_CNF
%
% calls:
%
%



is_win = strcmp( upper(computer) , 'PCWIN' );


% ButtonWidths

wwi = 32 + 6*i;  % Info Listbox

wwe = 16;        % MemoryEdit

wwt = 10;        % KennZahlenText  NameText

%--------------------------------------------------------------------------
% ButtonColors

cl0 = [ 1.0  1.0  1.0 ];
cl2 = [ 0.7  0.7  0.7 ];  % Slider

cl5 = [ 0.8  0.8  0.8 ];  % TagFrame


%--------------------------------------------------------------------------
% CallBackFunctions

CBEdt = 'EXAMPLE_clb_edt';
CBSld = 'EXAMPLE_clb_sld';


%--------------------------------------------------------------------------
% TagFrame

 str = { 'Standard'   'Oszillator'   'Planeten' };


 opt = struct( 'Frame'       , { []              } , ...
               'Standard'    , { EXAMPLE_gui_std }  );

 Frame = stdgui('Type' , 'tag_frame' , 'Color' , cl5, 'String' , str, ...
                'Option' , opt);

 Tag = struct( ...
               ...
  'Frame'  , { { 1  Frame } } );


%--------------------------------------------------------------------------
% KennZahlen

% Slider:  Format      Value  Min    Max   SliderStep

sld1 = {  '%5.2f'  [   5      1      50     2   10   ] };  % MemoryLimit


SepLeft  = stdgui('Type' , 'text' , 'String' , '');
   
TextEdit = stdgui('Type' , 'text' , 'String' , 'MemoryLimit');

SepText  = stdgui('Type' , 'text' , 'String' , '');

Edit1    = stdgui('Type' , 'edit' , 'Width' , wwe, 'Text' , ...
                   'GitterGröße  MByte' , 'CBFcn' , CBEdt );

Slider1  = stdgui('Type' , 'slider' , 'Color' , cl2 , 'Value' , sld1{2} ,  ...
                  'Userdata' , sld1{1} , 'CBFcn' , CBSld);

SepRight = stdgui('Type' , 'text' , 'String' , '');

SepFrame = stdgui('Type' , 'separator' , 'Width' , 0+i*imag(wwi) , ... 
                   'Color' , [0 0 0] , ...
                   'String' , '' , 'Style' , 'List' );

TextInfo = stdgui('Type' , 'text' , 'Width' , wwt , 'String' , 'KennZahlen' );
TextName = stdgui('Type' , 'text' , 'Width' , wwt , 'String' , '' );

    Info = stdgui('Type' , 'listbox' , 'Width' , wwi);

Par = struct( ...
  'SepLeft'    , { { 1      SepLeft  } } , ...
 'TextEdit'    , { { 2+0*i  TextEdit } } , ...
  'SepText'    , { { 2+1*i  SepText  } } , ...
     'Edit1'   , { { 2+2*i  Edit1    } } , ...
     'Slider1' , { { NaN    Slider1  } } , ...
  'SepRight1'  , { { 3      SepRight } } , ...
  'SepRight2'  , { { 4      SepRight } } , ...
  'Separator'  , { { 5      SepFrame } } , ...
  'SepInfo'    , { { 6      SepRight } } , ...
 'TextInfo'    , { { 7+0*i  TextInfo } } , ...
 'TextName'    , { { 7+1*i  TextName } } , ...
     'Info'    , { { 8      Info     } }  );


%--------------------------------------------------------------------------
% Output

Config = struct( 'Tag'    , {  Tag   }, ...
                 'Memory' , {  Par   }        );

