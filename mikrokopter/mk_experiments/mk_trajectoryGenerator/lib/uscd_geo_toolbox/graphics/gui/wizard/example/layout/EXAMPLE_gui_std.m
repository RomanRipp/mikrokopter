function Config = gui_std

% GUI_STD  GUI-Structure for Konfiguration-Parameter-Standard-Section
%
% called by: GUI_PAR


is_win = strcmp( upper(computer) , 'PCWIN' );


% ButtonWidths

wwe  = 14;           % ParameterEdit
wwo  = 32+6*is_win;  % Options PopUp


%--------------------------------------------------------------------------
% ButtonColors

cl0 = [ 1.0  1.0  1.0 ];
cl2 = [ 0.7  0.7  0.7 ];  % Slider


%--------------------------------------------------------------------------
% CallBackFunctions


CBEdt  = 'EXAMPLE_clb_edt';
CBSld  = 'EXAMPLE_clb_sld';
CBOpt  = '';


%--------------------------------------------------------------------------

% Slider:  Format      Value  Min    Max   SliderStep

sld1 = {  '%5.0f'  [ 1000  -5000  5000   100  500    ] };  % DrehZahl
sld2 = {  '%5.2f'  [ 2.5    0.1    100   0.1    1    ] };  % Vorschub
sld3 = {  '%4.2f'  [ 0.1    0.1     10   0.1    0.5  ] };  % GitterInkrement
sld4 = {  '%3.1f'  [ 1.0    0.1      5   0.1    0.5  ] };  % AuffaecherungsFaktor

% Option

str  =  { ' GitterInkrement manuell einstellen'
          ' GitterInkrement =  1  x min. DüsenDurchmesser'
          ' GitterInkrement = 1/2 x min. DüsenDurchmesser'   }; 

usd = sld3{2}(2) * [ 0  2  1 ];

%--------------------------------------------------------------------------


Edit1   = stdgui('Type' , 'edit' , 'Width' , wwe, 'Text' , 'DrehZahl  U/min' , 'CBFcn' , CBEdt);
Slider1 = stdgui('Type' , 'slider' , 'Color' , cl2, 'Value' , sld1{2}, 'Userdata' , sld1{1}, ...
                 'CBFcn' , CBSld);

Edit2   = stdgui('Type' , 'edit' , 'Width' , wwe, 'Text' , 'Vorschub  m/min' , 'CBFcn' , CBEdt);
Slider2 = stdgui('Type' , 'slider' , 'Color' , cl2, 'Value' , sld2{2}, 'Userdata' , sld2{1}, ...
                 'CBFcn' , CBSld);
 
Edit3   = stdgui('Type' , 'edit' , 'Width' , wwe, 'Text' , 'GitterInkrement  mm' , 'CBFcn' , CBEdt);
Slider3 = stdgui('Type' , 'slider' , 'Color' , cl2, 'Value' , sld3{2}, 'Userdata' , sld3{1}, ...
                 'CBFcn' , CBSld);
 
Edit4   = stdgui('Type' , 'edit' , 'Width' , wwe+2, 'Text' , 'AuffächerungsFaktor' , 'CBFcn' , CBEdt);
Slider4 = stdgui('Type' , 'slider' , 'Color' , cl2, 'Value' , sld4{2}, 'Userdata' , sld4{1}, ...
                 'CBFcn' , CBSld);
 
Par = struct( ...
 'Edit1'   , { { 1    Edit1   } } , ...
 'Slider1' , { { NaN  Slider1 } } , ...
 'Edit2'   , { { 2    Edit2   } } , ...
 'Slider2' , { { NaN  Slider2 } } , ...
 'Edit3'   , { { 3    Edit3   } } , ...
 'Slider3' , { { NaN  Slider3 } } , ...
 'Edit4'   , { { 4    Edit4   } } , ...
 'Slider4' , { { NaN  Slider4 } } );


%--------------------------------------------------------------------------

Text = stdgui( 'Type' , 'text' , 'String' , 'Option:' );
List = stdgui( 'Type' , 'popupmenu' , 'Width' , wwo , 'String' , str , ...
               'Value' , 2 , 'UserData' , usd , 'CBFcn' , CBOpt );

Opt = struct(  ...
  'Text'   , { { 1  Text } } , ...
  'List'   , { { 2  List } }        );


%--------------------------------------------------------------------------
% Output

Config = struct( 'Standard'  , {  Par } , ...
                 'Option'    , {  Opt }         );

