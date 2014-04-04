function Config = gui_geo

% GUI_GEO  GUI-Structure for Konifguration-Geometry-Section (Geometrie)
%
% called by: GUI_CNF
%


is_win = strcmp( upper(computer) , 'PCWIN' );

%--------------------------------------------------------------------------
% ButtonColors

cl0 = [ 1.0  1.0  1.0 ];
cl1 = [ 0.95  1.0   1.0 ]; 




%--------------------------------------------------------------------------
% CallBackFunctions

CBFcn = '';


%*******************************************
% Structures for UIControls


%--------------------------------------------------------------------------
% Tabular
 
str = ' Düsen-Konfiguration: ';


if is_win
    
  t1 = { 'Planet'
         'Radius  Winkel' };

  t2 = { 'Düse'
         'Radius Winkel   ø  Breite' };

  t3 = { 'Auff.-'
         'Faktor'  };

  t4 = { 'Rund / Flach'
         'StrahlDüse'    };

  f1 = ' %6.2f %6.1f ';
  f2 = ' %6.2f %6.1f %5.2f %5.2f ';
  f3 = '  %3.1f  ';
  f4 = 'char';


  ncol = [ 15  28  7  12 ];
      
else
    
  t1 = { 'Planet'
         'Radius   Winkel' };

  t2 = { 'Düse'
         ' Radius   Winkel    ø   Breite' };

  t3 = { 'Auff.-'
         'Faktor'  };

  t4 = { 'Rund / Flach'
         'StrahlDüse'    };

  f1 = ' %7.2f  %6.1f ';
  f2 = ' %7.2f  %6.1f  %5.2f  %5.2f ';
  f3 = '  %3.1f  ';
  f4 = 'char';

  ncol = [ 18  30  7  14 ];
      
end

 opt = {  ...
          ...  
  'ntab'  , size(ncol,2) , ...
  'nrow'  ,   8   , ...
  'ncol'  , ncol  , ...
  'title' ,  { t1  t2  t3  t4 }  , ...
  'form'  ,  { f1  f2  f3  f4 }  , ...
  'vedit' ,  { 'on'  'on'  'on'  'off'  } , ...
  'cedit' ,  'on'       , ...
  'cset'  ,  'on'       , ...
  'cins'  ,  'on'       , ...
  'cdel'  ,  'on'       , ...
  'creset',  'on'       , ...
  'cundo' ,  'off'      , ...
  'credo' ,  'off'      , ...
  'min'   ,   0         , ...
  'max'   ,   2         , ...
 'msgfcn' , 'main_dlg'        };


 Text = stdgui('Type' , 'text' , 'Width' , sum(ncol), ...
               'String' , str , 'Style' , 'list' , 'UserData' , str );

 List = stdgui('Type' , 'tab_list' , 'Color' , cl1, 'CBFcn' , CBFcn, 'Option' , opt);


% UserData = [ PlanetRadius  PlanetWinkel ...
%              DuesenRadius  DuesenWinkel DuesenDurchm DuesenBreite  Faktor ]

 Sort  = stdgui('Type' , 'pushbutton' , 'Color' , cl1, 'String' , 'Sort' , ...
               'UserData' , zeros(0,7), 'CBFcn' , CBFcn, 'Option' , 'Sortiere Liste');

 

Tab = struct( ...
  'Text'    , { { 1+0*i  Text } } , ...
  'Tabular' , { { 1+1*i  List } } , ...
  'Sort'    , { { NaN    Sort } }        );



%--------------------------------------------------------------------------
% Output

Config = struct( 'Duesen' , {  Tab  }  );

