function Config = gui_control 

% GUI_CONTROL  GUI-Structure for Sections
%
% called by:  GUI_CONTROL
%
% calls:  
%
%  GUI_DIR       Directory-Section (Ordner)
%  GUI_CNF   Konfiguration-Section
%  GUI_CALC    Calculation-Section (Berechnung)
%
%


%--------------------------------------------------------------------------
% Colors

cl5 = [ 0.9  0.9   0.9 ];  % TagFrame

%--------------------------------------------------------------------------

CBFcn = '';

%--------------------------------------------------------------------------
% TagFrame

 str = { 'Ordner'  'Konfiguration'   'Berechnung'  };


 opt = struct( 'Frame'       , { []         } , ...
               'Dir'         , { EXAMPLE_gui_dir    } , ...
               'Config'      , { EXAMPLE_gui_cnf    } , ...
               'Calculation' , { EXAMPLE_gui_calc   }         );

 Frame = stdgui('Type' , 'tag_frame' , 'Color' , cl5, 'String' , str, ...
                'CBFcn' , CBFcn, 'Option' , opt);

 Tag = struct( 'Frame' , { { 1  Frame } } );

%--------------------------------------------------------------------------

text1 = stdgui('Type','text','Width',60, ...
               'String','Example of an GUI, defined by EXAMPLE_GUI_MAIN, build with MKGUI');

text2 = stdgui('Type','text','Width',50, ...
               'String','Warning: Only some features work in this GUI!' );

%--------------------------------------------------------------------------
% Output

Config = struct(  'Text1' , { struct( 'Text' , { { 1  text1 } } ) } , ...
             'Separator'  , { NaN    } , ... 
                  'Text2' , { struct( 'Text' , { { 1  text2 } } ) } , ...
                  'Tag'   , {  Tag  }  );
