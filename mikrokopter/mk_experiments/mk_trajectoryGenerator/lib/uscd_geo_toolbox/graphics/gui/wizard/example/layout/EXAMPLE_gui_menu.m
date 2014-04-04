function Config = gui_menu

% GUI_MENU  GUI-Strcuture for Menu-Bar
%
% called by: GUI_CONTROL
%




%--------------------------------------------------------------------------
% CallBackFunctions

CBWindow = '';
CBSort   = '';
CBHelp   = '';
CBQuit   = '';
CBInvert = '';

   
%*************************************************************************
% UIMenus

Window = struct( 'Geometry'  , { { 'Geometrie'  CBWindow  } } , ...
                 'Grafic'    , { { 'Grafik'     ''        } } , ...
                 'Help'      , { { 'Hilfe'      CBWindow  } } , ...
                 'Messages'  , { { 'Messages'   ''        } }       );


Sort   = struct( 'Name'  , { { 'nach Name'  CBSort } } , ...
                 'Time'  , { { 'nach Zeit'  CBSort } }       );

Color = struct( 'Invert'    , { { 'Invertieren'   CBInvert } } );

Option = struct( 'Sort'  , { { 'Sortiere Ordner und Dateien ...'  Sort  } } , ...
                 'Color' , { { 'Farben ...'                       Color } }       );


Help   = struct( 'About' , { { 'Über ...'  CBHelp } } );


Quit   = struct( 'Quit'  , { { 'Program beenden'  CBQuit  } } , ...
                 'Exit'  , { { 'Matlab beenden'   CBQuit  } }       );


Blank = '  ';


   
%*************************************************************************
% Output

Config = struct( 'Window' , { { 'Window'    Window } } , ...
                 'Option' , { { 'Optionen'  Option } } , ...
                 'Blank1' , { { Blank       ''     } } , ...
                 'Help'   , { { 'Hilfe'     Help   } } , ...
                 'Blank2' , { { Blank       ''     } } , ...
                 'Finish' , { { 'Beenden'   Quit   } }       );
 
 

