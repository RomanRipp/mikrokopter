function Config = gui_cnf

% GUI_CNF  GUI-Structure for Konfiguration-Section
%
% called by: GUI_CONTROL
%
% calls:
%
%  GUI_GEO      Geometry-SubSection
%  GUI_PAR     Parameter-SubSection
%  GUI_INFO         Info-SubSection
%

% ButtonWidths

wwt  = 10;  % TextNames
wwf  = 24;  % FileLists
wwp1 =  6;  % PushButtons
wwp2 =  8;  % PushButtons
wwp3 = 10;  % PushButtons
wwl  = 32;  % Commentary


wwe = 12;  % Parameter
wwte = 10;

%--------------------------------------------------------------------------
% ButtonColors

cl0 = [ 1.0  1.0   1.0 ];
cl1 = [ 0.95  1.0   1.0 ];  % PushButton

cl3 = [ 1.0   1.0    0.9   ];   % Neu
cl4 = [ 1.0   0.95   0.95  ];   % Loeschen
cl5 = [ 0.85  0.85   0.85  ];   % TagFrame


%--------------------------------------------------------------------------
% CallBackFunctions

CBFile      = '';
CBComm      = '';
CBFrame     = '';


%--------------------------------------------------------------------------
% File


Text = stdgui('Type' , 'text' , 'String' , 'Datei');

Edit = stdgui('Type' , 'edit' , 'Width' , wwf, 'Color' , cl3, 'CBFcn' , CBFile,...
              'Option' , 'Neuer DateiName');

% UserData =  { Exist Name FullFile }
List = stdgui('Type' , 'popupmenu' , 'Width' , wwf, 'String' , {'Wähle Datei ...'} , ...
              'UserData' , cell(0,3), 'CBFcn' , CBFile,...
              'Option' , 'Wähle Konfigurations-Datei');

New  = stdgui('Type' , 'pushbutton' , 'Width' , wwp1, 'Color' , cl3, ...
              'String' , 'Neu' , ...
              'CBFcn' , CBFile, 'Option' , 'Neuer DateiName');

Load = stdgui('Type' , 'pushbutton' , 'Width' , wwp1, 'Color' , cl1, ...
              'String' , 'Laden' , ...
              'CBFcn' , CBFile, 'Option' , 'Lade gewählte Konfiguration');

Save = stdgui('Type' , 'pushbutton' , 'Width' , wwp2, 'Color' , cl1, ...
              'String' , 'Speichern' , ...
              'CBFcn' , CBFile, 'Option' , 'Sichere aktuelle Konfiguration');

Delete = stdgui('Type' , 'pushbutton' , 'Width' , wwp2, 'Color' , cl4, ...
                'String' , 'Löschen' , ...
               'CBFcn' , CBFile, 'Option' , 'Lösche gewählte Datei');

File = struct( ...
  'Text'   , { { 1  Text } } , ...
  'Edit'   , { { 2  Edit } } , ...
  'List'   , { { 2  List } } , ...
  'New'    , { { 3  New  } } , ...
  'Load'   , { { 4  Load } } , ...
  'Save'   , { { 5  Save } } , ...
  'Delete' , { { 6  Delete } }        );

%--------------------------------------------------------------------------
% Commentary

Text = stdgui('Type' , 'text' , 'String' , 'Kommentar');


% UserData = OldString

Edit = stdgui('Type' , 'edit' , 'Width' , wwl+2*i, 'UserData' , {''}, ...
              'Option' , 'Kommentar');


Apply = stdgui('Type' , 'pushbutton' , 'Width' , wwp3, 'Color' , cl1, ...
               'String' , 'Übernehmen' ,'CBFcn' , CBComm, ...
               'Option' , 'Übernehme Kommentar in aktuelle Konfiguration');

Reset = stdgui('Type' , 'pushbutton' , 'Width' , wwp2, 'Color' , cl4, ...
               'String' , 'Löschen' ,  'CBFcn' , CBComm, ...
               'Option' , 'Lösche Kommentar');

Comm = struct( ...
               ...
  'Text'   , { { 1  Text } } , ...
  'Edit'   , { { 2  Edit } } , ...
  'Apply'  , { { 3  Apply } } , ...
  'Reset'  , { { 4  Reset } }       );


%--------------------------------------------------------------------------
% TagFrame

 str = { 'Geometrie'  'Parameter'   'Information' };


 opt = struct( 'Frame'       , { []       } , ...
               'Geometry'    , { EXAMPLE_gui_geo  } , ...
               'Parameter'   , { EXAMPLE_gui_par  } , ...
               'Info'        , { EXAMPLE_gui_info  }      );


 Frame = stdgui('Type' , 'tag_frame' , 'Color' , cl5, 'String' , str, ...
                'CBFcn' , CBFrame, 'Option' , opt);

 Tag = struct( 'Frame' , { { 1  Frame } } );

%--------------------------------------------------------------------------
% Output

Config = struct( 'File'       , { File   } , ...
                 'Separator'  , { NaN    } , ... 
                 'Kommentar'  , { Comm   } , ...
                 'Tag'        , { Tag    }   );

