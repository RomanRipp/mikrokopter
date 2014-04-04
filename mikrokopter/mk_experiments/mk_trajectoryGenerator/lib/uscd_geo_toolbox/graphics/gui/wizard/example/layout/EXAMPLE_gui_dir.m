function Config = gui_dir 

% GUI_DIR  GUI-Structure for Directory-Section (Ordner)
%
% called by: GUI_CONTROL
%


is_win = strcmp( upper(computer) , 'PCWIN' );

%--------------------------------------------------------------------------
% ButtonWidths

wwt  = 10;  % TextNames
wwl  = 40;  % DirectoryLists
wwp0 =  3;  % PushButtons
wwp1 =  6;  % PushButtons
wwp2 =  8;

wwi = 90 - 15*is_win ;  % InfoListBox

wwf = wwi + ( 8 + 2*is_win ) * i;  % FileInfo

%--------------------------------------------------------------------------
% ButtonColors

cl0 = [ 1.0  1.0   1.0 ];
cl1 = [ 0.95  1.0   1.0 ];  % PushButton

cl3 = [ 1.0  1.0   0.9 ];  % Neu
cl4 = [ 1.0  0.95   0.95 ];  % Loeschen



%--------------------------------------------------------------------------
% CallBackFunctions

CBFcn = 'EXAMPLE_clb_dir';


%*******************************************
% Structures for UIControls


str1 = ' Inhalt: ';
str2 = ' Datei: ';

%--------------------------------------------------------------------------
% Arbeitsverzeichnis

Edit = stdgui('Type' , 'edit' , 'Width' , wwl, 'Color' , cl3, 'CBFcn' , CBFcn,...
              'Option' , 'Neuer OrdnerName', ... 
              'UserData' , [ 'Neuer Ordner in:' char(10) '  ' ] );


% UserData = [ 1  0  -1 ], History

List = stdgui('Type' , 'popupmenu' , 'Width' , wwl, 'UserData' , NaN, ...
              'CBFcn' , CBFcn, 'Option' , 'Wähle aktuellen Ordner');

Down  = stdgui('Type' , 'pushbutton' , 'Width' , wwp0, 'Color' , cl1, 'String' , '<<' , ...
              'CBFcn' , CBFcn, 'Option' , 'Verlasse gewählten Ordner (zurück)');

Up    = stdgui('Type' , 'pushbutton' , 'Width' , wwp0, 'Color' , cl1, 'String' , '>>' , ...
              'CBFcn' , CBFcn, 'Option' , 'Wechsel in unten gewählten Ordner');

New   = stdgui('Type' , 'pushbutton' , 'Width' , wwp1, 'Color' , cl3 ,  ...
              'String' , 'Neu' , 'CBFcn' , CBFcn , 'Option' , ...
              'Erzeuge neuen Ordner im aktuellen Ordner');

Delete = stdgui('Type' , 'pushbutton' , 'Width' , wwp2, 'Color' , cl4, ...
               'String' , 'Löschen' , ...
               'CBFcn' , CBFcn, 'Option' , 'Lösche aktuellen Ordner');
 
 Select = struct( ...
  'Edit'   , { { 1  Edit } } , ...
  'List'   , { { 1  List } } , ...
  'Down'   , { { 2  Down } } , ...
  'Up'     , { { 3  Up   } } , ...
  'New'    , { { 4  New  } } , ...
  'Delete' , { { 5  Delete } }        );

%--------------------------------------------------------------------------

tip = [ 'DoppelClick auf Ordner wechselt in diesen, ' char(10) ...
        'DoppelClick auf Datei zeigt Inhalt.' ];

Text = stdgui('Type' , 'text' , 'Width' , wwi, 'String' , str1,...
              'Style' , 'list' , 'UserData' , str1);


% UserData = { IsDir Name }

List = stdgui('Type' , 'listbox' , 'Width' , wwi+8*i, 'UserData' , { 0  '' }, ...
              'CBFcn' , CBFcn, 'Option' , tip);

 Dir = struct( ...
  'Text'   , { { 1+0*i  Text } } , ...
  'List'   , { { 1+1*i  List } }       );

%--------------------------------------------------------------------------

Text = stdgui('Type' , 'text' , 'Width' , wwi, 'String' , str2, 'Style' , 'list' , 'UserData' , str2);

List = stdgui('Type' , 'listbox' , 'Width' , wwf);

 File = struct( ...
  'Text'   , { { 1+0*i  Text } } , ...
  'List'   , { { 1+1*i  List } }       );

%--------------------------------------------------------------------------
% Output

Config = struct( 'Select'    , { Select }   , ...
                 'Dir'       , { Dir    }   , ...
                 'Separator' , { NaN    }   , ... 
                 'File'      , { File   }         );

