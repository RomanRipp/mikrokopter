function Config = gui_info

% GUI_INFO  GUI-Structure for Konifguration-Info-Section
%
% called by: GUI_CNF
%

is_win = strcmp( upper(computer) , 'PCWIN' );

%--------------------------------------------------------------------------
% ButtonWidths

wwl  =  40 - 5*is_win;

wwi = wwl + ( 12 + 2*is_win ) * i;

%--------------------------------------------------------------------------
% ButtonColors

cl0 = [ 1.0  1.0   1.0 ];


%--------------------------------------------------------------------------

str1 = ' Gewählte Datei: ';
str2 = ' Aktuelle Konfiguration: ';


%--------------------------------------------------------------------------
% SelectedFile & ActualFile

TextFile = stdgui('Type' , 'text' , 'Width' , wwl, 'String' , str1, ...
                  'Style' , 'list' , 'Userdata' , str1);
    File = stdgui('Type' , 'listbox' , 'Width' , wwi);

TextActual = stdgui('Type' , 'text' , 'Width' , wwl, 'String' , str2, ... 
                    'Style' , 'list' , 'UserData' , str2);
    Actual = stdgui('Type' , 'listbox' , 'Width' , wwi);


List = struct( ...
               ...
'TextFile'   , { { 1+0*i  TextFile   } } , ...
    'File'   , { { 1+1*i      File   } } , ...
'TextActual' , { { 2+0*i  TextActual } } , ...
    'Actual' , { { 2+1*i      Actual } }       );


%--------------------------------------------------------------------------
% Output

Config = struct( 'List'   , {  List   }  );

