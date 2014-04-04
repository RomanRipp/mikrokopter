function Config = gui_sel

% GUI_SEL  GUI-Structure for Select-SubSection (Auswahl)
%
% called by:  GUI_CALC   Calculation-Select-Section
%             GUI_OPT       Optimize-Select-Section


is_win = strcmp( upper(computer) , 'PCWIN' );


wwl = 30 -  5*is_win;  % SelectList
wwi = 50 - 10*is_win;  % InfoList

nrow = 15 + 2 * is_win;

wwt = wwi;

wwi =  wwi + ( nrow + 3 + 1*is_win )*i;

%--------------------------------------------------------------------------
% ButtonColors

cl0 = [ 1.0  1.0  1.0 ];

%--------------------------------------------------------------------------
% CallBackFunctions

CBFcn = '';

%--------------------------------------------------------------------------
% FileList

 str1 = ' Wähle Datei';
 str2 = ' Datei: ';


 opt = { 'RowNumber'  , nrow  , ...
         'ColNumber'  , wwl          };

% UserData = { Marked Name FullFile };

TextList = stdgui('Type' , 'text' , 'Width' , wwl, 'String' , str1, ...
                  'Style' , 'list' , 'UserData' , cell(0,3));
    List = stdgui('Type' , 'sel_list' , 'CBFcn' , CBFcn, 'Option' , opt);

TextInfo = stdgui('Type' , 'text' , 'Width' , wwt, 'String' , str2, ...
                  'Style' , 'list' , 'UserData' , str2);
    Info = stdgui('Type' , 'listbox' , 'Width' , wwi, 'UserData' , cell(0,2));

% UserData = { Name FullFile }


File = struct( ...
               ...
'TextSelect' , { { 1+0*i  TextList } } , ...
    'Select' , { { 1+1*i      List } } , ...
'TextInfo'   , { { 2+0*i  TextInfo } } , ...
    'Info'   , { { 2+1*i      Info } }       );


%--------------------------------------------------------------------------
% Output

Config = struct( 'File' , {  File  }   );
