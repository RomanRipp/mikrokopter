function  brwsini

% BRWSINI  Example for a BROWSER.INI-File of a Directory
%
%% Example for a BROWSER.INI-File of a Directory
%% =============================================
%%
%% BROWSER.INI defines an ItemStructure with Text,
%%              Files or SubDirectories.
%% 
%% All Items are Numbered, using '%3.3d'
%%
%%    Item 000       == Item  for this Directory
%% SubItem 001, ...  == Items for Children, 
%%                       Text, File or SubDirectory
%%
%% Note:  - If Item 000 exists, a single Item 
%%           with the Definitions from 000 will created.
%%          For each Children ( 001, ... ) a SubItem 
%%           will created.
%%
%%        - If Item 000 doesn't exist, a SubItem
%%           for each Children will created for
%%           the Item, which is defined in BROWSER.INI
%%           of the ParentDirectory.      
%%
%%
%% -----------------------------------------------------
%% Syntax: 
%%
%%###.<Field> :  <String>
%%
%%
%% -----------------------------------------------------
%% Fields for 000 or File/Text-SubItem
%%
%%###.Tag  :
%%###.Label: 
%%###.Titel:
%%###.Text :
%%###.File :
%%###.Note :
%%###.Matlab :
%%
%% Note:  
%%
%%  - Any of "Tag", "Label" or "Title"  must be defined and NOT empty.
%%
%%  - Any of "Text" or "File" of an SubItem ( 001, ... )
%%     must be defined and NOT empty.
%%
%%  - If the "Tag"   starts with '_', it will concatinated recursivly
%%     with the Tag  of the Parent.
%%
%%  - If the "Titel" starts with '_', it will concatinated
%%     with the Titel of the Parent, separated with ', '.
%%
%%  - The TAB-Character (char(9)) in the Text of an Item 
%%     (from "Text" and/or "File") will replaced with 8 Blanks.
%%
%%  - The '\n'-Characters in "Text" will replaced
%%      with the NewLine-Character (char(10)).
%%
%%  - If the "File" starts with 'MatlabHelpText', the Text of the Item will
%%      generated using HELP, the following String of "File" is the Input
%%      for HELP, example:  
%%       'MatlabHelpText  browser'  gets the Text from help('browser') and
%%       displays the Help for this File (BROWSER.M)
%%             
%%  - If the "File" starts with 'MatlabHelpWindow', the Matlab-Help-Window 
%%     will launched, using HELPWIN, the following String of "File" is the Input
%%     for HELPWIN, example:  
%%      'MatlabHelpWindow  browser'  calls helpwin('browser')
%%
%% -----------------------------------------------------
%% Fields for Directory-SubItem:
%%
%%###.Directory:
%%
%% Note:  The Directory must contains an own BROWSER.INI.
%%
%%
%% -----------------------------------------------------
%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%
%% Note: No Blanks at the Begin of a DefinitionLine !!!
%%
%
%
%% -------------------------
%% Item "Test"
%
%000.Tag  : Test
%000.Label: Test
%000.Titel: A Test
%000.Text : This is a Test
%
%
%% -------------------------
%% Definition for 4 SubItems
%
%%  --- File --
% 
%001.Tag  :  1
%001.Label: Test1
%001.Titel: Test for 1
%001.Text : Test: test1.txt
%001.File : test1.txt 
%
%% 'test1.txt' is an File in the Directory
%% "Text" can be undefined.
%
%%  --- Text --
% 
%002.Tag  :  2
%002.Label: Test2
%002.Titel: Test for 2
%002.Text : Test2, only Text
%
%
%%  --- HelpText --
%
%003.Tag  : Help
%003.Label: Help
%003.Titel: Matlab Help
%003.Text : 'Matlab Help'
%003.File : MatlabHelp
% 
%
%%  --- HelpWindow --
%
%004.Tag  : Help
%004.Label: HelpWindow
%004.Titel: Matlab HelpWindow
%004.Matlab: MatlabHelpWin
% 
%
%%  --- 2 Directories --
% 
%005.Directory: test3
%006.Directory: test4
%
%% 'test3' and 'test4'  are SubDirectories
%
