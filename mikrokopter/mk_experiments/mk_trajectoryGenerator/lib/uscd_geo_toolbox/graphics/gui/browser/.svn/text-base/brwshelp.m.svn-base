function [ Msg , s1 ] = brwshelp

% BRWSHELP  Example for Using BROWSER
%

 nl = char(10);  % NewLineCharacter

 EchoState = get(0,'Echo');
 EchoON    = strcmp(EchoState,'on');


 %***********************************************************
 % Wan't to open a BrowserFigure which 
 %  display's Information about the 3 M-Files:  

 Chapter = { 'browser'  'msg_logo'  'msg_list'  'textlogo' };


 %-----------------------------------------------------------
 % 3. SubItems per each Chapter

 Section = { 'About'  'Help'  'Listing' };

 % Press any Key to Continiue ...
  if EchoON, pause, end


 %-----------------------------------------------------------
 % The Default ItemStructure
  
  s0  = struct( 'Tag'      , { '' } , ...
                'Label'    , { '' } , ...
                'Titel'    , { '' } , ...
                'Text'     , { '' } , ...
                'File'     , { '' } , ...
	        'Note'     , { '' } , ...
	        'Matlab'   , { '' } , ...
                'Children' , { [] }      );

 % Press any Key to Continiue ...
  if EchoON, pause, end


 %-----------------------------------------------------------
 % The "Main"-Items (Chapter)

  n1 = size(Chapter,2);
  n2 = size(Section,2);

  s1 = s0(ones(n1,1));  % The    ItemStructure

  for ii = 1 : n1

    s1(ii).Tag   = upper(Chapter{ii});
    s1(ii).Label = upper(Chapter{ii});    

    %--------------------------------------------------------
    % The SubItems (Sections)

    s2 = s0(ones(n2,1)); % The SubItemStructure

    filename = [ Chapter{ii} '.m' ];
   
    helptext = help(filename);

    % Extract first HelpLine
          nn = find( double(helptext) == 10 );
    helptext = helptext(1:nn-1);
   
     for jj = 1 : n2

       s2(jj).Tag   = [ '__' Section{jj} ];
       s2(jj).Label = Section{jj};
       s2(jj).Titel = [ '_'  Section{jj} ];

       switch Section{jj}
         case 'About'

            s2(jj).Text = helptext;

         case 'Help'
         % Full HelpText

            s2(jj).File = ['MatlabHelp '  filename ];

            % Use 'MatlabHelp' in the File 
            %  to get the HelpText from HELP(FileName)

         case 'Listing'
         % Full Listing

            s2(jj).Text = [ nl  'Listing of '  filename  nl nl ];
            s2(jj).File = filename;

       end

     end
     % jj

     s1(ii).Children = s2;  % Add SubItem to Chapter
 
  end
  % ii


 %-----------------------------------------------------------
 % Add the TEXTLOGO-Item to MSG_LOGO

  s1(2).Children = cat( 1 , s1(2).Children , s1(4) );


 %-----------------------------------------------------------
 % HelpWin as 4. Chapter

  s1(4) = s0;

  s1(4).Tag    = 'HelpWin';
  s1(4).Label  = 'HelpWin';
  s1(4).Matlab = 'HelpWin';

  % Use 'MatlabHelpWindow' in the File to launch HELPWIN

  %----------------------------------------------------------
  % The Chapter as SubItems (Sections)

    s1(4).Children = s0(ones(n2,1)); % The SubItemStructure

     for jj = 1 : n1

       filename = [ Chapter{jj} '.m' ];

       s1(4).Children(jj).Tag    = [ '__' Chapter{jj} ];
       s1(4).Children(jj).Label  = Chapter{jj};
       s1(4).Children(jj).Titel  = [ '_'  Chapter{jj} ];
       s1(4).Children(jj).Matlab = [ 'HelpWin '  filename ];

       % Use [ 'HelpWin ' filename ] in the File 
       %  to launch HELPWIN  for  filename 
 
     end

     %-------------------------------------------------------
     % Add the TEXTLOGO-Item to MSG_LOGO

     s1(4).Children(2).Children = s1(4).Children(4);
     s1(4).Children(4)          = [];


 % Press any Key to Continiue ...
  if EchoON, pause, end


 %***********************************************************
 % Create a BrowserFigure

  Name = 'Browser M-Files';
  Logo = {'B' 'Browser' 'M-Files'};

  [ Msg , bfig , IH , IT ] = browser( 'New' , Name , Logo , s1 );

  if ~isempty(Msg)
     msgbox(Msg,'Error','warn');
     return 
  end

  set(bfig,'visible','on');

 % Play with the BrowserContents-Menu 
 % Press any Key to Continiue ...
  if EchoON, pause, end


 %-----------------------------------------------------------
 % DisActivate the Browser

  Msg = browser( 'DisActivate' , bfig );

  if ~isempty(Msg)
     msgbox(Msg,'Error','warn');
     return 
  end

 % Press any Key to Continiue ...
  if EchoON, pause, end


 %-----------------------------------------------------------
 % Activate the BROWSER-Help-Item

  Msg = browser( 'Activate' , bfig , 'BROWSER_Help' );

  if ~isempty(Msg)
     msgbox(Msg,'Error','warn');
     return 
  end

 % Press any Key to Continiue ...
  if EchoON, pause, end


 %***********************************************************
 % Create a New Figure, use MSG_LOGO

 fig = figure( 'name'    , 'MSG_LOGO' , ...
               'menubar' , 'none'     , ...
               'position' , get(bfig,'position')+[ 50 -100 0 0 ]       ); 

 msg_logo( fig , 'New' , { 'M' 'Msg_Logo' 'Example' } );


 %-----------------------------------------------------------
 % Add the MSG_LOGO-Item from the Browser, call them 'Help'

 [ Msg , hm ] = browser( 'Menu' , bfig , 'MSG_LOGO' , fig , 'Help' );

  if ~isempty(Msg)
     msgbox(Msg,'Error','warn');
     return 
  end

  figure(fig);

 % Use the New Menu to Activate a BrowserItem


