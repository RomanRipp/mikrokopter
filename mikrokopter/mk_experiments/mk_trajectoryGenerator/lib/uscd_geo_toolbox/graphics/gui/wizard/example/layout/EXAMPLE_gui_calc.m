function Config = gui_calc

% GUI_CALC  GUI-Structure for Calculation-Section (Berechnung)
%
% called by: GUI_CONTROL
%
% calls: 
%
%  GUI_SEL     Select-SubSection (Auswahl)
%  GUI_START    Start-SubSection (Start)
%


%--------------------------------------------------------------------------
% ButtonColors

cl5 = [ 0.85  0.85  0.85  ];   % TagFrame


%--------------------------------------------------------------------------
% TagFrame

 str = { 'Auswahl'  'Start' };


 opt = struct( 'Frame'     , { []          } , ...
               'Select'    , {  EXAMPLE_gui_sel    } , ...
               'Start'     , {  EXAMPLE_gui_start  }      );

 Frame = stdgui('Type' , 'tag_frame' , 'Color' , cl5, 'String' , str, 'Option' , opt);

 Tag = struct( 'Frame' , { { 1  Frame } } );

%--------------------------------------------------------------------------
% Output

Config = struct( 'Tag' , {  Tag  }   );
