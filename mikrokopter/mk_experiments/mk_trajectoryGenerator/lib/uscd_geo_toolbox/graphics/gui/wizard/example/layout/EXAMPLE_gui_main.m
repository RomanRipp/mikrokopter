function Config = gui_main 

% GUI_MAIN  GUI-Structure 
%
% called by: START_GUI
%
% calls:
%
%  GUI_MENU
%  GUI_LOGO
%  GUI_CONTROL
%


Config = struct( 'Menu'    , { EXAMPLE_gui_menu    } , ...
                 'Logo'    , { EXAMPLE_gui_logo    } , ...
                 'Control' , { EXAMPLE_gui_control }        );
