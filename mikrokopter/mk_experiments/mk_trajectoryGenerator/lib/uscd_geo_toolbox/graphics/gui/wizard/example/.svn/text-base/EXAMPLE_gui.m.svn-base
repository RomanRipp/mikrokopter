function [Msg,fig] = EXAMPLE_gui

% EXAMPLE_GUI   Example to create a GUI, using MAKE_GUI
%
% [Msg,FigureHandle] = EXAMPLE_GUI
%
% type: 
%
%  >> [Msg,fig] = EXAMPLE_gui;
%
%
% calls:  EXAMPLE_GUI_MAIN  to create the GUI-ConfigurationStructure
%         MAKE_GUI          to create the GUI
%


 Config = EXAMPLE_gui_main;  % GUI-ConfigurationStructure

 [Msg,fig] = mkgui( Config , 'GUI-Example' );

 if ~isempty(Msg)
    return
 end

%*******************************************************
% Disable TagButtons

try

  ud = get(fig,'userdata');

 par = ud.Children.Frame;

 fud = get( par , 'userdata' );


  hf = fud.Children.Control.Tag.Frame;
 
  tag_frame( hf , 'Enable' , 'off' );

end
%*******************************************************


 try
     set( fig , 'visible' , 'on' , ...
       'handlevisibility' , 'on' );

 catch
     return
 end


 par_dir  = fud.Children.Control.Tag.Dir;          %   Directory-Tag-Button
 par_cfg  = fud.Children.Control.Tag.Config;       %      Config-Tag-Button
 par_calc = fud.Children.Control.Tag.Calculation;  % Calculation-Tag-Button


%*****************************************************
% Setup for Directory-Tag (Ordner)

par = par_dir;   % Tag-Button == Parent for Directory-Section


  pfad = getenv('HOME');
  if isempty(pfad)
      pfad = pwd;
  end

  SortMode = 'a';  % alphabetical

%--------------------------------------------------------

ud = get( par , 'userdata' );

set( ud.Children.Select.Edit , 'visible' , 'off' );


%--------------------------------------------------------
ud.SortMode  = SortMode;
ud.Directory = '';

set( par , 'userdata' , ud );

%--------------------------------------------------------
% "1" ==> call of NEW_DIR,
%   updates Directory in Konfiguration-, CalculationSection



 % Change to one DOWN from pfad

   p = pfad;

  sp = size(p,2);

  if sp > 1

     fs = filesep;

     ind = ( 1 : sp-1*strcmp(p(sp),fs) );
      ii = find( double(p(ind)) ==  double(fs) );
                             
      if ~isempty(ii)

         p = p(1:max(ii));

      end

  end

  if strcmp(p,pfad);
    pfad = '';
  end

  EXAMPLE_clb_dir( par , 'Select' , 'NewDir' , 1 , p , pfad );


%*****************************************************
% Setup for Configuration-Geometry-Tabular

par = par_cfg;   % Tag-Button == Parent for Section

ud = get(par,'userdata');

par1 = ud.Children.Tag.Geometry;

ud1  = get( par1 , 'userdata' );

str = { '100 60'  '20  0  0.1'  '1.5'  'rund' };

Msg = tab_list(ud1.Children.Duesen.Tabular,'string',str);


%*****************************************************
% Setup for Calculation-Start-SelectList

par = par_calc;   % Tag-Button == Parent for Section

ud = get(par,'userdata');

par1 = ud.Children.Tag.Select;

ud1  = get( par1 , 'userdata' );

% Get Contents of Directory: "toolbox/gui/example/layout"

pfad = fileparts(which('EXAMPLE_gui_main'));

c = dircont(pfad,'a','*.m');

Msg = sel_list( ud1.Children.File.Select , 'Set' , c(:,4) );

%*******************************************************
% Enable TagButtons
 
try

  tag_frame( hf , 'Enable' , 'on' );

end

%*******************************************************
