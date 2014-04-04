function cnf = startcnf

% STARTCNF  Configuration for STARTMENU
%
% Returns Structure CNF with the Fields:
%
%  Name          Name of Session   
%  Label         Label of Session, 1 Word
%  Startup       Full Name of special STARTUP-Script
%  Option        Inputs for STARTUP-Script, CellArray !!!
%  WorkDir       Working Directory
%  EnhcDir       Additional Directories to add to SearchPath, 
%                  CellStringArray or CharacterArray
%  EnhcOpt       Option for additional Directories
%                   'begin' | 'end'
%
%
% see also:  STARTMENU
%
% example:
% 
% Use easily following Syntax to define the StartMenu-Items:
%
% %-------------------------------------------
%  cnf = [];   % Reset Output
%
% %-------------------------------------------
% % Local Definitions
%
%  HomePath    = getenv('HOME');
% 
%  ProjectPath = fullfile(HomePath,'project');
%  ToolboxPath = fullfile(HomePath,'matlab','toolbox');
% 
% %-------------------------------------------
% % 1. Item
%
%  n = prod(size(cnf));
%
%  cnf(n+1).Name     = 'Default';
%  cnf(n+1).Label    = 'Default';
%  cnf(n+1).Startup  = '';
%  cnf(n+1).Option   = '';
%  cnf(n+1).WorkDir  = HomePath;
%  cnf(n+1).EnhcDir  = ToolboxPath; 
%  cnf(n+1).EnhcOpt  = 'begin';
%
% %-------------------------------------------
% % 2. Item
%
%  n = prod(size(cnf));
%
%  cnf(n+1).Name     = 'Project Startup';
%  cnf(n+1).Label    = 'Project';
%  cnf(n+1).Startup  = fullfile(ProjectPath,'project_startup');
%  cnf(n+1).Option   = '';
%  cnf(n+1).WorkDir  = ProjectPath;
%  cnf(n+1).EnhcDir  = ToolboxPath; 
%  cnf(n+1).EnhcOpt  = 'begin';
%
% %-------------------------------------------
% % 3. Item
%
%  n = prod(size(cnf));
%
%   ...
%
%

%*************************************************************
% Reset Output
%-------------------------------------------------------------

cnf = [];

%*************************************************************
% Local Definitions
%-------------------------------------------------------------

ProjectPath = '/d1/project/matlab/';

HomePath    = getenv('HOME');

ToolboxPath = fullfile(HomePath,'matlab','toolbox');


%*************************************************************
% ConfigurationStructure
%-------------------------------------------------------------

n = prod(size(cnf));

cnf(n+1).Name     = 'Default';
cnf(n+1).Label    = 'Default';
cnf(n+1).Startup  = '';
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = fullfile(ProjectPath,'work');
cnf(n+1).EnhcDir  = ToolboxPath; 
cnf(n+1).EnhcOpt  = 'begin';

n = prod(size(cnf));

cnf(n+1).Name     = 'DataMan';
cnf(n+1).Label    = 'DataMan';
cnf(n+1).Startup  = '';
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = fullfile(ProjectPath,'dataman');
cnf(n+1).EnhcDir  = { ToolboxPath
                      fullfile(ProjectPath,'dataman') };
cnf(n+1).EnhcOpt  = '-begin';

n = prod(size(cnf));

cnf(n+1).Name     = 'Survey';
cnf(n+1).Label    = 'Survey';
cnf(n+1).Startup  = 'survey';
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = fullfile(ProjectPath,'instruments/matlab/instruments/survey');
cnf(n+1).EnhcDir  = { ToolboxPath };
cnf(n+1).EnhcOpt  = '-begin';

n = prod(size(cnf));

cnf(n+1).Name     = 'Topo';
cnf(n+1).Label    = 'Topo';
cnf(n+1).Startup  = '';
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = fullfile(ProjectPath,'topo');
cnf(n+1).EnhcDir  = { fullfile(ToolboxPath,'database')
                      fullfile(ToolboxPath,'graphics','map')
                      fullfile(ToolboxPath,'graphics','color')
                      fullfile(ProjectPath,'topo') };
cnf(n+1).EnhcOpt  = '-begin';

n = prod(size(cnf));

cnf(n+1).Name     = 'Map';
cnf(n+1).Label    = 'Map';
cnf(n+1).Startup  = '';
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = fullfile(ProjectPath,'map');
cnf(n+1).EnhcDir  = { ToolboxPath
                      fullfile(ProjectPath,'map') };
cnf(n+1).EnhcOpt  = '-begin';

n = prod(size(cnf));

cnf(n+1).Name     = 'GraphicsMap';
cnf(n+1).Label    = 'GraphicsMap';
cnf(n+1).Startup  = '';
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = fullfile(ProjectPath,'work');
cnf(n+1).EnhcDir  = { fullfile(ToolboxPath,'netcdf','mexcdf')
                      fullfile(ToolboxPath,'unix')
                      fullfile(ToolboxPath,'diskfun')
                      fullfile(ToolboxPath,'graphics','map')
                      fullfile(ToolboxPath,'database') };
cnf(n+1).EnhcOpt  = '-begin';

n = prod(size(cnf));

cnf(n+1).Name     = 'Start Hammelmann Planet-Simulation-GUI';
cnf(n+1).Label    = 'Hammelmann-Simulation';
cnf(n+1).Startup  = 'start_gui';
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = fullfile(ProjectPath,'hm','planet');
cnf(n+1).EnhcDir  = fullfile(ProjectPath,'hm','planet');
cnf(n+1).EnhcOpt  = '';

n = prod(size(cnf));

cnf(n+1).Name     = 'Hammelmann Drum-Simulation';
cnf(n+1).Label    = 'Hammelmann-Drum';
cnf(n+1).Startup  = '';
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = fullfile(ProjectPath,'hm','drum');
cnf(n+1).EnhcDir  = fullfile(ProjectPath,'hm','drum');
cnf(n+1).EnhcOpt  = '';

n = prod(size(cnf));

cnf(n+1).Name     = 'Cruise';
cnf(n+1).Label    = 'Cruise';
cnf(n+1).Startup  = fullfile(ProjectPath,'cruise','matlab','cruise_startup.m');
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = fullfile(ProjectPath,'cruise');
cnf(n+1).EnhcDir  = ''; 
cnf(n+1).EnhcOpt  = '';

n = prod(size(cnf));

cnf(n+1).Name     = 'Start Raytrace-GUI';
cnf(n+1).Label    = 'Raytrace';
cnf(n+1).Startup  = 'ray_start';
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = fullfile(ProjectPath,'raytrace');
cnf(n+1).EnhcDir  = fullfile(ProjectPath,'raytrace'); 
cnf(n+1).EnhcOpt  = '-begin';

n = prod(size(cnf));

cnf(n+1).Name     = 'Praktikum DataProcessing';
cnf(n+1).Label    = 'Praktikum';
cnf(n+1).Startup  = '';
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = fullfile(ProjectPath,'praktikum');
cnf(n+1).EnhcDir  = { fullfile(ProjectPath,'praktikum','tools') 
                      fullfile(ProjectPath,'praktikum','convert')
                      fullfile(ProjectPath,'praktikum','work')
                      fullfile(ToolboxPath,'data') };
cnf(n+1).EnhcOpt  = '-begin';

n = prod(size(cnf));

cnf(n+1).Name     = 'Instruments DataProcessing Development';
cnf(n+1).Label    = 'Instruments-Development';
cnf(n+1).Startup  = '';
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = fullfile(ProjectPath,'instruments');
cnf(n+1).EnhcDir  = { ToolboxPath
                      fullfile(ProjectPath,'instruments') };
cnf(n+1).EnhcOpt  = '-begin';

n = prod(size(cnf));

cnf(n+1).Name     = 'Topo Old';
cnf(n+1).Label    = 'Topo-Old';
cnf(n+1).Startup  = '';
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = fullfile(ProjectPath,'topo00');
cnf(n+1).EnhcDir  = { ToolboxPath
                      fullfile(ProjectPath,'topo00') };
cnf(n+1).EnhcOpt  = '-begin';

n = prod(size(cnf));

cnf(n+1).Name     = 'Diplom';
cnf(n+1).Label    = 'Diplom';
cnf(n+1).Startup  = '';
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = fullfile(ProjectPath,'diplom');
cnf(n+1).EnhcDir  = { ToolboxPath
                      fullfile(ProjectPath,'diplom') };
cnf(n+1).EnhcOpt  = '-begin';

n = prod(size(cnf));

cnf(n+1).Name     = 'Expo';
cnf(n+1).Label    = 'Expo';
cnf(n+1).Startup  = '';
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = fullfile(ProjectPath,'expo');
cnf(n+1).EnhcDir  = { ToolboxPath
                      fullfile(ProjectPath,'expo') };
cnf(n+1).EnhcOpt  = '-begin';

n = prod(size(cnf));

cnf(n+1).Name     = 'Thetis1';
cnf(n+1).Label    = 'Thetis1';
cnf(n+1).Startup  = '';
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = fullfile(ProjectPath,'thetis','thetis1');
cnf(n+1).EnhcDir  = { ToolboxPath
                      fullfile(ProjectPath,'map') 
                      fullfile(ProjectPath,'thetis','thetis1') };
cnf(n+1).EnhcOpt  = '-begin';

n = prod(size(cnf));

cnf(n+1).Name     = 'Instrument DataProcessing';
cnf(n+1).Label    = 'Instrument';
cnf(n+1).Startup  = '';
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = '/d0/cruise/';
cnf(n+1).EnhcDir  = { '/d0/cruise/matlab' };
cnf(n+1).EnhcOpt  = '-begin';

n = prod(size(cnf));

cnf(n+1).Name     = 'Turtle Startup';
cnf(n+1).Label    = 'Turtle';
cnf(n+1).Startup  = '';
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = fullfile(ProjectPath,'turtle');
cnf(n+1).EnhcDir  = { ToolboxPath
                       fullfile(ProjectPath,'turtle','prog')
                       fullfile(ProjectPath,'turtle','tools')
                       fullfile(ProjectPath,'turtle','nps')
                       fullfile(ProjectPath,'turtle','kalib')
                       fullfile(ProjectPath,'turtle','data')   }; 
cnf(n+1).EnhcOpt  = '-begin';

n = prod(size(cnf));

cnf(n+1).Name     = 'Start Tomolab Map';
cnf(n+1).Label    = 'Tomolab-Map';
cnf(n+1).Startup  = fullfile(ProjectPath,'tomolab','local','start','start.m');
cnf(n+1).Option   = {'map'};
cnf(n+1).WorkDir  = '';
cnf(n+1).EnhcDir  = ''; 
cnf(n+1).EnhcOpt  = '';

n = prod(size(cnf));

cnf(n+1).Name     = 'Start Tomolab Ocean-Toolbox';
cnf(n+1).Label    = 'Tomolab-Ocean';
cnf(n+1).Startup  = fullfile(ProjectPath,'tomolab','local','start','start.m');
cnf(n+1).Option   = {};
cnf(n+1).WorkDir  = '';
cnf(n+1).EnhcDir  = ''; 
cnf(n+1).EnhcOpt  = '';

n = prod(size(cnf));

cnf(n+1).Name     = 'Set Tomolab Path';
cnf(n+1).Label    = 'Tomolab-Path';
cnf(n+1).Startup  = fullfile(ProjectPath,'tomolab','local','start','start.m');
cnf(n+1).Option   = {'path'};
cnf(n+1).WorkDir  = '';
cnf(n+1).EnhcDir  = ''; 
cnf(n+1).EnhcOpt  = '';

n = prod(size(cnf));

cnf(n+1).Name     = 'Start Tomolab';
cnf(n+1).Label    = 'Tomolab';
cnf(n+1).Startup  = fullfile(ProjectPath,'tomolab','local','start','TomoStart.m');
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = '';
cnf(n+1).EnhcDir  = ''; 
cnf(n+1).EnhcOpt  = '';

n = prod(size(cnf));

cnf(n+1).Name     = 'Start DrifterAnimation';
cnf(n+1).Label    = 'Drifter';
cnf(n+1).Startup  = fullfile(ProjectPath,'drifter','start','drifter.m');;
cnf(n+1).Option   = '';
cnf(n+1).WorkDir  = '';
cnf(n+1).EnhcDir  = { ToolboxPath };
cnf(n+1).EnhcOpt  = '-begin';

n = prod(size(cnf));

cnf(n+1).Name     = 'Start Documentation';
cnf(n+1).Label    = 'Documentation';
cnf(n+1).Startup  = 'web';
cnf(n+1).Option   = '/d0/doc/index/index.html';
cnf(n+1).WorkDir  = ProjectPath;
cnf(n+1).EnhcDir  = ToolboxPath; 
cnf(n+1).EnhcOpt  = 'begin';

