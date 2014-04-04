% DataBase Toolbox # Contents of "toolbox/data/db"
%   
% The DataBase Toolbox provides an easier access to DataFiles
% in the RODB-Format, which are DataFiles with Header.
% 
% -----------------------------------------------------------
% Syntax for a Datafile:
% 
% ROOT/BASE/NAME/TYPE/NAME_???.###
% 
% - ROOT is the Main-Directory of the DataBase,
% - BASE is a Name of a Experiment or Region,
% - NAME is a Name of a Mooring or Cruise,
% - TYPE is a Name of an Instrument or DataType
% 
% The Name of each Datafile should start with NAME,
% the DataFiles of an TYPE should be numbered ("???").
% The Extension of a DatafileName should be equal to TYPE.
% 
% The Datafiles has a Header with informations about
% the dataset in the Structure:
%  
%    Variable1 = Value1
%    Variable2 = Value2
% 
% The last HeaderVariable should define the DataColumns,
% with the Variable "Columns", the VariableNames in the Value
% are separated by the ColonCharacter ":", example:
% 
%    Columns = P:T:S
% 
% VariableNames for Header and Data are defined in DBINIT.
% Other Variables are allowed.
% 
% The Data follows the Header in individiualy Columns.
% 
% -----------------------------------------------------------
% InfoFile:
% 
% Each Cruise or Mooring should have an InfoFile like:
% 
% ROOT/BASE/NAME/NAMEinfo.dat
% 
% The HeaderVariables of the InfoFile should specify
% the Project, Location, Date and Time, etc. 
% 
% The DataColumns for Moorings should be
% the InstrumentTypes, SerialNumbers and Depth, 
% for Cruises the Stations with Position and Date and Time.
% 
% -----------------------------------------------------------
% Inquire Data:
% 
% DataBase-Directory (ROOT) should be set by DBROOT.
% Inquire or modify other settings by DBSET.
% 
% Use DBFILE to check which Datafiles exist in the DataBase-Directory,
% and DBLOAD to load Data.
%   
%------------------------------------------------------------------------
%   
%   dbcheck    - Checks Files for Header- and DataVariables
%   dbfile     - Search for DataFiles
%   dbget      - Get Configuration for DataBase-Access
%   dbident    - Decode database variable keywords
%   dbinit     - Initialisation for Variable-KeyWords
%   dbload     - Load DataBase files
%   dbmsg      - Display Messages in Command or LogFile
%   dbread     - Read Variables from DataBaseFile
%   dbroot     - Set RootPath for DataBaseAccess
%   dbsave     - H2' , H1 , 'H2' , H2 , ... , 'D1' , D1 , 'D2' , D2 , ..
%   dbscan     - Check for DataBaseFiles and their Variables
%   dbset      - Set Configuration for DataBase-Access
%   dbtest     - Inquire InfoFiles (Mooring- and CruiseInfo)
%   head2val   - Converts between HeaderValue and HeaderString
%   splithead  - Splits a HeaderString in KeyWord / Value
%   splitunit  - Splits a HeaderValue in Value / Unit
%   
