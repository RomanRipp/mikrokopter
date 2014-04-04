function [nkey,nfrm,nlen,dbkey] = dbinit

% DBINIT  Initialisation for Variable-KeyWords
%
% [ NName , NFormat , NLength , Init ] = DBINIT
%
% Init: 'Name        Format   Length  Code'
%
%  a negative or ZERO Value for Code marks a HeaderVariable
%
%  a negative Value for Length marks a SeperatorCharacter!
%    use 99-CODE for Values ASCII-Codes equal or larger 100 !!!
%

nkey = 24;   % Length of KeyWord, See below !!!
nfrm = 12;   % Length of FormatString
nlen = 03;   % Length of LengthString

if nargout < 4
   return
end

% keyword array (first 24 characters) , its C-format (12 characters),
% the number of cols, which are necessary for the information and the
% keyword identifying number, negative values indicate header keywords
%
% For Keywords with Format "list" the LengthString 
%  is the negative ASCII-Code of the separator
%

dbkey = [ ...
%|<---      24      --->||<---12--->||3|            
'Type                    %s            0      0' ; ...  % Type of file
'Filename                %s            0     -1' ; ...  % (raw data) file name
'Shipname                %s            0     -2' ; ...  % ship name
'Cruise                  %s            0     -3' ; ...  % cruise identifier
'Mooring                 %s            0     -4' ; ...  % mooring identifer
'Model                   %s            0     -5' ; ...  % model identifer
'Project                 %s            0     -6' ; ...  % project(s)
'Experiment              %s            0     -7' ; ...  % name(s) of experiment
'Comment                 %s            0     -8' ; ...  % comments
'Author                  %s            0     -9' ; ...  % author of file
'Instrument              %s            0    -10' ; ...  % instrument identifier
'InstrumentType          %s            0    -10' ; ...  % instrument identifier
'InstrType               %s            0    -10' ; ...  % instrument identifier
'Serial_Number           %d            1    -11' ; ...  % serial number
'SerialNumber            %d            1    -11' ; ...  % serial number
'Station                 %s            0    -12' ; ...  % station number or name
'Profile                 %s            0    -13' ; ...  % profile identifier
'ProfileNumber           %d            1    -14' ; ...  % profile number
'BinNumber               %d            1    -15' ; ...  % bin number (adcp)
'CTD_File                %s            0    -20' ; ...  % CTD file name
'Meteo_File              %s            0    -21' ; ...  % meterological data file name
'WaterTemp               %g            1    -31' ; ...  % water temp (deg C)
'WaterDepth              %d            1    -32' ; ...  % water depth (m)
'Max_Depth               %d            1    -33' ; ...  % maximum profile depth (m)
'MaxDepth                %d            1    -34' ; ...  % maximum profile depth (m)
'Max_Press               %d            1    -41' ; ...  % maximum profile pressure (dbar)
'MaxPress                %d            1    -41' ; ...  % maximum profile pressure (dbar)
'InstrumentDepth         %g            1    -35' ; ...  % instrument depth (m)
'InstrDepth              %g            1    -35' ; ...  % instrument depth (m)
'Depth                   %g            1    -36' ; ...  % (instrument) depth (m)
'SensorSpacing           %g            1    -37' ; ...  % sensor spacing (m) (thermistor strings)
'LayerThickness          %g            1    -38' ; ...  % 
'Press_Offset            %g            1    -39' ; ...  % pressure offset (dbar)
'PressOffset             %g            1    -39' ; ...  % pressure offset (dbar)
'MagDev                  lon           1    -40' ; ...
'Mag_Dev                 lon           1    -40' ; ...
'Mag_Deviation           lon           1    -40' ; ...  % Magnetic Deviation
'MagDeviation            lon           1    -40' ; ...
'MagneticDeviation       lon           1    -40' ; ...
'MagneticDev             lon           1    -40' ; ...
'DepthStart              %g            1    -43' ; ...  % Bin 1 mid  (m) (adcp)
'DepthIntervall          %g            1    -44' ; ...  % Bin length (m) (adcp)
'Orientation             %s            0    -45' ; ...  % Orientation (up/down ...)
'BinDepth                %g            0    -46' ; ...  % DepthValues   of Bin
'SensorDepth             %g            0    -46' ; ...  % DepthValues   of Sensors
'SensorSerialNumber      %d            0    -47' ; ...  % SerialNumbers of Sensors
'LatentHeat              %g            1    -50' ; ...
'ThermalExpansion        %g            1    -51' ; ...
'HalineContraction       %g            1    -52' ; ...
'ReferenceTemperature    %g            1    -53' ; ...
'RefTemp                 %g            1    -53' ; ...
'ReferenceSalinity       %g            1    -54' ; ...
'RefSal                  %g            1    -54' ; ...
'ReferenceDensity        %g            1    -55' ; ...
'RefDens                 %g            1    -55' ; ...
'ReferencePressure       %g            1    -56' ; ...
'RefPress                %g            1    -56' ; ...
'Temp_Coef               %s            0    -60' ; ...
'TempCoef                %s            0    -60' ; ...
'TempCorrCoef            %s            0    -60' ; ...
'TempCorr                %s            0    -60' ; ...
'TempCorrMode            %s            0    -61' ; ...
'Cond_Coef               %s            0    -62' ; ...
'CondCoef                %s            0    -62' ; ...
'CondCorrCoef            %s            0    -62' ; ...
'CondCorr                %s            0    -62' ; ...
'CondCorrMode            %s            0    -63' ; ...
'Press_Coef              %s            0    -64' ; ...
'PressCoef               %s            0    -64' ; ...
'PressCorrCoef           %s            0    -64' ; ...
'PressCorr               %s            0    -64' ; ...
'PressCorrMode           %s            0    -65' ; ...
'Time_Coef               %s            0    -66' ; ...
'TimeCoef                %s            0    -66' ; ...
'TimeCorrCoef            %s            0    -66' ; ...
'TimeCorr                %s            0    -66' ; ...
'TimeCorrMode            %s            0    -67' ; ...
'Gravity                 %g            1    -70' ; ...
'Start_Date              date          3    -80' ; ...  % start date
'StartDate               date          3    -80' ; ...  % start date
'Date                    date          3    -80' ; ...  % (start) date
'Start_Time              time          2    -81' ; ...  % start time
'StartTime               time          2    -81' ; ...  % start time
'Time                    time          2    -81' ; ...  % (start) time
'End_Date                date          3    -82' ; ...  % end date
'EndDate                 date          3    -82' ; ...  % end date
'End_Time                time          2    -83' ; ...  % end time
'EndTime                 time          2    -83' ; ...  % end time
'Time_Step               %g            1    -84' ; ...
'TimeStep                %g            1    -84' ; ...
'SamplingInterval        %g            1    -84' ; ...
'Base_Date               date          2    -86' ; ...  % BaseDate for decimal Days
'BaseDate                date          2    -86' ; ...  % BaseDate for decimal Days
'Base                    date          2    -86' ; ...  % BaseDate for decimal Days
'Start_Lat               lat           1    -90' ; ...  % start latitude
'StartLat                lat           1    -90' ; ...  % start latitude
'Latitude                lat           1    -90' ; ...  % (start) latitude
'Start_Lon               lon           1    -91' ; ...  % start longitude
'StartLon                lon           1    -91' ; ...  % start longitude
'Start_Long              lon           1    -91' ; ...  % start longitude
'StartLong               lon           1    -91' ; ...  % start longitude
'Longitude               lon           1    -91' ; ...  % (start) longitude
'End_Lat                 lat           1    -92' ; ...  % end latitude
'EndLat                  lat           1    -92' ; ...  % end latitude
'End_Lon                 lon           1    -93' ; ...  % end longitude
'EndLon                  lon           1    -93' ; ...  % end longitude
'End_Long                lon           1    -93' ; ...  % end longitude
'EndLong                 lon           1    -93' ; ...  % end longitude
'Position                pos           2    -96' ; ...  % (start) lat  (start) lon
'Start_Pos               pos           2    -96' ; ...  % start lat  start lon
'StartPos                pos           2    -96' ; ...  % start lat  start lon
'End_Pos                 pos           2    -97' ; ...  % end lat  end lon
'EndPos                  pos           2    -97' ; ...  % end lat  end lon
'Columns                 list        -58   -100' ; ...  % Columns in data file
'Fields                  list        -58   -100' ; ...  % Columns in data file
'Columns cont.           list        -58   -101' ; ...  % Columns in data file
'Fields cont.            list        -58   -101' ; ...  % Columns in data file
'Units                   list        -58   -102' ; ...  % Units of Columns
'Units cont.             list        -58   -103' ; ...  % Units of Columns
'Sensors                 list        -58   -104' ; ...  % Sensors of Columns
'Sensors cont.           list        -58   -105' ; ...  % Sensors of Columns
'Recovery                list        -59   -201' ; ...  % Recovered Moorings
'Deployment              list        -59   -202' ; ...  % Deployed Moorings
'P                       data          0      1' ; ...  % pressure (dbar)
'D                       data          0      2' ; ...  % water depth (m)
'Z                       data          0      3' ; ...  % depth (m)
'MLD                     data          0      8' ; ...  % mixed layer depth (m)
'T                       data          0     20' ; ...  % temperature (deg C)
'WT                      data          0     19' ; ...  % wet bulb temperature (deg C)
'AT                      data          0     21' ; ...  % air temperature (deg C)
'AP                      data          0     22' ; ...  % air pressure (hPa)
'SST                     data          0     25' ; ...  % sea surface temperature (deg C)
'PT                      data          0     30' ; ...  % potential temp (deg C)
'TSGT                    data          0     35' ; ...  % TSG temperature (deg C)
'S                       data          0     41' ; ...  % salinity
'TSGS                    data          0     45' ; ...  % TSG salinity (PSU)
'C                       data          0     50' ; ...  % conductivity (mS/cm)
'TSGC                    data          0     55' ; ...  % TSG conductivity (mS/cm)
'O                       data          0     60' ; ...  % oxygen (ml/l)
'BO                      data          0     61' ; ...  % bottle oxygen (ml/l)
'DO                      data          0     67' ; ...  % oxygen difference (ml/l)
'ST                      data          0     70' ; ...  % sigma-t (kg/m^3)
'STH                     data          0     71' ; ...  % sigma-theta (kg/m^3)
'SV                      data          0     80' ; ...  % sound velocity (m/s)
'ID                      data          0     98' ; ...  % identifier
'NR                      data          0     99' ; ...  % number
'SN                      data          0    100' ; ...  % scan number (or serian number)
'SCN                     data          0    100' ; ...  % scan number
'BTL                     data          0    103' ; ...  % niskin bottle number
'PRF                     data          0    104' ; ...  % profile/cast number
'STN                     data          0    105' ; ...  % station
'OC                      data          0    110' ; ...  % oxygen current (micro A)
'OT                      data          0    111' ; ...  % oxygen temperature (deg C)
'DOC                     data          0    112' ; ...  % dOc/Dt (micro A/s)
'QS                      data          0    133' ; ...  % shortwave radiation (W/m^2)
'QL                      data          0    136' ; ...  % longwave radiation (W/m^2)
'QH                      data          0    137' ; ...  % latent heat flux (W/m^2)
'QB                      data          0    138' ; ...  % sensible heat flux (W/m^2)
'PR                      data          0    139' ; ...  % precipitation
'F11                     data          0    150' ; ...  % freon 11 (micro moles/kg)
'F12                     data          0    151' ; ...  % freon 12 (micro moles/kg)
'NO3                     data          0    182' ; ...  % nitrate (micro moles/l)
'NO2                     data          0    184' ; ...  % nitrite (micro moles/l)
'NN                      data          0    185' ; ...  % nitrate + nitrite (micro moles/l)
'PO4                     data          0    186' ; ...  % phosphate (micro moles/l)
'SI                      data          0    188' ; ...  % silicate (micro moles/l)
'CCL4                    data          0    190' ; ...  % ccl4, carbon tetrachloride (pico moles/kg)
'MTHCF                   data          0    191' ; ...  % methylchloroform (pico moles/kg)
'CHL_A                   data          0    200' ; ...  % chlorophyll_a (micro g/l)
'IC                      data          0    215' ; ...  % ice concentration (percent)
'CS                      data          0    300' ; ...  % current speed (cm/s)
'CD                      data          0    310' ; ...  % current direction (deg)
'U                       data          0    320' ; ...  % zonal current (cm/s)
'V                       data          0    321' ; ...  % meridional current (cm/s)
'W                       data          0    329' ; ...  % vertical velocity (cm/s)
'X                       data          0    370' ; ...  % x position (m)
'Y                       data          0    371' ; ...  % y position (m)
'L                       data          0    375' ; ...  % length, length scale (m)
'WS                      data          0    401' ; ...  % wind speed (m/s)
'WD                      data          0    410' ; ...  % wind direction (deg) (0-360)
'WU                      data          0    422' ; ...  % wind speed (m/s) (going E)
'WV                      data          0    423' ; ...  % wind speed (m/s) (going N)
'TX                      data          0    446' ; ...  % zonal wind stress (N/m^2)
'TY                      data          0    447' ; ...  % meridional wind stress (N/m^2)
'LAT                     data          0    500' ; ...  % latitude (degrees north)
'LON                     data          0    501' ; ...  % longitude (degrees east)
'YY                      data          0    601' ; ...  % year
'YYYY                    data          0    601' ; ...  % year
'MM                      data          0    602' ; ...  % month
'DD                      data          0    603' ; ...  % day
'HH                      data          0    604' ; ...  % hour (decimal)
'DINT                    data          0    610' ; ...  % Measurement Intervall [day]
'HINT                    data          0    611' ; ...  % Measurement Intervall [hour]
'MINT                    data          0    612' ; ...  % Measurement Intervall [min]
'SINT                    data          0    613' ; ...  % Measurement Intervall [sec]
'TIM                     data          0    626' ; ...  % time since start (s)
'SD                      data          0    700' ; ...  % standard deviation
'SK                      data          0    701' ; ...  % skewness
'PG                      data          0    721' ; ...  % ADCP percent good pings
'EV                      data          0    722' ; ...  % ADCP error velocity (cm/s)
'AMP                     data          0    723' ; ...  % ADCP amplitude
'SH                      data          0    913' ; ...  % Specific humidity (hPa)
'RH                      data          0    914' ; ...  % Relative humidity (%)
'AGC                     data          0   1202' ; ...  % ADCP Echo intensity (counts)
'HDG                     data          0   1215' ; ...  % ADCP heading (deg)
'PIT                     data          0   1216' ; ...  % ADCP pitch (deg)
'ROL                     data          0   1217' ; ...  % ADCP roll (deg)
'TLT                     data          0   1218' ; ...  % FSI ACM Tilt (deg)
'F113                    data          0   1703' ; ...  % freon 113 (pico moles/kg)
'TCO2                    data          0   1751' ; ...  % total carbon dioxide (micro moles)
'DIC                     data          0   1753' ; ...  % Dissolved Inorganic Carbon (micromoles/kg)
'ALK                     data          0   1756' ; ...  % alkalinity (micro moles/kg)
'T1                      data          0   9111' ; ...  % temp TC sensor 1 (deg C)
'T2                      data          0   9112' ; ...  % temp TC sensor 2 (deg C)
'T3                      data          0   9113' ; ...  % temp TC sensor 3 (deg C)
'T4                      data          0   9114' ; ...  % temp TC sensor 4 (deg C)
'T5                      data          0   9115' ; ...  % temp TC sensor 5 (deg C)
'T6                      data          0   9116' ; ...  % temp TC sensor 6 (deg C)
'T7                      data          0   9117' ; ...  % temp TC sensor 7 (deg C)
'T8                      data          0   9118' ; ...  % temp TC sensor 8 (deg C)
'T9                      data          0   9119' ; ...  % temp TC sensor 9 (deg C)
'T10                     data          0   9120' ; ...  % temp TC sensor 10 (deg C)
'T11                     data          0   9121' ; ...  % temp TC sensor 11 (deg C)
'T12                     data          0   9122' ; ...  % temp TC sensor 12 (deg C)
'PT1                     data          0   9131' ; ...  % pot temp TC sensor 1 (deg C)
'PT2                     data          0   9132' ; ...  % pot temp TC sensor 2 (deg C)
'PT3                     data          0   9133' ; ...  % pot temp TC sensor 3 (deg C)
'PT4                     data          0   9134' ; ...  % pot temp TC sensor 4 (deg C)
'PT5                     data          0   9135' ; ...  % pot temp TC sensor 5 (deg C)
'PT6                     data          0   9136' ; ...  % pot temp TC sensor 6 (deg C)
'PT7                     data          0   9137' ; ...  % pot temp TC sensor 7 (deg C)
'PT8                     data          0   9138' ; ...  % pot temp TC sensor 8 (deg C)
'PT9                     data          0   9139' ; ...  % pot temp TC sensor 9 (deg C)
'PT10                    data          0   9140' ; ...  % pot temp TC sensor 10 (deg C)
'PT11                    data          0   9141' ; ...  % pot temp TC sensor 11 (deg C)
'PT12                    data          0   9142' ; ...  % pot temp TC sensor 12 (deg C)
'DPT                     data          0   9143' ; ...  % Depth calc. with Moordyn (m) 
'                                              '       ];
