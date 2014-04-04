if strcmp(getenv('HOSTNAME'),'geo')

   dbroot /export/geo/data/moorings/data/ begin
   dbroot /export/geo/data/cruises/data/  begin

else %%% may be on CRONOS

   dbroot /po2data/moorings/data/ begin
   dbroot /po2data/cruises/data/  begin

end


% Inquire InfoFiles (Mooring- and CruiseInfo)
dbfile move:*::info.dat


% Inquire CTD-DataFiles like "cruise/ctd/asc/*.ctd"
dbfile move:*:ctd:asc/*.ctd


% Inquire CTD-DataFiles like "cruise/ctd/cruise_###.ctd
dbfile move:*:ctd


% Load CTD-Data from Atalante 2002-Cruise

d = dbload('move:ata2002:ctd:asc/*.ctd','Latitude:Longitude:P:T:S');


% Load CTD-Data from Atalante 2005-Cruise, Profiles 8 .. 10 

d = dbload('move:ata2005:ctd:[8 9 10]','Latitude:Longitude:P:T:S');


% Inquire MicroCat-Data

dbfile move:*:mc


% Load MicroCatData from Mooring V404_3, the first 5 Instruments

d = dbload('move:v404_3:mc:(1:5)','YY:MM:DD:HH:P:T:C')

t = datenum(d.YY,d.MM,d.DD)+d.HH/24;  % DayVector

if strcmp(getenv('HOSTNAME'),'geo')
   % Path to TIMEAXIS
   addpath  /export/geo/data/share/software/matlab/tools/toolbox/graphics/axis
end

figure , plot(t,d.T), timeaxis(gca);

