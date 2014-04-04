function Msg = ctd2cdf(ncfile,lon,lat,dpt,time,pres,temp,psal);

% CTD2CDF  Writes CTD-Data to NetCDF-File, using NetCDF-ModelFile "ctd.mnc"
%
%  Msg = CTD2CDF( NCFileName , Lon , Lat , Depth , Time , P , T , S );
%
%  Lon   Longitude    [  1 by NProfile ]                          LONGITUDE
%  Lat   Latitude     [  1 by NProfile ]                          LATITUDE 
%  Depth BottomDepth  [  1 by NProfile ]                          BOTTOM_DEPTH
%  Time  Time         [ NProfile by 6  ]  [ DD MM YY hh mm ss ]   DATE
%  P     Pressure     [ NZ by NProfile ]                          PRES
%  T     Temperature  [ NZ by NProfile ]                          TEMP
%  S     Salinity     [ NZ by NProfile ]                          PSAL
%
% Msg returns ErrorMessages, empty if all was succesfull.
%
% see also:  LOOK_CDF,  WRT_CDF
%


DefaultFile = 'ctd.mnc';

%  'DD/MM/YYYY hh:mm:ss'

DateForm   = [ '%2.2d/%2.2d/%4.0f %2.0f:%2.2d:%2.2d' ];


Msg = '';
nl  = char(10);

if nargin < 8
   Msg = 'Not enough InputArguments.'
   return
end

lat = lat(:)';
lon = lon(:)';
dpt = dpt(:)';

Np = size(lat,2);
Nz = size(pres,1);
   
if ( size(time,2) ~= 6 )
   if ( Np == size(time,2) )  &  ( size(time,1) == 6 )
     time = time';
   end
end

ok = ( ( Np == size(lon,2) )  &  ...
       ( Np == size(dpt,2) )  &  ...
       ( Np == size(time,1) )  &  ...
       ( Np == size(pres,2) )  &  ...
       ( Np == size(temp,2) )  &  ...
       ( Np == size(psal,2) )  &  ...
       ( Nz == size(pres,1) )  &  ...
       ( Nz == size(temp,1) )  &  ...
       ( Nz == size(psal,1) )         );

if ~ok
  Msg = 'Matrix Dimensions must be agree.';
end

if size(time,2) ~= 6
  Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
          'Time must have 6 Columns: [ DD MM YY hh mm ss ].' ];
end


if ~isempty(Msg)
   return
end


% create DateString

 Date = sprintf([ DateForm  char(10) ] , time' );

 % NewLine --> "';'"
 Date = strrep(Date,char(10),char([39  59  39]));
 Date = eval( char([  '{''' Date  '''}' ]) );

 Date = char(Date(1:Np));


dims = { 'N_DATE_TIME'   size(Date,2)  
         'mN_PROF'       Np
         'mN_ZLEV'       Nz            };

vars = { 'DATE'          Date'
         'LATITUDE'      lat  
         'LONGITUDE'     lon
         'BOTTOM_DEPTH'  dpt
         'PRES'          pres
         'TEMP'          temp
         'PSAL'          psal  };


Msg = wrt_cdf(ncfile,DefaultFile,dims,vars);

