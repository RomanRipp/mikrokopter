function d = read_coriolis(file);

% READ_CORIOLIS  Reads Data from Coriolis Ascii- or NetCDF-DataFile
%
% required Variables in Coriolis-NetCDF-File:
%
%   'PLATFORM_NUMBER' 'INST_REFERENCE'
%   'JULD' 'LATITUDE' 'LONGITUDE'
%   'DEPH'    'PRES'    'TEMP'    'PSAL'    'CHLT'
%   'DEPH_QC' 'PRES_QC' 'TEMP_QC' 'PSAL_QC' 'CHLT_QC'
%
% Code for Coriolis DataFiles:
% --------------------------------------------------------------
% Code   Meaning                           Comments
% --------------------------------------------------------------
% BA     BATHY
% BO     Bottle
% CD     CTD down cast
% CU     CTD up cast
% CT     CTD up or down cast              (Glider and CTD-Data)
% DB     DRIFTER
% MB     MBT
% PF     Profiling float
% TE     TESAC                            (evt. Mooring-Data)
% TK     TRACKOB
% TR     THERMISTOR CHAIN (delayed mode)
% XB     XBT
%
% required HeaderLines in Coriolis-ASCII-File:
%
%  WMO PLATFORM CODE : 64556
%  PLATFORM NAME : MERSEA GLIDER 1 SPRAY 2006
%  *DATE=12032006 TIME=1349 LAT=N47 38.35 LON=W008 28.18 DEPTH=  QC=1119
%  *NB PARAMETERS=04 RECORD LINES=00284
%  *Station number : 00009
%  *SURFACE SAMPLES=
%  *PRES   TEMP   PSAL   CHLT
%
% Returns Data-Structure:
%
%        day: [     1 x NPROF  double ]
%        lat: [     1 x NPROF  double ]
%        lon: [     1 x NPROF  double ]
%         qc: [ NPROF x 1      double ] ASCII only 
%         nr: [ NPROF x 1      double ] ASCII only, StationNumber
%         nn: [ NPROF x 1      double ] ASCII only, RecordLines
%       deph: [ NDEPH x NPROF  double ]
%       pres: [ NDEPH x NPROF  double ]
%       temp: [ NDEPH x NPROF  double ]
%       psal: [ NDEPH x NPROF  double ]
%       chlt: [ NDEPH x NPROF  double ]
%    deph_qc: [ NDEPH x NPROF  double ]
%    pres_qc: [ NDEPH x NPROF  double ]
%    temp_qc: [ NDEPH x NPROF  double ]
%    psal_qc: [ NDEPH x NPROF  double ]
%    chlt_qc: [ NDEPH x NPROF  double ]
%       name: { NTYPE x     1  cell   } "PLATFORM_NUMBER & INST_REFERENCE"
%         id: [     1 x NPROF  double ] [ 1 .. NTYPE ]
%
% Quality Flags: "0" No QC, "1" Correct, "2" Inconsistent, 
%                "3" Doubtful, "4" Bad, "5" Changed, 
%                "9" Missing, " " FillValue
%
% Negative Quality-Flag: Parameter set to NaN by READ_CO, outside valid Range
% 
% Negative Quality-Flag and nonzero imaginary Part for psal_qc (Salinity): 
%   Parameter set to NaN if Temperature out of Range.
% 
% required M-Files:
% 
% ASSIGN_CDF
% P2Z80, Z2P80
% STRWCMP, UNIQUED
% STRUCT/RNFIELD, STRUCT/INSFIELD
%
% BGREP, MGREP, READ_ASC, STR2VEC, STRWCMP
%

var0  = { 'PLATFORM_NUMBER' 'INST_REFERENCE' };
var1  = { 'JULD' 'LATITUDE' 'LONGITUDE' };   
var11 = { 'day' 'lat' 'lon' };
var12 = { 'qc'  'nr'  'nn'  };
var2  = { 'DEPH'    'PRES'    'TEMP'    'PSAL'    'CHLT' };
var3  = { 'DEPH_QC' 'PRES_QC' 'TEMP_QC' 'PSAL_QC' 'CHLT_QC' };
var32 = { 'name' 'id' 'info' };

txt_init = { 'DATE=' '*NB PARAMETERS=' '*Station number :' };

f = which(file);
if ~isempty(f)
    file = f;
end

if ~( exist(file,'file') == 2 )
    error('File not found.');
end

% ASSIGN Data from NetCDF-File
try
   [m,d] = assign_cdf(file,cat(2,var0,var1,var2,var3),NaN);
catch
    m = sprintf('Error using ASSIGN_CDF.\n%s',lasterr);
end

ok = ~isempty(m); 
if ok
   % ASSIGN_CDF failed, check for ASCII-File
   for t = txt_init
       ok = ~isempty(bgrep(t{1},file,'-d','-B2048','-M10'));
       if ~ok
           break
       end
   end
end

if ok
   m = '';
   try
      d = read_co_txt(file);
   catch
      m = sprintf('Error read ASCII-File.\n%s',lasterr);
   end
end

if ~isempty(m)
    error(m)
end

[p,n,e] = fileparts(file);

%*******************************************************************
if ok
%*******************************************************************
% Transform Ascii-Data

d = insfield(d,{'deph'},'pres');

d = insfield(d,lower(var2{end}),lower(var3));

def = NaN * d.temp;

for ii = 1 : size(var3,2)
    d = setfield(d,lower(var3{ii}),def);
end

d.id = NaN * ones(1,size(def,2));

%*******************************************************************
else
%*******************************************************************
% Transform NetCDF-Data

%-------------------------------------------------------------------
% Check for Unique PlatForm, Remove/Rename Fields

[nr,iy,nn,ix] = uniqued(d.PLATFORM_NUMBER,2+i);

d.name = cat( 1 , d.PLATFORM_NUMBER(:,iy) , ...
                  char(32*ones(size(iy))) , ...
                  d.INST_REFERENCE(:,iy)        );

d.name = cellstr(d.name');

d.id   = ix;

d =  rmfield(d,var0);
d =  rnfield(d,var1,var11);
d = insfield(d,var11{end},var12);
d =  rnfield(d,var2,lower(var2));
d =  rnfield(d,var3,lower(var3));

def = NaN * d.day;

for ii = 1 : size(var12,2)
    d = setfield(d,var12{ii},def);
end

d.day = d.day + datenum(1950,01,01);

%--------------------------------------------------------------
% Convert Flags for Quality: char --> double


f   = fieldnames(d);
iqc = find(strwcmp(f,'*_qc'));
for ii = iqc(:)'
    v = double(getfield(d,f{ii}));  % " 0..9" --> [ 32 48..57 ]
    v(find(v==32)) = NaN;           % FillValue
    v = v - 48;                     % "0..9" --> [48..57] --> [0..9]
    d = setfield(d,f{ii},v);
    if ~all( isnan(v(:)) | ( v(:) == 1 ) )
       nn = hist(v(:),0:9);
       fprintf('%s %s: ',n,f{ii});
       fprintf('%3.0f ',nn);
       fprintf('\n');
    end
end

%*******************************************************************
end
%*******************************************************************

%-------------------------------------------------------------------
% Check for DEPH / PRES

if isempty(d.deph) & ~isempty(d.pres)
   d.deph = p2z80(d.pres,ones(size(d.pres,1),1)*d.lat);
   d.deph_qc = NaN * d.deph;
end

if isempty(d.pres) & ~isempty(d.deph)
   d.pres = z2p80(d.deph,ones(size(d.deph,1),1)*d.lat);
   d.pres_qc = NaN * d.pres;
end

%-------------------------------------------------------------------
% Check TEMP and PSAL with Valid Range

   tlm = [  2   30 ];
   slm = [ 30   40 ]; 

   t2p = [];

if     strwcmp(n,'*_cis*')

   tlm = [  2   18 ];
   slm = [ 34.5 38 ];

   t2p = [   0 12
           100  8
           500  6
          1000  5 ];

elseif strwcmp(n,'*_pap*')

   tlm = [  2   18 ];
   slm = [ 34.5 38 ];

   t2p = [   0 15
          1000 12
          2000  5 ];

elseif strwcmp(n,'*_pap*')

   tlm = [ 10   30 ];
   slm = [ 36   40 ];

   t2p = [   0 30
           500 15
          1000 14
          2000 14 ];
end

ii = find( ( d.temp < tlm(1) ) | ( tlm(2) < d.temp ) );
if ~isempty(ii)
    d.temp(ii) = NaN;
    d.temp_qc(ii) = -d.temp_qc(ii);
    if ~isempty(d.psal)
        d.psal(ii) = NaN;
        d.psal_qc(ii) = -d.psal_qc(ii) + i;
    end
end

if ~isempty(d.psal)
    ii = find( ( d.psal < slm(1) ) | ( slm(2) < d.psal ) );
    if ~isempty(ii)
        d.psal(ii) = NaN;
        d.psal_qc(ii) = -d.psal_qc(ii);
    end
end

%-------------------------------------------------------------------
% Check CHLT in Depth

if ~isempty(d.chlt)
    ii = find( ( d.pres > 100 ) & ( d.chlt > 1 ) );
    if ~isempty(ii)
        d.chlt(ii) = NaN;
        d.chlt_qc(ii) = -d.temp_qc(ii);
    end
end

%%% round(max(d.PRES(:)))


%-------------------------------------------------------------------
dp = cat(1,NaN*d.pres(1,:),diff(d.pres,1,1));
ok = ( ( dp > 100 ) & ~isnan(d.temp) );
if any(ok(:))
   figure,plot(d.temp,d.pres,'.');
   title([n ' ' num2str(sum(ok(:))) ])
end

%-------------------------------------------------------------------
% InfoText

d.info = cell(size(d.name,1),1);

for ii = 1 : size(d.name,1)
    jj = find( d.id == ii );
    d.info{ii} = sprintf(' %s, %s,%7.2f,%7.2f,%6.2f,%6.2f, %s', ...
                           datestr(min(d.day(jj)),29) , ...
                           datestr(max(d.day(jj)),29) , ...
                           min(d.lon(jj)) , max(d.lon(jj)) , ...
                           min(d.lat(jj)) , max(d.lat(jj)) , ...
                          d.name{ii} );
end


%***********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function dat = read_co_txt(file)

% Reads Coriolis Ascii-DataFile for GliderData
%
%
% Required HeaderLines:
%
%  WMO PLATFORM CODE : 64556
%  PLATFORM NAME : MERSEA GLIDER 1 SPRAY 2006
%  *DATE=12032006 TIME=1349 LAT=N47 38.35 LON=W008 28.18 DEPTH=  QC=1119
%  *NB PARAMETERS=04 RECORD LINES=00284
%  *Station number : 00009
%
% Returns Data-Structure:
%
%     day: [ NPROF x 1      double ]
%     lat: [ NPROF x 1      double ]
%     lon: [ NPROF x 1      double ]
%      qc: [ NPROF x 1      double ]
%      nr: [ NPROF x 1      double ] StationNumber
%      nn: [ NPROF x 1      double ] RecordLines
%    pres: [ NDEPH x NPROF  double ]
%    temp: [ NDEPH x NPROF  double ]
%    psal: [ NDEPH x NPROF  double ]
%    chlt: [ NDEPH x NPROF  double ]
%    name: [     1 x 1      cell   ]
%
% required M-Files:
%
% MGREP, READ_ASC, STR2VEC, STRWCMP
%

f = which(file);
if ~isempty(f)
    file = f;
end

%**********************************************************************
% Inquire HeaderData:
%
%*DATE=12032006 TIME=1349 LAT=N47 38.35 LON=W008 28.18 DEPTH=  QC=1119
%*NB PARAMETERS=04 RECORD LINES=00284
%*Station number : 00009

%-------------------------------------------------------
% GREP for specific HeaderLines in File

nc = mgrep('WMO PLATFORM CODE : ' ,file,'-M10','-d');
nm = mgrep('PLATFORM NAME : '     ,file,'-M10','-d');


hd = mgrep('*DATE='           ,file,'-M10','-d');
nn = mgrep('*NB PARAMETERS='  ,file,'-M10','-d');
nr = mgrep('*Station number :',file,'-M10','-d');


is_south = strwcmp(hd(:,end),'*LAT=S*'); % True for Southern Latitude
is_west  = strwcmp(hd(:,end),'*LON=W*'); % True for Western Longitude

%-------------------------------------------------------
% Check for Codes and Name

if ~isempty(nc) & ~isempty(nm)
    nc = strrep(nc(:,end),'WMO PLATFORM CODE : ','');
    nm = strrep(nm(:,end),'PLATFORM NAME : ',' ');
    nc = char(cellstr(char(nc)));
    nm = char(cellstr(char(nm)));
    nm = cellstr([nc nm]);
else
    nm = {};
end

%-------------------------------------------------------
% Convert Strings to numeric

[m,hd] = str2vec(sprintf('%s\n',hd{:,end}));
   hd  = reshape(hd,7,length(hd)/7)';  % [ DDMMYYYY HHMM LatD LatM LonD LonM QC ]

[m,nn] = str2vec(sprintf('%s\n',nn{:,end}));
   nn  = reshape(nn,2,length(nn)/2)';  % [ NParam  RecLines ]

[m,nr] = str2vec(sprintf('%s\n',nr{:,end}));
   nr  = reshape(nr,1,length(nr)/1)';  % [ STatNr ]

%-------------------------------------------------------
% Get Position from HD

lat = ( 1 - 2*is_south ) .* ( hd(:,3) + hd(:,4)/60 ); 
lon = ( 1 - 2*is_west  ) .* ( hd(:,5) + hd(:,6)/60 );

%-------------------------------------------------------
% Get Date from HD

dd = floor(hd(:,1)/1e6);
mm = floor( ( hd(:,1) - dd*1e6 ) / 1e4 );
yy = hd(:,1) - 1e6*dd - 1e4 * mm;

hr = floor(hd(:,2)/1e2);
mn = hd(:,2) - hr*1e2;

%-------------------------------------------------------
% Create DataStructure with HeaderData 
%
% *PRES   TEMP   PSAL   CHLT

nd = size(hd,1);  % Number of Profiles

def = NaN * zeros(max(nn(:,2)),nd);  % Default Data-Matrice

dat = struct( 'day' , { datenum(yy,mm,dd,hr,mn,0*mn)' } , ...
              'lat' , { lat' }     , ...
              'lon' , { lon' }     , ...
              'qc'  , { hd(:,7)' } , ...
              'nr'  , { nr'      } , ...
              'nn'  , { nn(:,2)' } , ...
             'pres' , { def }      , ...
             'temp' , { def }      , ...
             'psal' , { def }      , ...
             'chlt' , { def }      , ...
             'name' , { nm }            );
              
%**********************************************************************
% Read numeric Data from File

[m,h,d] = read_asc(file);

ii = find( d(:,1) <= -9 );  % DummyLines at End of Profile

d(ii,:) = [];

if ~( size(d,1) == sum(nn(:,2)) )
    error(sprintf('Read %.0f Lines, expected %.0f',size(d,1),sum(nn(:,2))));
end

% Sort Variables into Structure

i0 = cumsum( cat(1,0,nn(:,2)) , 1 ); % StartIndex - 1

for ii = 1 : nd
    jj = ( 1 : nn(ii,2) );
    dat.pres(jj,ii) = d(i0(ii)+jj,1);
    dat.temp(jj,ii) = d(i0(ii)+jj,2);
    dat.psal(jj,ii) = d(i0(ii)+jj,3);
    dat.chlt(jj,ii) = d(i0(ii)+jj,4);
end
