function [adcp,cfg,ens]=rd_adcp(name,varargin);

% RD_ADCP  Read (raw binary) RDI ADCP files, 
%
%  ADCP=RD_ADCP(NAME) reads the raw binary RDI BB/Workhorse ADCP file NAME and
%  puts all the relevant configuration and measured data into a data structure 
%  ADCP (which is self-explanatory). This program is designed for handling data
%  recorded by moored instruments (primarily Workhorse-type but can also read
%  Broadband) and then downloaded post-deployment. For vessel-mount data I
%  usually make p-files (which integrate nav info and do coordinate transformations)
%  and then use RDPADCP. 
%
%  This current version does have some handling of VMDAS and WINRIVER output
%  files, but it is still 'beta'.
%
%  [ADCP,CFG]=RD_ADCP(...) returns configuration data in a
%  separate data structure.
%
%  Various options can be specified on input:
%  [..]=RD_ADCP(NAME,NUMAV) averages NUMAV ensembles together in the result.
%  [..]=RD_ADCP(NAME,NUMAV,NENS) reads only NENS ensembles (-1 for all).
%  [..]=RD_ADCP(NAME,NUMAV,[NFIRST NEND]) reads only the specified range
%   of ensembles. This is useful if you want to get rid of bad data before/after
%   the deployment period.
%
%  Note - sometimes the ends of files are filled with garbage. In this case you may
%         have to rerun things explicitly specifying how many records to read (or the
%         last record to read). I don't handle bad data very well.
%
%       - I don't read in absolutely every parameter stored in the binaries;
%         just the ones that are 'most' useful.
%
%  String parameter/option pairs can be added after these initial parameters:
%
%  'baseyear'    : Base century for BB/v8WH firmware (default to 2000).
%
%  'despike'    : [ 'no' | 'yes' | 3-element vector ]
%                 Controls ensemble averaging. With 'no' a simple mean is used 
%                 (default). With 'yes' a mean is applied to all values that fall 
%                 within a window around the median (giving some outlier rejection). 
%                 This is useful for noisy data. Window sizes are [.3 .3 .3] m/s 
%                 for [ horiz_vel vert_vel error_vel ] values. If you want to 
%                 change these values, set 'despike' to the 3-element vector.
%
% R. Pawlowicz (rich@ocgy.ubc.ca) - 17/09/99
%
% Revision Ch. Begler, 2007-2008
%

% R. Pawlowicz - 17/Oct/99 
%          5/july/00 - handled byte offsets (and mysterious 'extra" bytes) slightly better, Y2K
%          5/Oct/00 - bug fix - size of ens stayed 2 when NUMAV==1 due to initialization,
%                     hopefully this is now fixed.
%          10/Mar/02 - #bytes per record changes mysteriously,
%                      tried a more robust workaround. Guess that we have an extra
%                      2 bytes if the record length is even?
%          28/Mar/02 - added more firmware-dependent changes to format; hopefully this
%                      works for everything now (put previous changes on firmer footing?)
%          30/Mar/02 - made cfg output more intuitive by decoding things.
%                    - An early version of WAVESMON and PARSE which split out this
%                      data from a wave recorder inserted an extra two bytes per record.
%                      I have removed the code to handle this but if you need it see line 509
%         29/Nov/02  - A change in the bottom-track block for version 4.05 (very old!).
%         29/Jan/03  - Status block in v4.25 150khzBB two bytes short?
%         14/Oct/03  - Added code to at least 'ignore' WinRiver GPS blocks.
%         11/Nov/03  - VMDAS navigation block, added hooks to output
%                      navigation data.

init = '7F7F';    % Init-Bytes for Ensamble

num_av  =  1;     % Block filtering and decimation parameter (# ensembles to block together).
nens    = -1;     % Read all ensembles.
century = 2000;   % ADCP clock does not have century prior to firmware 16.05.
vels    = 'no';   % Default to simple averaging

lv=length(varargin);
if lv>=1 & ~isstr(varargin{1}),
  num_av=varargin{1}; % Block filtering and decimation parameter (# ensembles to block together).
  varargin(1)=[];
  lv=lv-1;
  if lv>=1 & ~isstr(varargin{1}),
    nens=varargin{1};
    varargin(1)=[];
    lv=lv-1;
  end;
end;

% Read optional args
while length(varargin)>0,
 switch varargin{1}(1:3),
	 case 'bas',
	   century = varargin{2};
	 case 'des',
	   if isstr(varargin{2}),
	    if strcmp(varargin{2},'no'), vels='no';
	    else vels=[.3 .3 .3]; end;
	   else
	    vels=varargin{2}; 
	   end;   
	 otherwise,
	   error(['Unknown command line option  ->' varargin{1}]);
   end;
   varargin([1 2])=[];
end;   	          	

%--------------------------------------------------------
% Info about M-File

mf = which(mfilename);

 d = dir(mf);

mfile = struct( 'path' , { fileparts(mf) } , ...
                'name' , { d.name } , ...
                'date' , { d.date } , ...
               'bytes' , { d.bytes }   );

%--------------------------------------------------------
% Check file information first

file = which(name);
if isempty(file)
   file = name;
end

if ~( exist(file,'file') == 2 )
    fprintf('ERROR******* Can''t find file %s\n',name);
    return;
end

d = dir(file);

bytes = d.bytes;

if ( bytes == 0 )
   fprintf('ERROR******* Empty file %s\n',file);
   return;
end

%--------------------------------------------------------
% Make Info

dfile = struct( 'path' , { fileparts(file) } , ...
                'name' , { d.name } , ...
                'date' , { d.date } , ...
               'bytes' , { d.bytes }      );


usr = getenv('USER'); if isempty(usr), usr = getenv('USERNAME'); end
hst = getenv('HOST'); if isempty(hst), hst = getenv('HOSTNAME'); end

typ = 'Unmodified';
if num_av > 1
   typ = sprintf('%.0fens avg.',num_av);
end

typ = sprintf('%s Data from binary RDI-DataFile',typ);

info = struct( 'author'  , { sprintf('%s@%s',usr,hst) }    , ...
               'date'    , { datestr(clock,0) }            , ...
               'type'    , { typ }  , ...
               'total'   , { NaN } , ...
               'read'    , { NaN*ones(1,2) } , ...
               'avg'     , { num_av } , ...
               'version' , { sprintf('MATLAB %s',version) } , ...
                 'mfile' , { mfile } , ...
                 'dfile' , { dfile } );


%--------------------------------------------------------
% Open File, Read first Ensamble

fprintf('\nOpening file %s\n\n',file);
fd = fopen(file,'r','ieee-le');

% Read first ensemble to initialize parameters

[ens,hdr,cfg,ok]=rd_buffer(fd,init,-2,0); % Initialize and read first two records

fseek(fd,0,'bof');              % Rewind
 
if (cfg.prog_ver<16.05 & cfg.prog_ver>5.999) | cfg.prog_ver<5.55,
  fprintf('**************Assuming that the century begins year %d *********** \n\n',century);
else
  century=0;  % century included in clock.  
end;

dats=datenum(century+ens.rtc(1,:),ens.rtc(2,:),ens.rtc(3,:),ens.rtc(4,:),ens.rtc(5,:),ens.rtc(6,:)+ens.rtc(7,:)/100);
t_int=diff(dats);
fprintf('Record begins at %s\n',datestr(dats(1),0));
fprintf('Ping interval appears to be %s\n',datestr(t_int,13));


% Estimate number of records (since I don't feel like handling EOFs correctly,
% we just don't read that far!)


% Now, this is a puzzle - it appears that this is not necessary in
% a firmware v16.12 sent to me, and I can't find any example for
% which it *is* necessary so I'm not sure why its there. It could be
% a leftoever from dealing with the bad WAVESMON/PARSE problem (now
% fixed) that inserted extra bytes.
% ...So its out for now.
%if cfg.prog_ver>=16.05, extrabytes=2; else extrabytes=0; end; % Extra bytes
extrabytes=0;

n0=fix(bytes/(hdr.nbyte+2+extrabytes));

info.total = n0;

fprintf('\nEstimating %d ensembles in this file.\n',n0);  

if length(nens)==1,
  if nens==-1
     nens = n0;
    fprintf('Reading all Ensembles');  
  else
     fprintf('Reading %d Ensembles',nens); 
  end; 
     nens = [ 1  nens ];
     nens = min(max(nens,1),n0);
else
  nens = min(max(nens,1),n0);
  fprintf('\nReading Ensembles %d-%d',nens); 
  fseek(fd,(hdr.nbyte+2+extrabytes)*(nens(1)-1),'bof');
end;

fprintf(', Reducing by a factor of %d\n',num_av);  


info.read = nens;


nens = diff(nens) + 1;


if num_av>1,
  if isstr(vels),
     fprintf('\n Simple mean used for ensemble averaging\n');
  else
     fprintf('\n Averaging after outlier rejection with parameters [%f %f %f]\n',vels);
  end;
end;
   
% Number of records after averaging.

n=fix(nens/num_av);

% Structure to hold all ADCP data 
% Note that I am not storing all the data contained in the raw binary file, merely
% things I think are useful.

switch cfg.sourceprog,
  case 'WINRIVER',
    adcp=struct('info',info,'name','adcp','config',cfg,'ok',zeros(1,n), ...
                'xtime',zeros(2,n),'mtime',zeros(1,n),'number',zeros(1,n),...
                'orientation',zeros(1,n),'heading',zeros(1,n),'pitch',zeros(1,n),'roll',zeros(1,n), ...
                'heading_std',zeros(1,n),'pitch_std',zeros(1,n),'roll_std',zeros(1,n),'depth',zeros(1,n),...
        	'temperature',zeros(1,n),'salinity',zeros(1,n),'ssp',zeros(1,n),...
        	'pressure',zeros(1,n),'pressure_std',zeros(1,n),'adc',zeros(8,n),...
        	'east_vel',zeros(cfg.n_cells,n),'north_vel',zeros(cfg.n_cells,n),'vert_vel',zeros(cfg.n_cells,n),...
        	'error_vel',zeros(cfg.n_cells,n),'corr',zeros(cfg.n_cells,4,n),...
        	'status',zeros(cfg.n_cells,4,n),'intens',zeros(cfg.n_cells,4,n),...
	        'bt_range',zeros(4,n),'bt_vel',zeros(4,n),...
            'nav_longitude',zeros(1,n),'nav_latitude',zeros(1,n));
  case 'VMDAS',
    adcp=struct('info',info,'name','adcp','config',cfg,'ok',zeros(1,n), ...
                'xtime',zeros(2,n),'mtime',zeros(1,n),'number',zeros(1,n),...
                'orientation',zeros(1,n),'heading',zeros(1,n),'pitch',zeros(1,n),'roll',zeros(1,n), ...
                'heading_std',zeros(1,n),'pitch_std',zeros(1,n),'roll_std',zeros(1,n),'depth',zeros(1,n),...
        	'temperature',zeros(1,n),'salinity',zeros(1,n),'ssp',zeros(1,n),...
        	'pressure',zeros(1,n),'pressure_std',zeros(1,n),'adc',zeros(8,n),......
        	'east_vel',zeros(cfg.n_cells,n),'north_vel',zeros(cfg.n_cells,n),'vert_vel',zeros(cfg.n_cells,n),...
        	'error_vel',zeros(cfg.n_cells,n),'corr',zeros(cfg.n_cells,4,n),...
        	'status',zeros(cfg.n_cells,4,n),'intens',zeros(cfg.n_cells,4,n),...
	        'bt_range',zeros(4,n),'bt_vel',zeros(4,n),...
	        'nav_smtime',zeros(1,n),'nav_emtime',zeros(1,n),...
	        'nav_slongitude',zeros(1,n),'nav_elongitude',zeros(1,n),...
	        'nav_slatitude',zeros(1,n),'nav_elatitude',zeros(1,n));
  otherwise 
    adcp=struct('info',info,'name','adcp','config',cfg,'ok',zeros(1,n), ...
                'xtime',zeros(2,n),'mtime',zeros(1,n),'number',zeros(1,n),...
                'orientation',zeros(1,n),'heading',zeros(1,n),'pitch',zeros(1,n),'roll',zeros(1,n), ...
                'heading_std',zeros(1,n),'pitch_std',zeros(1,n),'roll_std',zeros(1,n),'depth',zeros(1,n),...
        	'temperature',zeros(1,n),'salinity',zeros(1,n),'ssp',zeros(1,n),...
        	'pressure',zeros(1,n),'pressure_std',zeros(1,n),'adc',zeros(8,n),......
        	'east_vel',zeros(cfg.n_cells,n),'north_vel',zeros(cfg.n_cells,n),'vert_vel',zeros(cfg.n_cells,n),...
        	'error_vel',zeros(cfg.n_cells,n),'corr',zeros(cfg.n_cells,4,n),...
        	'percent',zeros(cfg.n_cells,4,n),'intens',zeros(cfg.n_cells,4,n),...
		    'bt_range',zeros(4,n),'bt_vel',zeros(4,n));
end;


% Calibration factors for backscatter data

clear global ens

len = (hdr.nbyte+2+extrabytes);

% Use an Number which has not any DimLength of a Variable
n_expand = 13;   

% Loop for all records
for k=1:n,

  % Gives display so you know something is going on...
    
  if rem(k,50)==0,  fprintf('%d\n',k*num_av);end;

  pos = ftell(fd); % 

  if mod(pos,len) == 0
     fprintf('.');
  else
     pos = len*round(pos/len);
     fseek(fd,pos,'bof');
     fprintf('*');
  end

  % Read an ensemble
  
try
  [ens,hd,cf,ok]=rd_buffer(fd,init,num_av,k);
catch
  fprintf('%s\n',lasterr);
  ens = [];
  hd  = [];
  cf  = [];
  ok  = [];
end

  if ~isstruct(ens), % If aborting...
    fprintf('Only %d records found..suggest re-running RD_ADCP using this parameter\n',(k-1)*num_av);
    fprintf('(If this message preceded by a POSSIBLE PROGRAM PROBLEM message, re-run using %d)\n',(k-1)*num_av-1);
    break;
  end;

  %% Check for valid Year !!!
  yy = century + ens.rtc(1,:);
  kk = find( ( 1980 <= yy ) & ( yy <= 2100 ) );
  if isempty(kk)
     yy = ens.rtc(1,:);
     kk = find( ( 1980 <= yy ) & ( yy <= 2100 ) );
     if isempty(kk)
        kk = find( yy <= 100 );
        if isempty(kk)
           kk = (1:size(ens.rtc,2));
        end
     end
  end

  dats = datenum( yy , ens.rtc(2,:) , ens.rtc(3,:) , ...
                  ens.rtc(4,:),ens.rtc(5,:),ens.rtc(6,:)+ens.rtc(7,:)/100);

  adcp.ok(k) = ok;
  
  adcp.xtime(:,k)     = [ min(dats) ; max(dats) ];  
  adcp.mtime(k)       =  median(dats(kk));  
  adcp.number(k)      =ens.number(1);
  adcp.orientation(k) =cf.orientation;
  adcp.heading(k)     =mean(ens.heading);
  adcp.pitch(k)       =mean(ens.pitch);
  adcp.roll(k)        =mean(ens.roll);
  adcp.heading_std(k) =mean(ens.heading_std);
  adcp.pitch_std(k)   =mean(ens.pitch_std);
  adcp.roll_std(k)    =mean(ens.roll_std);
  adcp.depth(k)       =mean(ens.depth);
  adcp.temperature(k) =mean(ens.temperature);
  adcp.salinity(k)    =mean(ens.salinity);
  adcp.ssp(k)         =mean(ens.ssp);
  adcp.pressure(k)    =mean(ens.pressure);
  adcp.pressure_std(k)=mean(ens.pressure_std);
  adcp.adc(:,k)       = mean(ens.adc,2);

  if isstr(vels),
    adcp.east_vel(:,k)    =nmean(ens.east_vel ,2);
    adcp.north_vel(:,k)   =nmean(ens.north_vel,2);
    adcp.vert_vel(:,k)    =nmean(ens.vert_vel ,2);
    adcp.error_vel(:,k)   =nmean(ens.error_vel,2);
  else
   adcp.east_vel(:,k)    =nmedian(ens.east_vel  ,vels(1),2);
   adcp.north_vel(:,k)   =nmedian(ens.north_vel,vels(1),2);
   adcp.vert_vel(:,k)    =nmedian(ens.vert_vel  ,vels(2),2);
   adcp.error_vel(:,k)   =nmedian(ens.error_vel,vels(3),2);
  end;
  
  adcp.corr(:,:,k)      =nmean(ens.corr,3);        % added correlation RKD 9/00
  
  adcp.intens(:,:,k)   =nmean(ens.intens,3);

  adcp.percent(:,:,k)  =nmean(ens.percent,3);   
  
  adcp.bt_range(:,k)   =nmean(ens.bt_range,2);
  adcp.bt_vel(:,k)     =nmean(ens.bt_vel,2);

  switch cfg.sourceprog,
    case 'WINRIVER',
     adcp.nav_longitude(k)=nmean(ens.slongitude);
     adcp.nav_latitude(k)=nmean(ens.slatitude);  
   case 'VMDAS',
     adcp.nav_smtime(k)   =ens.smtime(1);
     adcp.nav_emtime(k)   =ens.emtime(end);
     adcp.nav_slatitude(k)=ens.slatitude(1);
     adcp.nav_elatitude(k)=ens.elatitude(end);
     adcp.nav_slongitude(k)=ens.slongitude(1);
     adcp.nav_elongitude(k)=ens.elongitude(end);
  end;

  % Expand Variables to Match Number of Records!

  if k == n_expand
     for f = fieldnames(adcp)'
         vf = getfield(adcp,f{1});
         si = size(vf);
         ii = ( si == n_expand );
         if any(ii)
            ii = max(find(ii));
            si(ii) = n - n_expand;
            adcp = setfield(adcp,f{1},cat(ii,vf,NaN*ones(si)));
         end
     end
  end

end;  

fprintf('\n');
fclose(fd);

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function hdr = rd_hdr(fid,init);

% Read config data

cfgid = fread(fid,1,'uint16');

if ~( cfgid == hex2dec(init) )
    error(['File ID is ' dec2hex(cfgid) ' not 7F7F - data corrupted or not a BB/WH raw file?']);
end;

hdr = rd_hdrseg(fid);

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [hdr,nbyte] = rd_hdrseg(fid);

% Reads a Header

nbyte = fread(fid,1,'int16'); 

        fseek(fid,1,'cof');

ndat = fread(fid,1,'int8');

offs = fread(fid,ndat,'int16');

hdr = struct('nbyte',{nbyte} , ...
             'dat_offsets',{offs} );

nbyte = 4 + ndat*2;

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function cfg = rd_cfg(fid);

% Read config data

cid = fread(fid,1,'uint16');

if ~( cid == hex2dec('0000') )
    warning(['Fixed header ID ' cid 'incorrect - data corrupted or not a BB/WH raw file?']);
end 

cfg = rd_fixseg(fid);


%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function opt = getopt(val,varargin);

% Returns one of a list (0=first in varargin, etc.)

if val+1 > length(varargin)
   opt = 'unknown';
else
   opt = varargin{val+1};
end
   			

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [cfg,nbyte] = rd_fixseg(fid);

% Reads the configuration data from the fixed leader


pos = ftell(fid);

cfg.name       ='wh-adcp';
cfg.sourceprog = 'instrument';  % default - depending on what data blocks are
                              % around we can modify this later in rd_buffer.
cfg.prog_ver   = fread(fid,1,'uint8')+fread(fid,1,'uint8')/100;

if any( fix(cfg.prog_ver) == [ 4  5 ] )
    cfg.name='bb-adcp';
elseif any( fix(cfg.prog_ver) == [ 8 16 ] )
    cfg.name='wh-adcp';
else
    cfg.name='unrecognized firmware version'   ;    
end;    

cc = fread(fid,2,'uint8');  % Coded stuff

cfg.config          =[ dec2base(cc(2),2,8) '-' dec2base(cc(1),2,8) ];
 cfg.beam_angle     =getopt(bitand(cc(2),3),15,20,30);
 cfg.beam_freq      =getopt(bitand(cc(1),7),75,150,300,600,1200,2400);
 cfg.beam_pattern   =getopt(bitand(cc(1),8)==8,'concave','convex'); % 1=convex,0=concave
% cfg.orientation   =getopt(bitand(cc(1),128)==128,'down','up');   % 1=up,0=down
 cfg.orientation    = ( bitand(cc(1),128) == 128 );    % 1=up,0=down

cfg.simflag        =getopt(fread(fid,1,'uint8'),'real','simulated'); % Flag for simulated data
fseek(fid,1,'cof'); 
cfg.n_beams        =fread(fid,1,'uint8');
cfg.n_cells        =fread(fid,1,'uint8');
cfg.pings_per_ensemble=fread(fid,1,'uint16');
cfg.cell_size      =fread(fid,1,'uint16')*.01;	 % meters
cfg.blank          =fread(fid,1,'uint16')*.01;	 % meters
cfg.prof_mode      =fread(fid,1,'uint8');         %
cfg.corr_threshold =fread(fid,1,'uint8');
cfg.n_codereps     =fread(fid,1,'uint8');
cfg.min_pgood      =fread(fid,1,'uint8');
cfg.evel_threshold =fread(fid,1,'uint16');
cfg.time_between_ping_groups=sum(fread(fid,3,'uint8').*[60 1 .01]'); % seconds

cc      =fread(fid,1,'uint8');                                % Lots of bit-mapped info
  cfg.coord=dec2base(cc,2,8);
  cfg.coord_sys      =getopt(bitand(bitshift(cc,-3),3),'beam','instrument','ship','earth');
  cfg.use_pitchroll  =getopt(bitand(cc,4)==4,'no','yes');  
  cfg.use_3beam      =getopt(bitand(cc,2)==2,'no','yes');
  cfg.bin_mapping    =getopt(bitand(cc,1)==1,'no','yes');

cfg.xducer_misalign=fread(fid,1,'int16')*.01;    % degrees
cfg.magnetic_var   =fread(fid,1,'int16')*.01;	% degrees
cfg.sensors_types  ='-cdhprst';
cfg.sensors_src    =dec2base(fread(fid,1,'uint8'),2,8);
cfg.sensors_avail  =dec2base(fread(fid,1,'uint8'),2,8);
cfg.bin1_dist      =fread(fid,1,'uint16')*.01;	% meters
cfg.xmit_pulse     =fread(fid,1,'uint16')*.01;	% meters
cfg.water_ref_cells=fread(fid,2,'uint8');
cfg.fls_target_threshold =fread(fid,1,'uint8');
fseek(fid,1,'cof');
cfg.xmit_lag       =fread(fid,1,'uint16')*.01; % meters


if cfg.prog_ver >= 8.14,  % Added CPU serial number with v8.14
  cfg.serialnum = fread(fid,8,'uint8');
end;

if cfg.prog_ver >= 8.24,  % Added 2 more bytes with v8.24 firmware
  cfg.sysbandwidth = fread(fid,2,'uint8');
end;

if cfg.prog_ver >= 16.05,                      % Added 1 more bytes with v16.05 firmware
  cfg.syspower  = fread(fid,1,'uint8');
end;

% It is useful to have this precomputed.

cfg.ranges = cfg.bin1_dist+[0:cfg.n_cells-1]' * cfg.cell_size;

%%% if cfg.orientation==1, cfg.ranges=-cfg.ranges; end

nbyte = ftell(fid) - pos;

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ens,hdr,cfg,ok] = rd_buffer(fd,init,num_av,ens_num);

ok = 1;

% To save it being re-initialized every time.

global ens hdr

% A fudge to try and read files not handled quite right.

global FIXOFFSET SOURCE

% If num_av < 0 we are reading only 1 element and initializing
if num_av < 0 | isempty(ens),

  FIXOFFSET=0;   
  SOURCE=0;  % 0=instrument, 1=VMDAS, 2=WINRIVER

  n   = abs(num_av);
  pos = ftell(fd);
  hdr = rd_hdr(fd,init);
  cfg = rd_cfg(fd);

  fseek(fd,pos,'bof');
  clear global ens
  global ens
 
 ens=struct('number',zeros(1,n),'rtc',zeros(7,n),'BIT',zeros(1,n),'ssp',zeros(1,n),'depth',zeros(1,n),'pitch',zeros(1,n),...
            'roll',zeros(1,n),'heading',zeros(1,n),'temperature',zeros(1,n),'salinity',zeros(1,n),...
            'mpt',zeros(1,n),'heading_std',zeros(1,n),'pitch_std',zeros(1,n),...
            'roll_std',zeros(1,n),'adc',zeros(8,n),'error_status_wd',zeros(1,n),...
            'pressure',zeros(1,n),'pressure_std',zeros(1,n),...
            'east_vel',zeros(cfg.n_cells,n),'north_vel',zeros(cfg.n_cells,n),'vert_vel',zeros(cfg.n_cells,n),...
            'error_vel',zeros(cfg.n_cells,n),'intens',zeros(cfg.n_cells,4,n),'percent',zeros(cfg.n_cells,4,n),...
            'corr',zeros(cfg.n_cells,4,n),'status',zeros(cfg.n_cells,4,n),'bt_range',zeros(4,n),'bt_vel',zeros(4,n),...
            'smtime',zeros(1,n),'emtime',zeros(1,n),'slatitude',zeros(1,n),...
	        'slongitude',zeros(1,n),'elatitude',zeros(1,n),'elongitude',zeros(1,n),...
	        'flags',zeros(1,n));
  num_av=abs(num_av);
end;

k=0;
while k<num_av,
   
   
   id1=dec2hex(fread(fd,1,'uint16'));

   if ~strcmp(id1,'7F7F'),
	if isempty(id1),  % End of file
	 ens=-1;
	 return;
	end;    
        error(['Not a workhorse/broadband file or bad data encountered: ->' id1]); 
   end;

   startpos=ftell(fd)-2;  % Starting position.
   
   
   % Read the # data types.
   [hdr,nbyte]=rd_hdrseg(fd);      
   byte_offset=nbyte+2;

   % Read all the data types.

  len = length(hdr.dat_offsets);

  ids = -1 * ones(len,1);

  %*******************************************************************************
  for n = 1 : len
  %*******************************************************************************

    fseek(fd,startpos+hdr.dat_offsets(n),'bof');
  
    id=fread(fd,1,'uint16');

    ok = ( ok & ~any( id == ids ) );

    ids(n) = id;

    % handle all the various segments of data. Note that since I read the IDs as a two
    % byte number in little-endian order the high and low bytes are exchanged compared to
    % the values given in the manual.
    %

if 0 %%% ens_num >= 5690
   fprintf(1,'\n%4.0f/%2.0f:  "%s" #%4.4d *%4.4d >%6.6d ...',ens_num,n,dec2hex(id,4),nbyte,byte_offset,ftell(fd));
end

    %*****************************************************************
    switch dec2hex(id,4),           
    %*****************************************************************

    %----------------------------------------------------------------- 
    case '0000',   % Fixed leader
    %----------------------------------------------------------------- 

      [cfg,nbyte]=rd_fixseg(fd);
      nbyte=nbyte+2;
      
    %----------------------------------------------------------------- 
    case '0080'   % Variable Leader
    %----------------------------------------------------------------- 

      k=k+1;
      ens.number(k)         =fread(fd,1,'uint16');
      ens.rtc(:,k)          =fread(fd,7,'uint8');
      ens.number(k)         =ens.number(k)+65536*fread(fd,1,'uint8');
      ens.BIT(k)            =fread(fd,1,'uint16');
      ens.ssp(k)            =fread(fd,1,'uint16');
      ens.depth(k)          =fread(fd,1,'uint16')*.1;   % meters
      ens.heading(k)        =fread(fd,1,'uint16')*.01;  % degrees
      ens.pitch(k)          =fread(fd,1,'int16')*.01;   % degrees
      ens.roll(k)           =fread(fd,1,'int16')*.01;   % degrees
      ens.salinity(k)       =fread(fd,1,'int16');       % PSU
      ens.temperature(k)    =fread(fd,1,'int16')*.01;   % Deg C
      ens.mpt(k)            =sum(fread(fd,3,'uint8').*[60 1 .01]'); % seconds
      ens.heading_std(k)    =fread(fd,1,'uint8');     % degrees
      ens.pitch_std(k)      =fread(fd,1,'int8')*.1;   % degrees
      ens.roll_std(k)       =fread(fd,1,'int8')*.1;   % degrees
      ens.adc(:,k)          =fread(fd,8,'uint8');
      nbyte=2+40;

      if strcmp(cfg.name,'bb-adcp'),
      
          if cfg.prog_ver>=5.55,
              fseek(fd,15,'cof'); % 14 zeros and one byte for number WM4 bytes
	          cent=fread(fd,1,'uint8');            % possibly also for 5.55-5.58 but
	          ens.rtc(:,k)=fread(fd,7,'uint8');    % I have no data to test.
	          ens.rtc(1,k)=ens.rtc(1,k)+cent*100;
	          nbyte=nbyte+15+8;
		  end;
          
      elseif strcmp(cfg.name,'wh-adcp'), % for WH versions.		

          ens.error_status_wd(k)=fread(fd,1,'uint32');
          nbyte=nbyte+4;;

	      if cfg.prog_ver>=8.13,  % Added pressure sensor stuff in 8.13
                  fseek(fd,2,'cof');   
                  ens.pressure(k)       =fread(fd,1,'uint32');  
                  ens.pressure_std(k)   =fread(fd,1,'uint32');
	          nbyte=nbyte+10;  
	      end;

	      if cfg.prog_ver>8.24,  % Spare byte added 8.24
	          fseek(fd,1,'cof');
	          nbyte=nbyte+1;
	      end;

	      if cfg.prog_ver>=16.05,   % Added more fields with century in clock 16.05
	          cent=fread(fd,1,'uint8');            
	          ens.rtc(:,k)=fread(fd,7,'uint8');   
	          ens.rtc(1,k)=ens.rtc(1,k)+cent*100;
	          nbyte=nbyte+8;
	      end;
      end;
  	      
    %----------------------------------------------------------------- 
    case '0100',  % Velocities
    %----------------------------------------------------------------- 

      vels=fread(fd,[4 cfg.n_cells],'int16')'*.001;     % m/s
      ens.east_vel(:,k) =vels(:,1);
      ens.north_vel(:,k)=vels(:,2);
      ens.vert_vel(:,k) =vels(:,3);
      ens.error_vel(:,k)=vels(:,4);
      nbyte=2+4*cfg.n_cells*2;
      
    %----------------------------------------------------------------- 
    case '0200',  % Correlations
    %----------------------------------------------------------------- 

      ens.corr(:,:,k)   =fread(fd,[4 cfg.n_cells],'uint8')';
      nbyte=2+4*cfg.n_cells;
      
    %----------------------------------------------------------------- 
    case '0300',  % Echo Intensities  
    %----------------------------------------------------------------- 

      ens.intens(:,:,k)   =fread(fd,[4 cfg.n_cells],'uint8')';
      nbyte=2+4*cfg.n_cells;

    %----------------------------------------------------------------- 
    case '0400',  % Percent good
    %----------------------------------------------------------------- 

      ens.percent(:,:,k)   =fread(fd,[4 cfg.n_cells],'uint8')';
      nbyte=2+4*cfg.n_cells;
   
    %----------------------------------------------------------------- 
    case '0500',  % Status
    %----------------------------------------------------------------- 

         % Note in one case with a 4.25 firmware SC-BB, it seems like
         % this block was actually two bytes short!
      ens.status(:,:,k)   =fread(fd,[4 cfg.n_cells],'uint8')';
      nbyte=2+4*cfg.n_cells;

    %----------------------------------------------------------------- 
    case '0600', % Bottom track
    %----------------------------------------------------------------- 

                 % In WINRIVER GPS data is tucked into here in odd ways, as long
                 % as GPS is enabled.
      if SOURCE==2,
          fseek(fd,2,'cof');
          long1=fread(fd,1,'uint16');
          fseek(fd,6,'cof');           
          cfac=180/2^31;
          ens.slatitude(k)  =fread(fd,1,'int32')*cfac;
      else    
          fseek(fd,14,'cof'); % Skip over a bunch of stuff
      end;    
      ens.bt_range(:,k)=fread(fd,4,'uint16')*.01; %
      ens.bt_vel(:,k)  =fread(fd,4,'int16');
      if SOURCE==2,
          fseek(fd,12+2,'cof');
          ens.slongitude(k)=(long1+65536*fread(fd,1,'uint16'))*cfac;
          if ens.slongitude(k)>180, ens.slongitude(k)=ens.slongitude(k)-360; end;
          fseek(fd,71-33-16,'cof');
          nbyte=2+68; 
      else    
          fseek(fd,71-33,'cof');
          nbyte=2+68;
      end;    
      if cfg.prog_ver>=5.3,    % Version 4.05 firmware seems to be missing these last 11 bytes.
       fseek(fd,78-71,'cof');  
       ens.bt_range(:,k)=ens.bt_range(:,k)+fread(fd,4,'uint8')*655.36;
       nbyte=nbyte+11;
       if cfg.prog_ver>=16,   % RDI documentation claims these extra bytes were added in v 8.17
           fseek(fd,4,'cof');  % but they don't appear in my 8.33 data.
           nbyte=nbyte+4;
       end;
      end;
     
% The raw files produced by VMDAS contain a binary navigation data
% block. 
      
    %----------------------------------------------------------------- 
    case '2000',  % Something from VMDAS.
    %----------------------------------------------------------------- 

      cfg.sourceprog='VMDAS';
      SOURCE=1;
      utim  =fread(fd,4,'uint8');
      mtime =datenum(utim(3)+utim(4)*256,utim(2),utim(1));
      ens.smtime(k)     =mtime+fread(fd,1,'uint32')/8640000;
      fseek(fd,4,'cof');
      cfac=180/2^31;
      ens.slatitude(k)  =fread(fd,1,'int32')*cfac;
      ens.slongitude(k) =fread(fd,1,'int32')*cfac;
      ens.emtime(k)     =mtime+fread(fd,1,'uint32')/8640000;
      ens.elatitude(k)  =fread(fd,1,'int32')*cfac;
      ens.elongitude(k) =fread(fd,1,'int32')*cfac;
      fseek(fd,12,'cof');
      ens.flags(k)      =fread(fd,1,'uint16');	
      fseek(fd,30,'cof');
      nbyte=2+76;
       
% The following blocks come from WINRIVER files, they aparently contain
% the raw NMEA data received from a serial port.
%
% Note that for WINRIVER files somewhat decoded data is also available
% tucked into the bottom track block.
    
    %----------------------------------------------------------------- 
    case '2100', % $xxDBT  (Winriver addition) 38
    %----------------------------------------------------------------- 

      cfg.sourceprog='WINRIVER';
      SOURCE=2;
      str=fread(fd,38,'uchar')';
      nbyte=2+38;

    %----------------------------------------------------------------- 
    case '2101', % $xxGGA  (Winriver addition) 94 in maanual but 97 seems to work
    %----------------------------------------------------------------- 

      cfg.sourceprog='WINRIVER';
      SOURCE=2;
      str=fread(fd,97,'uchar')';
      nbyte=2+97;
      
    %----------------------------------------------------------------- 
    case '2102', % $xxVTG  (Winriver addition) 45
    %----------------------------------------------------------------- 

      cfg.sourceprog='WINRIVER';
      SOURCE=2;
      str=fread(fd,45,'uchar')';
      nbyte=2+45;
      
    %----------------------------------------------------------------- 
    case '2103', % $xxGSA  (Winriver addition) 60
    %----------------------------------------------------------------- 

      cfg.sourceprog='WINRIVER';
      SOURCE=2;
      str=fread(fd,60,'uchar')';
      nbyte=2+60;

    %----------------------------------------------------------------- 
    case '2104',  %xxHDT or HDG (Winriver addition) 38
    %----------------------------------------------------------------- 

      cfg.sourceprog='WINRIVER';
      SOURCE=2;
      str=fread(fd,38,'uchar')';
      nbyte=2+38;
      
      
        
    %----------------------------------------------------------------- 
    case '0701', % Number of good pings
    %----------------------------------------------------------------- 

      fseek(fd,4*cfg.n_cells,'cof');
      nbyte=2+4*cfg.n_cells;
    
    %----------------------------------------------------------------- 
    case '0702', % Sum of squared velocities
    %----------------------------------------------------------------- 

      fseek(fd,4*cfg.n_cells,'cof');
      nbyte=2+4*cfg.n_cells;

    %----------------------------------------------------------------- 
    case '0703', % Sum of velocities      
    %----------------------------------------------------------------- 

      fseek(fd,4*cfg.n_cells,'cof');
      nbyte=2+4*cfg.n_cells;

% These blocks were implemented for 5-beam systems

    %----------------------------------------------------------------- 
    case '0A00', % Beam 5 velocity (not implemented)
    %----------------------------------------------------------------- 

      fseek(fd,cfg.n_cells,'cof');
      nbyte=2+cfg.n_cells;

    %----------------------------------------------------------------- 
    case '0301', % Beam 5 Number of good pings (not implemented)
    %----------------------------------------------------------------- 

      fseek(fd,cfg.n_cells,'cof');
      nbyte=2+cfg.n_cells;

    %----------------------------------------------------------------- 
    case '0302', % Beam 5 Sum of squared velocities (not implemented)
    %----------------------------------------------------------------- 

      fseek(fd,cfg.n_cells,'cof');
      nbyte=2+cfg.n_cells;
             
    %----------------------------------------------------------------- 
    case '0303', % Beam 5 Sum of velocities (not implemented)
    %----------------------------------------------------------------- 

      fseek(fd,cfg.n_cells,'cof');
      nbyte=2+cfg.n_cells;
             
    %----------------------------------------------------------------- 
    case '020C', % Ambient sound profile (not implemented)
    %----------------------------------------------------------------- 

      fseek(fd,4,'cof');
      nbyte=2+4;
             
    %----------------------------------------------------------------- 
    otherwise,
    %----------------------------------------------------------------- 
      
      fprintf('Unrecognized ID code: %s\n',dec2hex(id,4));
      nbyte=2;

      
    %*****************************************************************
    end
    %*****************************************************************
   
    % here I adjust the number of bytes so I am sure to begin
    % reading at the next valid offset. If everything is working right I shouldn't have
    % to do this but every so often firware changes result in some differences.

    byte_offset=byte_offset+nbyte;   
      
%% fprintf(1,' #%4.4d *%4.4d >%6.6d\n\n',nbyte,byte_offset,ftell(fd));

    if n<length(hdr.dat_offsets),
      if hdr.dat_offsets(n+1)~=byte_offset
        %fprintf('%.0f/%s: Adjust location by %d\n',ens_num,dec2hex(id,4),hdr.dat_offsets(n+1)-byte_offset);
        fseek(fd,hdr.dat_offsets(n+1)-byte_offset,'cof');
        ok = 0;
      end;	
      byte_offset = hdr.dat_offsets(n+1); 
    end;

  %*******************************************************************************
  end
  %*******************************************************************************

  % Now at the end of the record we have two reserved bytes, followed
  % by a two-byte checksum = 4 bytes to skip over.

  readbytes=ftell(fd)-startpos;
  offset=(hdr.nbyte+2)-byte_offset; % The 2 is for the checksum
  if offset ~=4 & FIXOFFSET==0, 
	fprintf('\n*****************************************************\n');
    fprintf('Adjust location by %d (readbytes=%d, hdr.nbyte=%d)\n',offset,readbytes,hdr.nbyte);
    fprintf(' NOTE - THIS IS A PROGRAM PROBLEM, POSSIBLY FIXED BY A FUDGE\n');
    fprintf('        PLEASE REPORT TO rich@ocgy.ubc.ca WITH DETAILS!!\n');
    fprintf('        ATTEMPTING TO RECOVER...\n');
    fprintf('******************************************************\n');
    FIXOFFSET=offset-4;
  end;  
  fseek(fd,4+FIXOFFSET,'cof'); 
   
  % An early version of WAVESMON and PARSE contained a bug which stuck an additional two
  % bytes in these files, but they really shouldn't be there 
  %if cfg.prog_ver>=16.05,    
  %	  fseek(fd,2,'cof');
  %end;
  	   
end;

% Blank out stuff bigger than error velocity
% big_err=abs(ens.error_vel)>.2;
big_err=0;
	
% Blank out invalid data	
ens.east_vel(ens.east_vel==-32.768 | big_err)=NaN;
ens.north_vel(ens.north_vel==-32.768 | big_err)=NaN;
ens.vert_vel(ens.vert_vel==-32.768 | big_err)=NaN;
ens.error_vel(ens.error_vel==-32.768 | big_err)=NaN;




%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function y=nmedian(x,window,dim);

% Copied from median but with handling of NaN different.

if nargin==2, 
  dim = min(find(size(x)~=1)); 
  if isempty(dim), dim = 1; end
end

siz = [size(x) ones(1,dim-ndims(x))];
n = size(x,dim);

% Permute and reshape so that DIM becomes the row dimension of a 2-D array
perm = [dim:max(length(size(x)),dim) 1:dim-1];
x = reshape(permute(x,perm),n,prod(siz)/n);

% Sort along first dimension
x = sort(x,1);
[n1,n2]=size(x);

if n1==1,
 y=x;
else
  if n2==1,
   kk=sum(finite(x),1);
   if kk>0,
     x1=x(max(fix(kk/2),1));
     x2=x(max(ceil(kk/2),1));
     x(abs(x-(x1+x2)/2)>window)=NaN;
   end;
   x = sort(x,1);
   kk=sum(finite(x),1);
   x(isnan(x))=0;
   y=NaN;
   if kk>0,
    y=sum(x)/kk;
   end;
  else
   kk=sum(finite(x),1);
   ll=kk<n1-2;
   kk(ll)=0;x(:,ll)=NaN;
   x1=x(max(fix(kk/2),1)+[0:n2-1]*n1);
   x2=x(max(ceil(kk/2),1)+[0:n2-1]*n1);

   x(abs(x-ones(n1,1)*(x1+x2)/2)>window)=NaN;
   x = sort(x,1);
   kk=sum(finite(x),1);
   x(isnan(x))=0;
   y=NaN+ones(1,n2);
   if any(kk),
    y(kk>0)=sum(x(:,kk>0))./kk(kk>0);
   end;
  end;
end; 

% Permute and reshape back
siz(dim) = 1;
y = ipermute(reshape(y,siz(perm)),perm);


%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function y=nmean(x,dim);

% R_NMEAN Computes the mean of matrix ignoring NaN
%         values
%   R_NMEAN(X,DIM) takes the mean along the dimension DIM of X. 
%

kk=finite(x);
x(~kk)=0;

if nargin==1, 
  % Determine which dimension SUM will use
  dim = min(find(size(x)~=1));
  if isempty(dim), dim = 1; end
end;

if dim>length(size(x)),
 y=x;              % For matlab 5.0 only!!! Later versions have a fixed 'sum'
else
  ndat=sum(kk,dim);
  indat=ndat==0;
  ndat(indat)=1; % If there are no good data then it doesn't matter what
                 % we average by - and this avoid div-by-zero warnings.

  y = sum(x,dim)./ndat;
  y(indat)=NaN;
end;

























