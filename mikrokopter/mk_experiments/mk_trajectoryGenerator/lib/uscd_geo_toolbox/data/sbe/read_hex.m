function [msg,head,data,init,fig] = read_hex(file,varargin)

% READ_HEX  Reads SBE-Ascii-DataFiles with hexadecimal Values (HEX)
%
% [ Msg , Head , Data , Init ] = READ_HEX( FileName )
%
% Adjusted for SeaCat-SBE16plus-DataFiles, verify the Subfunctions
%  which are listen below with the CalibrationSheets for correct formulas.
% 
% Head =  HeaderStructure from READ_SBE
%
% Data = [ Day  T  C  P  S  [V0]  [V1]  [V2]  [V3]  [SBE38T] ]
%
% INIT = { VarID  Unit  Tag  Format }
%
%            Variable           Unit    SubFunc    Coefficients from HEX-File
%
% Day        Matlab DATENUM    [day]    hex2time
% T          Temperature       [degC]   hex2temp   TA0, TA1, TA2, TA3
% C          Conductivity      [mS/cm]  hex2cond   G, H, I, J, CTCOR, CPCOR
% P          Pressure                   hex2pres   !!! not complete yet !!!
%                                       hex2qres   !!! not complete yet !!!
% V0 .. V3   ExternalVoltages  [mV]     hex2volt   
%
% SBE38T     SBE38-Temperature [degC]   hex2sbe38  ???
% 
% To convert the external voltages of additional Sensors, 
%  add/edit the following Coefficients to the HEX-File before !!!
%
% Seapoint   Chl_a             [ug/l]   hex2chla   GAIN, SENS
% SBE43      Oxygen      [%] | [ml/l]   hex2oxyg   SOC, BOC, VOFFSET, TCOR, PCOR
%
%
% Use additional String-Inputs to specify the external Sensors 
%  and their Order in the DataFile:
%
% READ_HEX( FileName , EXT1 , EXT2 , ... )
%
%   SBE43:    'O' | 'X' | 'Oxygen'                   (PREF required, see below)
%   Seapoint: 'F' | 'C' | 'Fluorometer | 'Chlorphyll'
%
% For not specified additional Sensors the Voltage [mV] of the Channel is returned.
% See also: READ_SBE, SBE_COND
%
% Valid Limits for the Data can be given by a Structure as last Input:
%
% READ_HEX( FileName , ... , LIM )
%
% Fields of Structure LIM:
%
%    pres: [Pmin Pmax]  dbar
%    temp: [Tmin Tmax]  degC
%    cond: [Cmin Cmax]  mS/sm
%    salt: [Smin Smax]  PSU
%    chla: [Fmin Fmax]  ug/l  Seapoint Fluorecense
%    oxyg: [Omin Omax]  %     Percent of OxygenSaturation
%
% Values outside this Limits are set to NaN. 
% Check the begin of this M-File for the default Limits!!!
%
%----------------------------------------------------------------------------
% Compressibility Compensation of Sea-Bird Conductivity Sensors
%
% For Instruments without a build-in pressure sensor you can set
%  the reference pressure for the compressibility compensation of 
%  the conductivity sensor:  READ_HEX( FileName , PREF , ... )
%
% By default the setup reference pressure from the header or ZERO is used.
%
% To correct the conductivity by a known pressure use the formula:
%
% C = C * ( 1 + CTCOR * T + CPCOR  * PREF ) / ( 1 + CTCOR * T + CPCOR  * P )
%
% defaults:  CTCOR = 3.25e-06;  CPCOR = -9.57e-08
%
% Reference: 
%
% Seabird Electonics Application Note 10,
% "Compressibility Compensation of Sea-Bird Conductivity Sensors",
% http://www.seabird.com/application_notes/AN10.htm
% 
%----------------------------------------------------------------------------
%
% A NONZERO IMAGinary part of PREF plots the TimeSeries!!!
%
% [ Msg , Header , Data , HeadINI , FIG ] = READ_HEX( FileName , PREF+i , ... )
%
%----------------------------------------------------------------------------
% Structure of HEX-DataLine (SeaCat SBE 16plus):
%
% tttttt    6  Temperature A/D-counts
% cccccc    6  CondFreq          = val / 256
%
% vvvv      4  PressTempCompVolt = val / 13.107  % PressureType 1 (StrainGauge)
%
% pppppp    6  PressFreq         = val / 256     % PressureType 3 (Quartz)
% vvvv      4  PressTempCompVolt = val / 13.107  % PressureType 3 (Quartz) 
%
% vvvv      4  Ext0Volt          = val / 13.107  % [mV]
% vvvv      4  Ext1Volt          = val / 13.107
% vvvv      4  Ext2Volt          = val / 13.107
% vvvv      4  Ext3Volt          = val / 13.107
%
% tttttt    6  SBE38Temp         = val / 100 - 10
%
% ssssssss  8  seconds since 01.01.1980
%
%          20 .. 40/46/52 Counts
%

% Default Limits !!!

lim = struct( 'pres' , [ 0 6000] , ...
              'temp' , [ 0  30 ] , ...
              'cond' , [ 25 60 ] , ...
              'salt' , [ 30 40 ] , ...
              'chla' , [ 0  10 ] , ...
              'oxyg' , [ 0 150 ]       );


fig  = [];
init = cell(0,3);

%***************************************************************
% Get Header from DataFile

data = [];

[msg,head,ini] = read_sbe(file,i);

if ~isempty(msg) | ( nargout < 3 )
    return
end

c = head.proc;

frm = 'raw HEX';

ii = strcmp(c(:,1),'output format');
if any(ii)
   ii = find(ii);
   of = c{ii(end),2};
   if ~strcmp(of,frm)
       msg = sprintf('Invalid OutputFormat "%s", required: "%s".',of,frm);
       return
   end
end
 
%***************************************************************
% Check other Inputs

Nin = prod(size(varargin));

pref = [];
ext  = '';

cc  = 'OF';  % !!! Valid add. Sensors !!!!!!

msg = cell(0,1);

if Nin > 0
   if isstruct(varargin{Nin})
      [m,lim] = structcmp(lim,varargin{end},0);
      if ~isempty(m)
          warning('Invalid Structure for Limits.')
      end
      Nin = Nin - 1;
   end
end


for ii = 1 : Nin
    v = varargin{ii};
    is_num = isnumeric(v);
    is_str = chkstr(v);
    m = '';
    if ~( is_num | is_str )
        m = 'Add. inputs must be single numerics or strings.';
    elseif ~isempty(v)
        if is_num 
           if ~( prod(size(v)) == 1 )
               m = 'Value for PREF must be a single numeric.';
           elseif ~( isfinite(v) & ( v >= 0 ) )
               m = 'Value for PREF must be positive finite.';
           elseif ~isempty(pref)
               m = 'Duplicate PREF in Inputs.';
           else
               pref = v; 
           end
        elseif is_str
           v1 = upper(v(1)); 
           if v1 == 'X', v1 == 'O'; end
           if v1 == 'C', v1 == 'F'; end
           if ~any( v1 == cc )
               m = sprintf('Unknown Sensor "%s", known: %s.',v,cc);
           else
               ext = cat(2,ext,v1);
           end
        end
    end
    if ~isempty(m)
        msg = cat(1,msg,{sprintf('%.0f: %s',ii+1,m)});
    end
end

ext = cat(2,ext,' ');

next = sum(~(ext == ' '));

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    msg = sprintf('Invalid Inputs.\n%s',msg);
    return
end

%***************************************************************
% Check for ReferencePressure in Header

if isempty(pref)
   pref = 0;
   ii = ( strcmp( lower(c(:,1)) , 'reference pressure' ) | ...
          strcmp( lower(c(:,1)) , 'refpress' )         ); 
   if any(ii)
      ii = find(ii);
      ii = ii(:)';
      for jj = ii(end:-1:1)
          p  = sscanf(c{jj,2},'%f');
          if prod(size(p)) == 1
             pref = p;
             break
          end
      end
   end
end

ShowResult = ~( imag(pref) == 0 );

pref = real(pref);

%***************************************************************
% Check for Sensors by Coefficients in Header

tok = hex2temp(c);
cok = hex2cond(c);
pok = hex2pres(c);
qok = hex2qres(c);
fok = hex2chla(c);
xok = hex2oxyg(c);

fok = ( fok * any( ext == 'F' ) );
xok = ( xok * any( ext == 'O' ) );

m = cell(0,1);

if tok == -1, m = cat(1,m,{'Uncomplete Coeff for Temp.'}); end
if cok == -1, m = cat(1,m,{'Uncomplete Coeff for Cond.'}); end
if pok == -1, m = cat(1,m,{'Uncomplete Coeff for Press.'}); end
if qok == -1, m = cat(1,m,{'Uncomplete Coeff for Quartz.'}); end
if fok == -1, m = cat(1,m,{'Uncomplete Coeff for Fluorometer.'}); end
if xok == -1, m = cat(1,m,{'Uncomplete Coeff for Oxygen.'}); end

if pok & qok, cat(1,m,{'Found 2 PressureSensors.'}); end

if any( ext == 'F' ) & ( fok == 0 )
   m = cat(1,m,{'Can''t find Coeff for Fluorometer.'}); 
end

if any( ext == 'O' ) & ( xok == 0 )
   m = cat(1,m,{'Can''t find Coeff for Oxygen.'}); 
end

if ~isempty(m)
    msg = sprintf('%s\n',m{:});
    msg = sprintf('Invalid Header.\n%s',msg);
    return
end

%***************************************************************
% Read DataLines

n0 = 6*tok + 6*cok + 4*pok + 10*qok + 8;

[msg,d] = readhxl(file,head.bytes,n0);

if isempty(msg)
   if isempty(d)
      msg = 'No Data.';
   end
end

if ~isempty(msg)
    return
end

%***************************************************************
% Get Parameters
%
% tttttt    6  Temperature A/D-counts
% cccccc    6  CondFreq          = val / 256
%
% vvvv      4  PressTempCompVolt = val / 13.107  % PressureType 1 (StrainGauge)
%
% pppppp    6  PressFreq         = val / 256     % PressureType 3 (Quartz)
% vvvv      4  PressTempCompVolt = val / 13.107  % PressureType 3 (Quartz) 
%
% vvvv      4  Ext0Volt          = val / 13.107
% vvvv      4  Ext1Volt          = val / 13.107
% vvvv      4  Ext2Volt          = val / 13.107
% vvvv      4  Ext3Volt          = val / 13.107
%
% tttttt    6  SBE38Temp         = val / 100 - 10
%
% ssssssss  8  seconds since 01.01.1980
%
%          20 .. 48 / 54 Counts
%


  n = size(d,2);

  m = n - n0;   % Residue 

sok = ( ( m >= 6 ) & ( mod(m,4) == 2 ) );  % True for SBE38Temp

  m = m - 6*sok;

vok = max( 0 , floor(m/4) );  % Number of ExternalVoltage

fprintf(1,'Sensors: ');
if tok, fprintf(1,'T '); end
if cok, fprintf(1,'C '); end
if pok, fprintf(1,'P '); end
if qok, fprintf(1,'Q '); end

if vok, 
   for ii = 1 : vok
       fprintf(1,'V%.0f',ii-1);
       if next >= ii
          if     fok & ( ext(ii) == 'F' )
                 fprintf(1,'(Seapoint Fluorom.)');
          elseif xok & ( ext(ii) == 'O' )
                 fprintf(1,'(SBE43 Oxygen)');
          end
       end
       fprintf(1,' ');
   end
end
if sok, fprintf(1,'SBE38T '); end

fprintf(1,'\n');

m =  6*tok + 6*cok + 4*pok + 10*qok + 4*vok + 6*sok + 8;

if ~( n == m );
    msg = 'Invalid Configuration found.';
end

nn = 5 + vok + sok;

data = NaN * zeros(size(d,1),nn);

init = cell(nn,4);
init(:) = {''};

%------------------------------------------------------
% DateNum

data(:,1) = hex2time( d( : , n-8+1 : n ) );

init(1,:) = { 'TIM'  'day'  'Time'  '%.6f' };

z = 0; % HEX-Counter

%------------------------------------------------------
% Temp

t = 0;   % Temp for Press and Cond

init(2,[1 2]) = { 'T'  'degC' };

if tok
   data(:,2) = hex2temp( c , d(:,z+(1:6)) );
   data(:,2) = lim_check(data(:,2),lim.temp);
   init(2,[3 4]) = { 'Temp'  '%6.3f' };
   t = data(:,2);
   z = z + 6;
end

%------------------------------------------------------
% Press before Cond !!!

zc = z;             % Save Counter for Cond
z  = z + 6 * cok;   % Set  Counter for Press and following

if pok
   data(:,4) = hex2pres( c , d(:,z+(1:4)) , t );
   pfrm = '%6.1f';
   z = z + 4;
elseif qok
   data(:,4) = hex2qres( c , d(:,z+6+(1:4)) , d(:,z+(1:6)) , t );
   pfrm = '%7.2f';
   z = z + 10;
end

app = '';

init(4,[1 2]) = { 'P'  'dbar' };

if pok | qok
   data(:,4) = lim_check(data(:,4),lim.pres);
   init(4,[3 4]) = { 'Pres'  pfrm };
   p = data(:,4); 
else
   p = pref; 
   fprintf(1,'Use ReferencePressure: %.1f dbar\n',pref);
   app = sprintf('_%.0f',pref);
end

%------------------------------------------------------
% Cond

cnd = 0;

init(3,[1 2]) = { 'C'  'mS/cm' };

if cok
  data(:,3) = hex2cond( c , d(:,zc+(1:6)) , t , p );
  data(:,3) = lim_check(data(:,3),lim.cond);
  init(3,[3 4]) = { 'Cond' '%7.4f' };
  cnd = data(:,3);
end

%------------------------------------------------------
% Salinity

sal = 0;

init(5,[1 2]) = { 'S'  'PSU' };
 
if tok & cok

   sal = sal78( p , t*1.00024 , cnd , 42.914 , 0 );

   sal = lim_check(sal,lim.salt);

   data(:,5) = sal;

   init(5,[3 4]) = { ['Salt' app]  '%7.4f' };

end

%------------------------------------------------------
% External Voltage

for ii = 1 : vok 
    jj = 5 + ii;
    data(:,jj) = hex2volt(d(:,z+(1:4)));
    lab = sprintf('EXT%.0f',ii);
    uni = 'mV';
    tag = lab;
    frm = '%7.2f';
    if next >= ii
       if     fok & ( ext(ii) == 'F' )
          data(:,jj) = hex2chla( c , data(:,jj) );
          data(:,jj) = lim_check(data(:,jj),lim.chla);
          lab = 'CHL_A'; uni = 'ug/l'; tag = 'Chla'; frm = '%6.2f';
       elseif xok & ( ext(ii) == 'O' ) & tok & cok
          %%% oxs = oxsat(sal,t); uni = 'ml/l'; frm = '5.2f';
          oxs = 100; uni = '%'; frm = '%5.1f';
          data(:,jj) = hex2oxyg( c , data(:,jj) , t , p , oxs );
          data(:,jj) = lim_check(data(:,jj),lim.oxyg);
          lab = 'O'; tag = 'Oxyg';
       end
    end
    init(jj,:) = { lab  uni  tag  frm };
    z = z + 4;
end


%------------------------------------------------------
% SBE38-Temperature

if sok
   data(:,nn) = hex2sbe38(d);
   data(:,nn) = lim_check(data(:,nn),lim.temp);
   init(nn,:) = { 'T'  'degC'  'Temp38'  '%6.3f' };
end

%***************************************************************

if ShowResult
   fig = sc_plot(data,c,init,vok,sok,pref);
   set(fig,'visible','on');
end


%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% SubFunctions to convert the HEX-Values
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%*******************************************************************

function d = hex2time(d)

d = datenum(1980,01,01) + hx2dc(d)/3600/24;

%*******************************************************************

function d = hex2temp(c,d)

[ok,c] = checkcoeff(c,{'TA0' 'TA1' 'TA2' 'TA3'}); % TOFFSET

if nargin < 2
   d = ok;
   return
end

if ~ok
    d = NaN * ones(size(d,1),1);
    return
end

d = ( hx2dc(d) - 524288 ) / 1.6e7;

d = ( d*2.9e9 + 1.024e8 ) ./ ( 2.048e4 - d*2e5 );

d = log(d);

d = ( c.TA0 + c.TA1 * d + c.TA2 * d.^2 + c.TA3 * d.^3 );

d = 1 ./ d - 273.15;

%*******************************************************************

function d = hex2cond(c,d,t,p)

% SBE 52-MP: 
%
% f = instrument frequency (kHz) * (1.0 + WBOTC * t) 0.5 / 1000.00
%

[ok,c] = checkcoeff(c,{'G' 'H' 'I' 'J' 'CTCOR' 'CPCOR' }); % CSLOPE

if nargin < 2
   d = ok;
   return
end

if ~ok
    d = NaN * ones(size(d,1),1);
    return
end

if nargin < 3
   t = 0;
elseif all(isnan(t))
   t = 0;
end

if nargin < 4
   p = 0;
elseif all(isnan(p))
   p = 0;
end

d = hx2dc(d) / 256 / 1000;

d = ( c.G + c.H * d.^2 + c.I * d.^3 + c.J * d.^4 );

d = 10 * d ./ ( 1 + c.CTCOR * t + c.CPCOR  * p );

%*******************************************************************

function d = hex2pres(c,d,t)

% Strain Gauge

p = { 'PA0' 'PA1' 'PA2'              ...
      'PTEMPA0' 'PTEMPA1' 'PTEMPA2'  ...
      'PTCA0' 'PTCA1' 'PTCA2'        ...
      'PTCB0' 'PTCB1' 'PTCB2'             };

[ok,c] = checkcoeff(c,p);

if nargin < 2
   d = ok;
   return
end

if ~ok
    d = NaN * ones(size(d,1),1);
    return
end

if nargin < 3
   t = 0;
elseif all(isnan(t))
   t = 0;
end

d = hx2dc(d) / 13.107;

warning('HEX2PRES: Strain Gauge Pressure Conversion not complete yet.')

%*******************************************************************

function d = hex2qres(c,d,f,t)

% Quartz
%
% pppppp    6  PressFreq         = val / 256     % PressureType 3 (Quartz)
% vvvv      4  PressTempCompVolt = val / 13.107  % PressureType 3 (Quartz) 


p = { 'PC1' 'PC2' 'PC3'       ...
      'PD1' 'PD2'             ...
      'PT1' 'PT2' 'PT3' 'PT4'      }; %  PSLOPE POFFSET EXTFREQSF

[ok,c] = checkcoeff(c,p);

if nargin < 2
   d = ok;
   return
end

if ~ok
    d = NaN * ones(size(d,1),1);
    return
end

[ok,ext] = checkcoeff(c,'EXTFREQSF');

if nargin < 4
   t = 0;
elseif all(isnan(t))
   t = 0;
end

d = hx2dc(d) / 13.107;
f = hx2dc(f) / 256;

warning('HEX2QRES: Quartz Pressure Conversion not complete yet.')

%*******************************************************************

function d = hex2volt(d)

d = hx2dc(d) / 13.107;  % [mV]

%*******************************************************************

function d = hex2sbe38(d)

d = hx2dc(d) / 100 - 10;

%*******************************************************************

function d = hex2chla(c,d)

% Seapoint Fluorometer
%
% Gain 	Sens [V/µg/l]   Range [µg/l]
%
% 30x 	1.0 	5
% 10x 	0.33 	15
% 3x 	0.1 	50
% 1x 	0.033 	150
%

[ok,c] = checkcoeff(c,{'GAIN' 'SENS' }); % CSLOPE

if nargin < 2
   d = ok;
   return
end

if ~ok
    d = NaN * ones(size(d,1),1);
    return
end

if size(d,2) > 1
   d = hex2volt(d);  % [mV]
end

if isempty(c.SENS)
   c.SENS = c.GAIN / 30;
end

d = d / 1000 / c.SENS;

%*******************************************************************

function d = hex2oxyg(c,d,t,p,oxs)

% SBE43 Oxygen Sensor 
%
% Oxygen[perc] = SOC*(V+VOFF) * exp(T*TCOR) * exp(P*PCOR) * 100
% Oxygen[ml/l] = SOC*(V+VOFF) * exp(T*TCOR) * exp(P*PCOR) * OXSAT(T,S)
%
%


Nin = nargin;

[ok,c] = checkcoeff(c,{'SOC' 'BOC' 'VOFFSET' 'TCOR' 'PCOR'});

if Nin < 2
   d = ok;
   return
end

if ~ok
    d = NaN * ones(size(d,1),1);
    return
end

if size(d,2) > 1
   d = hex2volt(d);  % [mV]
end

if Nin < 3
   return
end

if Nin < 4
   p = 0;
elseif all(isnan(p))
   p = 0;
end

if Nin < 5
   oxs = 1;
end

d = ( d/1000 + c.VOFFSET ) * c.SOC;

d = d .* exp(t*c.TCOR) .* exp(p*c.PCOR);

d = d .* oxs ;


%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,d] = checkcoeff(c,p)

% 1 if all ok; -1 if any but not all ok

d = p([1 1],:);

d(2,:) = { {[]} };

d = struct(d{:});

n = prod(size(p));

ok = zeros(1,n);

c(:,1) = upper(c(:,1));
p      = upper(p);

for ii = 1 : n

    jj = strcmp( c(:,1) , p{ii} );
    ok(ii) = any(jj);
 
    if ok(ii)
       jj = min(find(jj));
        v = eval(c{jj,2},'NaN');
       ok(ii) = ( isnumeric(v) & ( prod(size(v)) == 1 ) );
       if ok(ii)
          ok(ii) = isfinite(v);
       end 
       if ok(ii)
          d = setfield( d , p{ii} , v );
       end 
    end

end

ok = 2 * all(ok) - any(ok);

%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [x,ii] = lim_check(x,lim);

ii = [];

if isempty(x)
   return
end

ii = ( ( x < lim(1) ) | ( lim(2) < x ) );

if ~any(ii)
    ii = [];
    return
end

ii = find(ii);

x(ii) = NaN;

%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function d = hx2dc(d);

if isempty(d)
   return
end

n = size(d,2);

for ii = 1 : (n-1)

    d(:,n) = d(:,n) + d(:,ii) * 16^(n-ii);

end

d = d(:,n);

%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,d] = readhxl(file,bytes,n0);

% Read HEX-Lines

msg = '';
d   = [];

%*******************************************

fid = fopen(file,'r');

fseek(fid,bytes,'bof');

if feof(fid)
   return
end

d = fread(fid,'char');

fclose(fid);

%*******************************************
% Basics check with NewLine

%-----------------------------------
% Check NewLines, remove following

d(find(d==13)) = 10;

d = cat( 1 , d , 10 , 10 ); % Append NL

ii = find( d == 10 );
jj = find( diff(ii,1,1) == 1 );
d(ii(jj)) = [];

if all( d == 10 )
   d = [];
   return
end

%-----------------------------------
% Remove leading NewLines

while d(1) == 10
      d = d(2:end);
end

%*******************************************
% Check Characters

nfrm = '%6.0f ';  % Format for LineNumbers, Blank at End !!!

% Make UpperCase Array

aa = double('aA');

d = d - abs(diff(aa)) * ( d >= max(aa) );

%-----------------------------------
% Check for valid Characters [ 0 .. 9 | A .. F ]
 
ok = ( ( ( double('0') <= d ) & ( d <= double('9') ) ) | ...
       ( ( double('A') <= d ) & ( d <= double('F') ) )       );

ii = ( d == 10 );

if ~all( ok | ii )

    if ~any(ok)
        msg = 'No valid Characters.';
        return
    end

    ok = ( ok | ii );

    %-----------------------------------
    % Get Lines

    i0 = cat( 1 , 0 , find(ii) );

     n = diff(i0) - 1;

    i0 = i0(find(i0<size(d,1)));

    ii = zeros(size(d));

    ii(i0+1) = 1;

    ii = cumsum(ii);

    %-----------------------------------
    % Get Mean of LineOk

    [ii,ok] = grpmean(ii,ok);

     ok = ( ok == 1 );  % All Characters Ok

    if ~any(ok)
        msg = 'No valid Lines.';
        return
    end

    %-----------------------------------
    % Remove Lines wich Mean is less then ONE

    fprintf(1,'Remove %.0f Lines with invalid Characters: ',sum(~ok));

    jj = find(~ok);
    nj = prod(size(jj));

    if nj <= 6
       fprintf(1,nfrm,jj);
       fprintf(1,'\n');
    else
       fprintf(1,'\n');
       frm = nfrm(ones(nj,1),:);
       frm([(10:10:nj) nj],end) = char(10);
       frm = frm'; frm = frm(:)';
       fprintf(1,frm,jj);
    end

    ok = find(ok);

    i0 = i0(ok) + 1;
    lg =  n(ok) + 1;

    ok = grp2ind(i0,lg);

     d = d(ok);

end

%*******************************************
% Check Length of Line

i0 = cat( 1 , 0 , find( d == 10 ) );

 n = diff(i0) - 1;

 m = ( n - n0 );

ok = ( ( mod(m,4) == 2 ) & ( m >= 6 ) );
ok = ( ( mod(m,4) == 0 ) | ok );
ok = ( ( m >= 0 ) & ok );

if ~all(ok)

    if ~any(ok)
        msg = 'No Lines with valid Length.';
        return
    end

    %-----------------------------------

    fprintf(1,'Remove %.0f Lines with invalid Length: ',sum(~ok));

    jj = find(~ok);
    nj = prod(size(jj));

    if nj <= 6
       fprintf(1,nfrm,jj);
       fprintf(1,'\n');
    else
       fprintf(1,'\n');
       frm = nfrm(ones(nj,1),:);
       frm([(10:10:nj) nj],end) = char(10);
       frm = frm'; frm = frm(:)';
       fprintf(1,frm,jj);
    end

    ok = find(ok);

    i0 = i0(ok) + 1;
    lg =  n(ok) + 1;

    ok = grp2ind(i0,lg);

     d = d(ok);

    %-----------------------------------

    i0 = cat( 1 , 0 , find( d == 10 ) );

     n = diff(i0) - 1;

end

%-----------------------------------
% Check for Common Length

 m = n(1);

if ~all( n == m )

    [ni,m,nn] = grpmean(n,n);

    [nn,ii] = max(nn);

    m = n(ii);

    ok = ( n == m );

    fprintf(1,'Remove %.0f Lines with uncommon length of %.0f: ',sum(~ok),m);

    jj = find(~ok);
    nj = prod(size(jj));

    if nj <= 6
       fprintf(1,nfrm,jj);
       fprintf(1,'\n');
    else
       fprintf(1,'\n');
       frm = nfrm(ones(nj,1),:);
       frm([(10:10:nj) nj],end) = char(10);
       frm = frm'; frm = frm(:)';
       fprintf(1,frm,jj);
    end

    ok = find(ok);

    i0 = i0(ok) + 1;
    lg =  n(ok) + 1;

    ok = grp2ind(i0,lg);

    d = d(ok);

end

%*******************************************
% Reshape to Array

d(find( d == 10 )) = [];

d = reshape(d,m,size(d,1)/m);

d = permute(d,[2 1]);

n = size(d,2);

%-----------------------------------
% Make Decimal Digits from HEX

d = d - ( double('A') - double('9') - 1 ) * ( d >= double('A') );

d = d - double('0');

%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function fig = sc_plot(d,c,init,vok,sok,pref);

% Data = [ Day  T  C  P  S  [V0]  [V1]  [V2]  [V3]  [SBE38T] ]

base = datevec(d(1,1));
base = [ base(1)  01  01 ];

d(:,1) = d(:,1) - datenum(base);

nd = size(d,2);

ini = cell(0,3);  % { Column  Label  Tag }


i38 = []; if sok, i38 = nd; end
 
iv  = []; if vok, iv = 5 + ( 1 : vok ); end

ind = [ 4 2 i38 3 5 iv ];  % [ P T T38 C S Ext ]

ind = ind( find(~strcmp(init(ind,3),'')) );

n = size(ind,2);

fig = figure('units'            , 'pixels'   , ...
             'position'         , [ 100 35 200 200 ] , ...
             'paperorientation' , 'portrait' , ...
             'paperunits'       , 'inches'   , ...
             'color'            , [ 1  1  1 ] , ... 
             'menubar'          , 'none'  , ...
             'renderer'         , 'painters' , ...
             'visible'          , 'off'             );
 
papsi = get(fig,'papersize');

pappos = [0.8 1.2 papsi-2*[ 0.8  1.2 ] ];

set(fig,'paperposition',pappos)

wysiwyg
drawnow

%-----------------------------------------------------

xl = d([1 end],1)';

xl = xl + [ -1  1 ] * 0.01 * diff(xl);

xl = [ floor(xl(1))  ceil(xl(2)) ];

%-----------------------------------------------------
%  Determine AxesPosition

axe = zeros(n,1);

fs  = 10;  % FontSize

fh = 1/72 * fs / pappos(4);  % FontHeight normalized
vd = 2.0 * fh;      % vertical Distance
v0 = 4.0 * fh;      % Bottom
v1 = 3.0 * fh;      % Top

% Normalized AxeHeight

hh = ( 1 - v0 - v1 - (n-1)*vd ) / n .* ones(1,n);


for ii = 1 : n

    pos = [ 0.10  v0+sum(hh(ii+1:n))+(n-1-(ii-1))*vd  0.85  hh(ii) ];

    tag = upper(init{ind(ii),3});

    axe(ii) = axes('parent'     , fig    , ...
                   'position'   , pos    , ...
                   'xgrid'      , 'on'   , ...
                   'ygrid'      , 'on'   , ...
                   'xlim'       ,  xl    , ...
                   'layer'      , 'top'  , ...
                   'tickdir'    , 'in'   , ...
                   'box'        , 'on'   , ... 
                   'fontsize'   , fs     , ...
                   'fontweight' , 'bold' , ...
                   'tag'        , tag , ...
                   'nextplot'   , 'add'  , ...
                   'color'      , 'none'        );

    plot(d(:,1),d(:,ind(ii)),'.','tag',tag);

    lab = sprintf('%s   [%s]',init{ind(ii),[3 2]});

    set( get(axe(ii),'ylabel') , 'string' , lab , ...
         'fontsize' , fs , 'fontweight' , 'bold' );

    if strcmp(tag(1:4),'PRES')
       set(axe(ii),'ydir','reverse');
    end

end

is_time = ~isempty(which('timeaxis'));
is_zoom = ~isempty(which('axe_zoom'));

if is_time
   if is_zoom
      timeaxis(axe(n:-1:1),'x',10,base,'zoom');
   else
      timeaxis(axe(n:-1:1),'x',10,base);
   end
elseif is_zoom
   axe_zoom(axe,'new');
   axe_zoom(axe(n),'Help','off')
else
   set( get(axe(n),'xlabel') , 'string' , ...
        [ 'day of ' int2str(base(1)) ] , ...
        'fontsize'  , fs , ...
        'fontweight', 'bold', ...
        'interpreter' , 'none'   )
end


snr = [];
ii = strwcmp(upper(c(:,1)),{'TEMP*SN' 'COND*SN'});
if any(ii)
   ii = find(ii);
   for jj = ii(:)'
       [m,v] = str2vec(c{ii,2});
       if isnumeric(v) & ( prod(size(v)) == 1 )
          snr = v;
          break
       end
   end
end


ht = get( axe(1),'title');

tt = ['Seacat ',num2str(snr),'   Pref ',sprintf('%.1f',pref),' dbar'];

set( ht , 'string' , tt , 'tag' , 'TITLE' , ...
     'fontsize' , fs+1 , 'fontweight' , 'bold' );


%-----------------------------------------------------
% Determine correct FigurePosition in Screen


figpos = get(fig,'position');    % [ Left Bottom Width Height ] 
scr_si = get(0,'screensize');

FigMaxHeight = scr_si(4) - 100;

figpos(4) = figpos(4) + ( FigMaxHeight - figpos(4) ) * ...
                        ( FigMaxHeight < figpos(4) );

figpos(2) = scr_si(4) - figpos(4) - 60;

set(fig,'position',figpos)


%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ost=oxsat(s,t)
%OXSAT Oxygen saturation.
%  OST = OXSAT(S,T) return the oxygen saturation in ml/l as a function of
%  salinity s (PSS-78) and in situ temperature t (øC).
%
%  References:
%  Weiss, R. F., The solubility of nitrogen, oxygen and argon in water
%    and seawater. Deep-Sea Res., 17, 721-735, 1970.

%  Christian Mertens, IfM Kiel
%  $Revision: 1.0 $ $Date: 1996/01/05 13:20:22 $


a = [-173.4292 249.6339 143.3483 -21.8492];
b = [-0.033096 0.014259 -0.0017]; 

t = 0.01*(t + 273.15);
ost = a(1) + a(2)./t + a(3)*log(t) + a(4)*t + s.*(b(1) + (b(2) + b(3)*t).*t);
ost = exp(ost);

%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function X = sal78(P,T,CND,C1535,M)

% SAL78  Converts between Salinity and Conductivity
%
% SAL = SAL78( P , T68 , CND , C1535 , 0 )
% CND = SAL78( P , T68 , SAL , C1535 , 1 )
%
% FUNCTION SERVES FOR TWO PURPOSES (last Input):
%	1:	CONVERT CONDUCTIVITY TO SALINITY (0)
%	2:	CONVERT SALINITY TO CONDUCTIVITY (1)
%
%	ALGORITHMS RECOMMENDED BY JPOTS USING THE 1978
%	PRACTICAL SALINITY SCALE(IPSS-78) AND IPTS-68
%	FOR TEMPERATURE.
%
%       NOTE: The Conversion from IPTS-68 to ITS90 is:
%              T90 = 0.99976 * T68
%              T68 = 1.00024 * T90
%
% SAL78 COMPUTES EITHER 
%     THE CONDUCTIVITY RATIO    if C1535 = 1.0
%  or THE ABSOLUTE CONDUCTIVITY if C1535 = 42.914 = C(T=15,S=35)
%
%	UNITS:
%		PRESSURE        P         DBARS
%		TEMPERATURE     T         DEG.C.
%		SALINITY        S         NSU
%
%	RETURNS ZERO FOR CND < 0.0005  AND  last 0
%	RETURNS ZERO FOR SAL < 0.02    AND  last 1
%
% CHECKVALUES:
%
%           SAL78 = 1.888091
%       FOR   SAL =    40 NSU
%               T =    40 DEG
%               P =   10000 DBARS
%
%             CND = SAL78(10000,40,40,1,1) = 1.888091
%
%           SAL78 = 39.99999
%       FOR   CND = 1.888091
%               T =      40 DEG C.
%               P =   10000 DBARS
%
%             SAL = SAL78(10000,40,1.888091,1,0) = 39.999996
%


P = P/10 ;

if nargin == 3
	C1535 = 42.914 ;
	M = 0 ;
end

%ZERO SALINITY TRAP
if M == 0 
     zerocnd = find(CND < 5e-4) ; 
else
     zerocnd = find(CND < 0.02) ;
end

%SELECT BRANCH FOR SALINITY (M=0) OR CONDUCT.(M=1)
DT = T - 15.0 ;

if M == 0

     %CONVERT CONDUCTIVITY TO SALINITY
     R = CND/C1535 ;
     RT = R./(rt35(T).*(1.0 + c(P)./(b(T) + a(T).*R))) ;
     RT = sqrt(abs(RT)) ;
     %SALINITY RETURN
     X = sal(RT,DT) ;
     X(zerocnd) = zeros(size(zerocnd)) ;
else
     %CONVERT SALINITY TO CONDUCTIVITY
     %FIRST APPROXIMATION
     RT = sqrt(CND/35.0) ;
     SI = sal(RT,DT) ;
     [m,n] = size(CND) ;
     for i=1:m
          for j=1:n
               N = 0 ;
               DELS = 1 ;
               %ITERATE (MAX 10 ITERAT.) TO INVERT SAL POLYNOMIAL
               %FOR sqrt(RT)
               while (DELS > 1.e-4 & N < 10)
                    RT(i,j) = RT(i,j) + (CND(i,j) - SI(i,j))/dsal(RT(i,j),DT(i,j)) ;
                    SI(i,j) = sal(RT(i,j),DT(i,j)) ;
                    N = N + 1 ;
                    DELS = abs(SI(i,j) - CND(i,j)) ;
               end
          end
     end
     %COMPUTE CONDUCTIVITY RATIO
     RTT = rt35(T) .* RT .* RT ;
     AT = a(T) ;
     BT = b(T) ;
     CP = c(P) ;
     CP = RTT.*(CP + BT) ;
     BT = BT - RTT.*AT ; 
     R = sqrt(abs(BT.*BT+4.0*AT.*CP)) - BT ;
     %CONDUCTIVITY RETURN
     X = 0.5*C1535*R./AT ;
     X(zerocnd) = zeros(size(zerocnd)) ;
end


%**************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function D = dsal(XR,XT)
%DSAL	Required by SAL78.

%	10-Mar-93, C. Mertens, IfM, Kiel

	D=((((13.5405*XR-28.1044).*XR+42.2823).*XR+50.7702).*XR ...
                  -0.1692)+(XT./(1.0+0.0162*XT)).*((((-0.0720*XR+0.2544).*XR ...
                  -0.1125).*XR-0.0132).*XR-0.0056) ;

%**************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function rt = rt35(XT)
%RT35	Required by SAL78.

%	10-Mar-93, C. Mertens, IfM Kiel

	rt=(((1.0031E-9*XT-6.9698E-7).*XT+1.104259E-4).*XT ...
             +2.00564E-2).*XT+0.6766097 ;


%**************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function S = sal(XR,XT)
%SAL	Required by SAL78.

%	10-Mar-93, C. Mertens, IfM Kiel

	S = ((((2.7081*XR-7.0261).*XR+14.0941).*XR+25.3851).*XR ...
     	    -0.1692).*XR+0.0080+(XT./(1.0+0.0162*XT)).*(((((-0.0144*XR+ ...
     	     0.0636).*XR-0.0375).*XR-0.0066).*XR-0.0056).*XR+0.0005) ;



%**************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function x = a(XT)
%ACOEF	Required by SAL78.

%	10-Mar-93, C. Mertens, IfM Kiel

x=-3.107E-3*XT+0.4215 ;

%**************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function x = b(XT)
%BCOEF	Required by SAL78.

%	10-Mar-93, C. Mertens, IfM Kiel

x = (4.464E-4*XT+3.426E-2).*XT+1.0 ;

%**************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function x = c(XP)
%CCOEF	Required by SAL78.

%	10-Mar-93, C. Mertens, IfM Kiel

x = ((3.989E-12*XP-6.370E-8).*XP+2.070E-4).*XP ;
