function fig = stpl_adcp(file,lev,varargin)

% STPL_ADCP  StickPlot of RDI-ADCP-Data, proceeded by RDI2MAT
%
% Fig = STPL_ADCP(File,Level,Title,Smooth,Unit,Maximum)
%
%-----------------------------------------------------------------------------
% Basic Inputs:
%
% File    FileName of RODB-MAT-File for ADCP-Data, see RDI2MAT
%
% Level   BinLevels to show (IndexVector)
%          EMPTY ==> all Bins will shown
%
%-----------------------------------------------------------------------------
% Optional Inputs (accepted in any Order):
%
% Titel   Titel of Grafic, default: 'Moor|Stat, Type #SNR, Depth, Smooth'
%
% Smooth  SmoothWindow [hour] for Data, using MEANIND1
%          the smoothed Data will lowsampled!!!
%
%          i*Window + Real  |  -Window
%
%          Imaginary Window:  A POSitive Real part will use a  
%                               LEFTside window only. 
%                             (looking backward, aligned at begin to origin)
%
%                              A NEGative Real part will use a 
%                               RIGHTside window only.
%                             (looking foreward, aligned at end to origin)
%
%                              If the Real part is ZERO, 
%                               the values at Begin and End 
%                               will aligned to the origin, 
%                               a negative Value of Window prevents this.
%
%           Negative Window    equal to -i * |Window|, see above
%
%
% Unit, Maximum    positive single numerics, see STPL
%
%-----------------------------------------------------------------------------
% 
% see also:  STPL, MEANIND1, RDI2MAT
%


Nin = nargin;

%*************************************************************************
% Check File and Level

if Nin < 1
   error('Input File is missing.')
end

if ~( ischar(file) & ~isempty(file)  & ...
      ( prod(size(file)) == size(file,2) ) )
    error('File must be a String.');
end

if ~( exist(file,'file') == 2 )
    error(sprintf('File "%s" doesn''t exist.',file))
end

if Nin < 2
   lev = [];
end

%*************************************************************************
% Check Other: Titel, Unit, MaxInt, Smooth

vin = varargin;

titel = -1;
smt   = 0;

if ~isempty(vin)

    n = prod(size(vin));
 
    ok = zeros(n,1);

    for ii = 1 : n
        v = vin{ii};
        if ( ischar(v) & ( prod(size(v)) == size(v,2) ) )
           if ~isempty(v)
               titel = v;
           end
        elseif ~( isnumeric(v) & ( prod(size(v)) == 1 ) )
           error(sprintf('Invalid %.0f Input.',ii+2));
        elseif ~( imag(v) == 0 )
           smt = imag(v) + i*real(v);
        elseif  v < 0
           smt = v;
        else
           ok(ii) = 1;
        end
    end

    vin = vin(find(ok));

    if prod(size(vin)) > 2
       error('Too many InputArguments for STPL.')
    end

end


%*************************************************************************
% Data

load(file,'dat');

n = size(dat.z,1);

if ~isempty(lev)
    if ~all( ( 1 <= lev ) & ( lev <= size(dat.u,1) ) )
        error(sprintf('Level exeeds Matrix Dimensions: %.0f',size(dat.u,1)))
    end
else
    lev = ( 1 : n );
end

dat.u = dat.u(lev,:);
dat.v = dat.v(lev,:);
dat.z = dat.z(lev);

%*************************************************************************
% Smooth

if ~isequal(smt,0)

    dt = median(diff(dat.day)) * 24;

    rs = real(smt) / dt;    % hour --> samples
    rs = sign(rs) * ( 2 * floor( abs(rs) / 2 ) + 1 );

    int = rs + i*imag(smt);

    rs  = rs * dt;    % samples --> hour

    dat.u = meanind1(dat.u',int,'cos')';
    dat.v = meanind1(dat.v',int,'cos')';

    % LowSample

    smp = floor( abs(real(abs(int))) / 3 );
    ind   = ( 1 : smp : size(dat.day,2) );

    dat.u = dat.u(:,ind);
    dat.v = dat.v(:,ind);

    dat.day = dat.day(:,ind);

end

%*************************************************************************
% Make CellArrays, Sort by Depth

dat.u = num2cell(dat.u,2);
dat.v = num2cell(dat.v,2);

[dat.z,si] = sort(dat.z);

dat.u = dat.u(si);
dat.v = dat.v(si);

n = prod(size(dat.z));

%*************************************************************************
% Time

if ~isequal(dat.base([2 3]),[1 1])
    dat.day = dat.day + datenum(dat.base(1),dat.base(2),dat.base(3));
    dat.base([2 3]) = 1;
    dat.day = dat.day - datenum(dat.base(1),dat.base(2),dat.base(3));
end

dat.day = { dat.day };
dat.day = dat.day(ones(1,n));

%************************************************************************
% Label

for ii = 1 : n
    lab{ii} = sprintf('%4.0fm',dat.z(ii));
end

%************************************************************************
% Title

if isequal(titel,-1)

   titel = dat.typ;

   if isempty(titel)
      titel = dat.sys;
   end

   if isempty(titel)
      titel = 'ADCP';
   end

   if ~isempty(dat.moor)
       titel = sprintf('%s, %s',dat.moor,titel);
   elseif ~isempty(dat.stat)
       titel = sprintf('%s,  %s',dat.stat,titel);
   end

   if ~isempty(dat.snr)
       titel = sprintf('%s #%.0f',titel,dat.snr);
   end

   if ~isempty(dat.dpt)
       titel = sprintf('%s,  %.0fm',titel,dat.dpt);
   end

   if ~isequal(smt,0)
       titel = sprintf('%s,  %.0fh lowpass',titel,abs(rs));
   end

end

%*************************************************************************

[fig,axe] = stpl(dat.day,dat.u,dat.v,dat.base(1),vin{:},lab);

ht = get(axe,'title');

set(ht,'string',titel, ...
        'fontweight','bold','interpreter','none','fontsize',12)

