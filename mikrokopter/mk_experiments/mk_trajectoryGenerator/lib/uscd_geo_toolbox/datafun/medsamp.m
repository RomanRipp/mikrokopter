function [msg,dat,day,org] = medsamp(dat,ini,smp,lim);

% MEDSAMP  Lowsampling of TimeSeries using median filter
%
% [Msg,Dat,Day] = MEDSAMP(Data,Init,Sample,Limits);
%
% Data  Data, one TimeSeries per Column, equal timestep by Init
% 
% Init  Initialisation for Time of Data
%       [ StartDay[day]  Intervall[sec] ]
%
% Sample  Parameter for LowSampling
%         [ LowSampleintervall[min]  [DayStart]  [DayEnd] ]
%
%         Unknown or NaN for [DayStart]  [DayEnd] 
%           use from Start- and EndTime of TimeSeries
%
% Limits  Optional Limits for Data  [ Min Max ]
%         Values Outside are set to NaN
%
%--------------------------------------------------------------
% Return the TimeVector only
%
% [Msg,Day] = MEDSAMP(Sample);
%
% [Msg,Day,Interv,OrgDay] = MEDSAMP(NData,Init,Sample);
%
%--------------------------------------------------------------
% 
% MEDSAMP(ASC_File,Init,Sample,lim) load the Data from ASCII-File, using READ_ASC
%
% MEDSAMP(CSA_File,[],Sample,lim)  load the Data from CSA-File, using READ_CSA
%               requested Header in CSA-Files: start_time, samp_interv, npts
%
%--------------------------------------------------------------
%

Nin  = nargin;
Nout = nargout;

msg = [];
day = [];
org = [];

is_org = ( Nout == 4 );

get_time = ( Nin == 1 );

if get_time
   smp = dat;
   dat = [];
   ini = NaN * ones(1,2);
end

if ~get_time & ( Nin < 3 )
   msg = 'Not enough Input Arguments.';
   return
end

if Nin < 4
   tlm = [];
end

if Nin < 5
   lim = [];
end

%************************************************************
% Check Inputs

msg = cell(0,1);

%------------------------------------------------------------
% Check Data

if chkstr(dat,1)

   file = dat;

   m = '';

   if isempty(ini)

      ini = { 'start_time'  [ 3  6 ]
              'samp_interv'   1      
              'npts'          1      };
      val = ini;
      val(:,2) = { {NaN} };
      val = permute(val,[2 1]);
      val = struct(val{:});
      try
         [m,head,cols,dat] = read_csa(file);
      catch
          m = sprintf('Error call READ_CSA.\n%s',lasterr);
      end

      if isempty(m)
         [ok,val] = gethead(head,val,ini);
         if ~ok
             m = sprintf('%s, ',ini{:,1});
             m = sprintf('Requested HeaderVariables: %s',m(1:(end-2)));
         elseif ~( val.npts == size(dat,1) )
             m = 'Invalid Size of Data.';
         else
            ini = [ val.start_time  val.samp_interv ];
         end
      end

   else

      try
         [m,head,dat] = read_asc(file);
      catch
          m = sprintf('Error call READ_ASC.\n%s',lasterr);
      end

   end

   if ~isempty(m)
       msg = cat(1,msg,{sprintf('Invalid DataFile "%s".\n%s',file,m)});
       dat = [];
   end

elseif ~isempty(dat)

   if ~isnumeric(dat) & ( ndims(dat) == 2 )
       msg = cat(1,msg,{'Data must be a 2D-Numeric.'});
   end

   get_time = ( prod(size(dat)) == 1 );

   if get_time
      get_time = ( ( mod(dat,1) == 0 ) & ( dat > 0 ) );
   end

end

%------------------------------------------------------------
% Check Init

if ~( isnumeric(ini) & ( prod(size(ini)) == 2 ) )
    msg = cat(1,msg,{'Init must be a 2-Element numeric.'});
end

%------------------------------------------------------------
% Check Sample

p = prod(size(smp));

if ~( isnumeric(smp) & any( p == [ 1  2  3 ] ) )
    msg = cat(1,msg,{'Sample must be a 1-3-Element numeric.'});
else
    smp = smp(:)';
    smp = cat( 2 , smp(:)' , NaN*ones(1,3-p) );
    tlm = smp([2 3]);
    smp = smp(1);
end

%------------------------------------------------------------
% Check Limits

if ~( isnumeric(smp) & any( prod(size(lim)) == [ 0  2 ] ) )
    msg = cat(1,msg,{'Limit must be empty or a 2-Element numeric.'});
elseif ~isempty(lim)
    if diff(lim) <= 0
       warning('Decreasing Limits.')
    end
end

%------------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    msg = sprintf('Invalid Inputs.\n%s',msg);
    return
end

msg = '';

%************************************************************
% Get Time

is_org = ( is_org | get_time );

s = size(dat); p = prod(s);

is_flip = ( isequal(s,[1 p]) & ( p > 1 ) );

if is_flip
   dat = permute(dat,[2 1]);
end

%------------------------------------------------------------
% Parameter

nd = size(dat,1);
nc = size(dat,2);

if get_time & ( nd*nc == 1 )
   nd = dat;
   nc = 1;
end

dlm = ini([1 1]);  % day              % StartEndDay
itv = ini(2);      % sec              % Intervall

if isnan(tlm(1)), tlm(1) = dlm(1); end

base = floor(tlm(1));                 % BaseDay

tlm = tlm - base;
dlm = dlm - base;

tlm = tlm * 24 * 3600;                % day ---> sec
smp = smp * 60;                       % min ---> sec

dlm = dlm * 24 * 3600;                % day ---> sec

dlm(2) = dlm(2) + ( nd - 1 ) * itv;   % Original EndTime

%------------------------------------------------------------
% Set NaN-Values of requested TimeLimits 
%                to  original TimeLimits

if isnan(tlm(2)), tlm(2) = dlm(2); end

if smp < 3*itv
   warning('Intervall to lowsample less then 3 times original.');
   smp = 3 * itv;
end

% Start/End in range of new Samples

op = [ -1  1 ];

vlm = op .* floor( (op.*dlm)/smp );
wlm = op .* floor( (op.*tlm)/smp );

% MedianWindow

ni = 2 * floor( smp/itv / 2 ) + 1;    % Odd Number of orig samples in new Intervall
ni = ni + 2 * ( ni*itv < smp );       % Raise by 2 if new Intervall not covered

%------------------------------------------------------------------
% Check Time

if any(isnan(tlm))
   nt = 0;
else
   nt  = diff(wlm) + 1;
   nt  = max(0,nt);
end

if nt == 0
   day = zeros(0,1);
else
   day = ( wlm(1) : wlm(2) )' * smp  / 24 / 3600 + base;
end

if is_org
   if nd > 0
      org = ( dlm(1) + ( ( 1 : nd )' - 1 ) * itv ) / 24 / 3600 + base;
   else
      org = zeros(0,nc);
   end
end

if get_time
   dat = day;
   day = ni;
   if is_flip
      dat = permute(dat,[2 1]);
   end
   return
end

if isempty(dat)
   if isempty(wlm)
      dat = [];
   elseif isempty(ini)
      dat = day;
   else
      dat = NaN * ones(nt,size(dat,2));
   end
   if is_flip
      dat = permute(dat,[2 1]);
      day = permute(day,[2 1]);
      org = permute(org,[2 1]);
   end
   return
end

%************************************************************
% Check for Limits and NaN

%--------------------------------------------------
% Check Data with Limits

if ~isempty(lim)
    lm = mean(lim);
    dl = diff(lim) / 2;
    dat( find( abs(dat-lm) > dl ) ) = NaN;
end

%--------------------------------------------------
% Check Data for NaN's

ii  = isnan(dat);

isn = any(ii,1);
aln = all(ii,1);

if all(aln)
   dat = NaN * zeros(nt,nc);
   return
elseif any(isn)
   ii = find(ii);
   dat(ii) = NaN;
end

%************************************************************

if is_org
   org = cat( 2 , dat , org );
end

n2 = ( ni - 1 ) / 2;                           % Half Intervall without center

dd = ( vlm * smp  - dlm ) / itv;               % Nr of overlapping orig samples at begin and End

vlm(1) = vlm(1) + ( dd(1) <  n2 );            % Start a new sample later 
vlm(2) = vlm(2) - ( dd(2) > -n2 );            % End   a new sample before

vlm(1) = max(vlm(1),wlm(1));                   % Check with TimeStart
vlm(2) = min(vlm(2),wlm(2));                   % Check with TimeEnd

if vlm(1) > vlm(2)
   dat = NaN * zeros(nt,nc);
   if is_flip
      dat = permute(dat,[2 1]);
      day = permute(day,[2 1]);
      org = permute(org,[2 1]);
   end
   return
end

%--------------------------------------------------

nn = diff(vlm) + 1;

if nd < nt
   dat = cat( 1 , dat , NaN * ones(nt-nd,nc) );
    nd = nt;
end

% CenterIndex in orig samples

ii = round( ( (vlm(1):vlm(2))*smp - dlm(1) ) / itv ) + 1;
ii = ii(:);

i0 = vlm(1) - wlm(1);
jj = max( 0 , i0 ) + ( 1 : nn );

%%% fprintf(1,'%3.0f ',n2,dd,ii([1 end]),day(jj([1 nn]))*24*3600);

n2 = ( -n2 : n2 );

for ic = find(~aln)

    kk = ( ic - 1 ) * nd;

    if isn(ic)
       dat(jj,ic) = mednan(dat(ii(:,ones(ni,1))+n2(ones(nn,1),:)+kk),2);
    else
       dat(jj,ic) = median(dat(ii(:,ones(ni,1))+n2(ones(nn,1),:)+kk),2);
    end

end

if i0 > 0
   dat(1:i0,:) = NaN;
end

if nt < nd
   dat = dat(1:nt,:);
end

if is_flip
   dat = permute(dat,[2 1]);
   day = permute(day,[2 1]);
   org = permute(org,[2 1]);
end

%-------------------------------------------------------------
% Check correct Time-relation
%-------------------------------------------------------------
if 0
%-------------------------------------------------------------

tt = ( dlm(1) + ( ( 1 : nd )' - 1 ) * itv ) / 24 / 3600 + base;
tt = median(tt(ii(:,ones(ni,1))+n2(ones(nn,1),:)),2);

disp('Median Time'), disp(datestr(tt(1:3),0))
disp(' ')
disp('Output Time'), disp(datestr(day(i0+(1:3)),0))

%-------------------------------------------------------------
end
%-------------------------------------------------------------


%*************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,val] = gethead(head,val,ini)

% GETHEAD  Check, Transform Header from CSA-Files

n = size(ini,1);

ok = zeros(n,1);

for ii = 1 : n
    [ok(ii),v,uni] = headval(head,ini{ii,[1 2]});
    if isequal(ok(ii),1)
        switch ini{ii,1}
          case 'npts'
               val.npts = v;
          case 'start_time'
               v = cat( 2 , v , zeros(1,6-size(v,2)) );
               val.start_time = datenum(v(1),v(2),v(3),v(4),v(5),v(6));
          case 'samp_interv'
               scl = [ 1 60 60 24 ];
               kk  = 1 * strwcmp(uni,'*sec*') + ...
                     2 * strwcmp(uni,'*min*') + ...
                     3 * strwcmp(uni,'*hour*') + ...
                     4 * strwcmp(uni,'*day*');
               if ( kk == 0 )
                   fprintf(1,'\nWarning: Can''t get IntervallUnits, assume seconds\n');
               end
               kk = max(1,kk); 
               val.samp_interv = v * prod(scl(1:kk));  % Seconds !!!
        end
    else
        if isnan(ok(ii))
           m = 'Invalid';
        elseif ok(ii) == -1
           m = 'Multiple';
        else
           m = 'No';
        end
        fprintf(1,'\n%s %s.',m,ini{ii,1});
        break
    end
end

ok = all( ok == 1 );

%*************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,val,uni] = headval(head,prop,len)

val = [];
uni = '';

ok = isequal(len,-1);
if ok
   return
end

jj = strcmp(lower(head(:,1)),lower(prop));

ok = any(jj) - 2 * ( sum(jj) > 1 );

if ~( ok == 1 )
    return
end

jj = find(jj);

uni = head{jj,2};
val = head{jj,3};

if ~( isnumeric(head{jj,3}) & any(size(head{jj,3},2) == len ) )
     ok = NaN;
end

%*************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = chkstr(str,opt)


% CHKSTR  Checks Input for String
%
%  ok = chkstr(str,Option)
%
%  Option ~= 0 ==> only true for nonempty Strings
%
%   default:   Option == 0  (true for empty Strings)
%

 
if nargin < 2
   opt = 0;
end

ok = ( strcmp( class(str) , 'char' )      & ...
       ( prod(size(str)) == size(str,2) ) & ...
       ( isequal(opt,0) | ~isempty(str) )         );

%*************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function y = mednan(x,dim)

% MEDNAN   Median value, ignores NaN
%
%   For vectors, MEDNAN(X) is the median value of the elements in X.
%   For matrices, MEDNAN(X) is a row vector containing the median
%   value of each column.  For N-D arrays, MEDNAN(X) is the median
%   value of the elements along the first non-singleton dimension
%   of X.
%
%   MEDNAN(X,DIM) takes the median along the dimension DIM of X.
%
%   Example: If X = [0 1 2
%                    3 4 5]
%
%   then mednan(X,1) is [1.5 2.5 3.5] and median(X,2) is [1
%                                                         4]
%
%   See also MEDIAN, MEAN, STD, MIN, MAX, COV.

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.11 $  $Date: 1997/11/21 23:23:56 $

if nargin==1, 
  dim = min(find(size(x)~=1)); 
  if isempty(dim), dim = 1; end
end
if isempty(x), y = []; return, end

siz = [size(x) ones(1,dim-ndims(x))];
n = size(x,dim);

% Permute and reshape so that DIM becomes the row dimension of a 2-D array
perm = [dim:max(length(size(x)),dim) 1:dim-1];
x = reshape(permute(x,perm),n,prod(siz)/n);

% Sort along first dimension
x = sort(x,1);

y  = NaN * ones(1,size(x,2));

%-----------------------------------------------------------
% Check for NaN's

ii = any(isnan(x),1);

jj = find(~ii);  % NoNan
ii = find( ii);  % IsNan

%-----------------------------------------------------------
% Median of NoNaN - Columns

if ~isempty(jj)
    if rem(n,2) % Odd number of elements along DIM
       y(jj) = x((n+1)/2,jj);
    else                % Even number of elements along DIM
       y(jj) = (x(n/2,jj) + x((n/2)+1,jj))/2;
    end
end

%-----------------------------------------------------------
% Median of IsNaN - Columns

if ~isempty(ii)

    % Run over Columns
    for jj = ii

        kk = ~isnan(x(:,jj));
        n  = sum(kk,1);
        r  = rem(n,2);  
        r  = r - 1 * ( n == 0 );
   
        kk = find(kk);

       if     r == 0  % Odd number of elements along DIM

              y(jj) = (x(kk(n/2),jj) + x(kk(n/2+1),jj))/2;

       elseif r == 1  % Even number of elements along DIM

              y(jj) = x(kk((n+1)/2),jj);     

       end

    end

end

%-----------------------------------------------------------

siz(dim) = 1;

y = ipermute(reshape(y,siz(perm)),perm);
