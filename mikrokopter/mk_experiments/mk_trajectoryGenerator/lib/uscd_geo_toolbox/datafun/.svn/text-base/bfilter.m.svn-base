function [varargout] = bfilter(varargin)

% BFILTER  Butterworth digital and analog filter (BUTTER)
%
% [B,A] = BFILTER( N , Wn , Mode , Option )
%
% Returns the Filter coefficients for an Nth order Butterworth filter
%   in length N+1 vectors B (numerator) and A (denominator). 
%   The coefficients are listed in descending powers of z. 
%   The cutoff or band frequency Wn must be 0.0 < Wn < 1.0, 
%   with 1.0 corresponding to half the sample rate.
%
%   Mode = 'lp' | 'low'  | 'lowpass'   lowpass   default
%          'hp' | 'high' | 'highpass'  highpass 
%          'bp' | 'band' | 'bandpass'  bandpass; Wn = [W1 W2]
%          'bs' | 'stop' | 'bandstop'  bandstop; Wn = [W1 W2]
%
%   Option = 'z'  design digital Filter, default, ( 0 < Wn < 1 )
%   Option = 's'  design analog Butterworth filters. In this case, 
%                  Wn is in [rad/s]  and it can be greater than 1.
%
%   When used with three left-hand arguments, as in
%   [Z,P,K] = BFILTER(N,Wn,...), the zeros and poles are returned in
%   length N column vectors Z and P, and the gain in scalar K. 
%
%   When used with four left-hand arguments, as in
%   [A,B,C,D] = BFILTER(N,Wn,...), state-space matrices are returned.
%
%
% [ Y , ... ] = BFILTER( X , N , Wn , ... )
%
% Filters the data in vector X, using a Zero-phase forward 
% and reverse digital filter (FILTFILT). The 2-directional filter
% will applied from both sides of the vector to avoid an overshooting
% at the begin.
%
% forward digital filter:
%
%     y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                      - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
%
%  The length of the input X must be more than three times
%   the filter order N.
%
%   After filtering in the forward direction, the filtered sequence is then 
%   reversed and run back through the filter; Y is the time reverse of the 
%   output of the second filtering operation.  The result has precisely zero 
%   phase distortion and magnitude modified by the square of the filter's 
%   magnitude response.  Care is taken to minimize startup and ending 
%   transients by matching initial conditions.
%
% see also: MFILTER, MEANIND1, MEANIND2, NOISE
%           BUTTER, FILTFILT (SignalProcessingToolbox)
%


%---------------------------------------------------------------------
% SubFunctions from Signal and Matlab-Toolboxes
%---------------------------------------------------------------------
%
% FILTFILT (FILT2) Zero-phase forward and reverse digital filtering.
% y = filtfilt(b,a,x)
%
% BUTTER (BUTTCOEF) Butterworth digital and analog filter design.
%
%   [Msg,B,A]     = BUTTCOEF(N,Wn,Mode,Option)
%   [Msg,Z,P,K]   = BUTTCOEF(N,Wn,Mode,Option)
%   [Msg,A,B,C,D] = BUTTCOEF(N,Wn,Mode,Option)
%
% IIRCHK  Parameter checking for BUTTER, CHEBY1, CHEBY2, and ELLIP.
% [btype,analog,errStr] = iirchk(Wn,varargin)
%
% BUTTAP Butterworth analog lowpass filter prototype.
% [z,p,k] = buttap(n)
%
% BILINEAR Bilinear transformation with optional frequency prewarping.
% [zd, pd, kd, dd] = bilinear(z, p, k, fs, fp, fp1)
%
% LP_TRANS Lowpass to low/band/high/stop - pass Transformation
% [at,bt,ct,dt] = lp_transf(mode,a,b,c,d,wo,bw)
%
% LP2LP Lowpass to lowpass analog filter transformation.
% LP2BP Lowpass to bandpass analog filter transformation.
% LP2HP Lowpass to highpass analog filter transformation.
% LP2BS Lowpass to bandstop analog filter transformation.
%
% ABCDCHK Checks dimensional consistency of A,B,C,D matrices.
% [msg,A,B,C,D] = abcdchk(A,B,C,D)
%
% POLY Convert roots to polynomial.
% c = poly(x)
%
% TZERO  Transmission zeros of LTI systems.
% [z,gain] = tzero(a,b,c,d)
%
% SS2TF  State-space to transfer function conversion.
% [num, den] = ss2tf(a,b,c,d,iu)
%
% SS2ZP  State-space to zero-pole conversion.
% [z,p,k] = ss2zp(a,b,c,d,iu)
%
% TF2SS  Transfer function to state-space conversion.
% [a,b,c,d] = tf2ss(num, den)
%
% ZP2SS  Zero-pole to state-space conversion.
% [a,b,c,d] = zp2ss(z,p,k)
%
% CPLXPAIR Sort numbers into complex conjugate pairs.
% z = cplxpair(x,tol,dim)
%
% SHIFTDIM Shift dimensions.
% [b,nshifts] = shiftdim(x,n)
%
%---------------------------------------------------------------------

Nin  = nargin;
Nout = nargout;

varargout = cell(1,Nout);

if Nin <1
   exit
end

l = prod(size(varargin{1}));
 
v = ( l > 1 );  % True for VectorInput

if Nin < 2+v
   error('Not enough InputArguments.')
end

[msg,varargout{(1:min(Nout-v,4))+v}] = buttcoef(varargin{(1+v):end});

if ~isempty(msg)
    error(msg)
elseif ~v
    return
end

if Nout == 3
   [b,a] = varargout{[2 3]}
else
   [msg,b,a] = buttcoef(varargin{(1+v):end});
end

si = size(varargin{1});

nd = size(si,2);

prm = [];
if ( si(1) == 1 )
   prm = cat(2,[2 1],(3:nd));
   varargin{1} = permute(varargin{1},prm);
   si = si(prm);
end

n = si(1);

rsh = ( nd > 2 );
if rsh
   varargin{1} = reshape(varargin{1},[ n  prod(si(2:nd)) ]);
end

int = (n+1)/2 + [ -1  1 ] * 1/varargin{3};

w = ( ( 1 : n )' - int(2) ) / ( int(2) - int(1) ); 
w = min(max(w,-1),0);
w = ( 1 + cos( pi * w ) ) / 2;

%%% w = winint((1:n)',int,'cos',1);

o = ones(1,size(varargin{1},2));

% Filter twice backward

varargout{1} = filt2(b,a,varargin{1}(n:-1:1,:));

varargout{1} = varargout{1}(n:-1:1,:);

% Add foreward twice filter

varargout{1} =    w(:,o)  .* filt2(b,a,varargin{1}) + ...
               (1-w(:,o)) .* varargout{1};

if rsh
   varargout{1} = reshape(varargout{1},si);
end
   
if ~isempty(prm)
   varargin{1} = ipermute(varargin{1},prm);
end
    

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function y = filt2(b,a,x)
%FILT2 Zero-phase forward and reverse digital filtering.
%   Y = FILT2(B, A, X) filters the data in vector X with the filter described
%   by vectors A and B to create the filtered data Y.  The filter is described 
%   by the difference equation:
%
%     y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                      - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
%
%
%   After filtering in the forward direction, the filtered sequence is then 
%   reversed and run back through the filter; Y is the time reverse of the 
%   output of the second filtering operation.  The result has precisely zero 
%   phase distortion and magnitude modified by the square of the filter's 
%   magnitude response.  Care is taken to minimize startup and ending 
%   transients by matching initial conditions.
%
%   The length of the input x must be more than three times
%   the filter order, defined as max(length(b)-1,length(a)-1).
%
%   Note that FILT2 should not be used with differentiator and Hilbert FIR
%   filters, since the operation of these filters depends heavily on their
%   phase response.
%
%   See also FILTER.

%   References: 
%     [1] Sanjit K. Mitra, Digital Signal Processing, 2nd ed., McGraw-Hill, 2001
%     [2] Fredrik Gustafsson, Determining the initial states in forward-backward 
%         filtering, IEEE Transactions on Signal Processing, pp. 988--992, April 1996, 
%         Volume 44, Issue 4

%   Author(s): L. Shure, 5-17-88
%   revised by T. Krauss, 1-21-94
%   Initial Conditions: Fredrik Gustafsson
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.7.4.2 $  $Date: 2004/12/26 22:15:49 $

    error(nargchk(3,3,nargin))
    if (isempty(b) || isempty(a) || isempty(x))
        y = [];
        return
    end

    [m,n] = size(x);
    if (n>1) && (m>1)
        y = x;
        for i=1:n  % loop over columns
           y(:,i) = filt2(b,a,x(:,i));
        end
        return
        % error('Only works for vector input.')
    end
    if m==1
        x = x(:);   % convert row to column
    end
    len = size(x,1);   % length of input
    b = b(:).';
    a = a(:).';
    nb = length(b);
    na = length(a);
    nfilt = max(nb,na);

    nfact = 3*(nfilt-1);  % length of edge transients

    if (len<=nfact),    % input data too short!
        error('Data must have length more than 3 times filter order.');
    end

% set up filter's initial conditions to remove dc offset problems at the 
% beginning and end of the sequence
    if nb < nfilt, b(nfilt)=0; end   % zero-pad if necessary
    if na < nfilt, a(nfilt)=0; end
% use sparse matrix to solve system of linear equations for initial conditions
% zi are the steady-state states of the filter b(z)/a(z) in the state-space 
% implementation of the 'filter' command.
    rows = [1:nfilt-1  2:nfilt-1  1:nfilt-2];
    cols = [ones(1,nfilt-1) 2:nfilt-1  2:nfilt-1];
    data = [1+a(2) a(3:nfilt) ones(1,nfilt-2)  -ones(1,nfilt-2)];
    sp = sparse(rows,cols,data);
    zi = sp \ ( b(2:nfilt).' - a(2:nfilt).'*b(1) );
% non-sparse:
% zi = ( eye(nfilt-1) - [-a(2:nfilt).' [eye(nfilt-2); zeros(1,nfilt-2)]] ) \ ...
%      ( b(2:nfilt).' - a(2:nfilt).'*b(1) );

% Extrapolate beginning and end of data sequence using a "reflection
% method".  Slopes of original and extrapolated sequences match at
% the end points.
% This reduces end effects.
    y = [2*x(1)-x((nfact+1):-1:2);x;2*x(len)-x((len-1):-1:len-nfact)];

% filter, reverse data, filter again, and reverse data again
    y = filter(b,a,y,zi*y(1));
    y = y(length(y):-1:1);
    y = filter(b,a,y,zi*y(1));
    y = y(length(y):-1:1);

% remove extrapolated pieces of y
    y([1:nfact len+nfact+(1:nfact)]) = [];

    if m == 1
        y = y.';   % convert back to row if necessary
    end


%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,a,b,c,d] = buttcoef(n,Wn,varargin)

% BUTTCOEF Butterworth digital and analog filter design.
%
%   [B,A] = BUTTCOEF(N,Wn) designs an Nth order lowpass digital
%   Butterworth filter and returns the filter coefficients in length 
%   N+1 vectors B (numerator) and A (denominator). The coefficients 
%   are listed in descending powers of z. The cutoff frequency 
%   Wn must be 0.0 < Wn < 1.0, with 1.0 corresponding to 
%   half the sample rate.
%
%   If Wn is a two-element vector, Wn = [W1 W2], BUTTER returns an 
%   order 2N bandpass filter with passband  W1 < W < W2.
%
%   [B,A] = BUTTCOEF(N,Wn,Mode,Option) 
%
%   Mode = 'lp' | 'low'  | 'lowpass'   lowpass   default
%          'hp' | 'high' | 'highpass'  highpass 
%          'bp' | 'band' | 'bandpass'  bandpass; Wn = [W1 W2]
%          'bs' | 'stop' | 'bandstop'  bandstop; Wn = [W1 W2]
%
%   Option = 'z'  design digital Filter, default, ( 0 < Wn < 1 )
%   Option = 's'  design analog Butterworth filters.  
%                  In this case, Wn is in [rad/s]  and it can be greater than 1.0.
%
%   When used with three left-hand arguments, as in
%   [Z,P,K] = BUTTCOEF(...), the zeros and poles are returned in
%   length N column vectors Z and P, and the gain in scalar K. 
%
%   When used with four left-hand arguments, as in
%   [A,B,C,D] = BUTTCOEF(...), state-space matrices are returned.
%
%   See also BUTTORD, BESSELF, CHEBY1, CHEBY2, ELLIP, FREQZ, FILTER.

%   Author(s): J.N. Little, 1-14-87
%   	   J.N. Little, 1-14-88, revised
%   	   L. Shure, 4-29-88, revised
%   	   T. Krauss, 3-24-93, revised
%   (C) 1988-2004 The MathWorks, Inc.
%   $Revision: 1.8.4.3 $  $Date: 2004/10/18 21:07:33 $

%   References:
%     [1] T. W. Parks and C. S. Burrus, Digital Filter Design,
%         John Wiley & Sons, 1987, chapter 7, section 7.3.3.

Nout = nargout;

a = [];
b = [];
c = [];
d = [];

msg = cell(0,1);

ok = ( isnumeric(n) & ( prod(size(n)) == 1 ) );
if ok
   ok = ( ( n >= 0 ) & ( mod(n,1) == 0 ) );
end
if ~ok
    msg = cat(1,msg,{'Filter order must be a positive Integer.'});
elseif n > 500
    msg = cat(1,msg,{'Filter order too large.'})
end

[btype,analog,errStr] = iirchk(Wn,varargin{:});

if ~isempty(errStr)
    msg = cat(1,msg,{errStr});
end

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    return
end

msg = '';

%------------------------------------------------------------------

% step 1: get analog, pre-warped frequencies
if ~analog,
    fs = 2;
    Wn = 2*fs*tan(pi*Wn/fs);
end

% step 2: convert to low-pass prototype estimate

Bw = [];

if any(strcmp(btype,{'bp' 'bs'}))
   Bw = Wn(2) - Wn(1);
   Wn = sqrt(Wn(1)*Wn(2)); % center frequency
end


% step 3: Get N-th order Butterworth analog lowpass prototype
[z,p,k] = buttap(n);

% Transform to state-space
[a,b,c,d] = zp2ss(z,p,k);

% step 4: Transform to lowpass, bandpass, highpass, or bandstop of desired Wn
[a,b,c,d] = lp_transf(btype,a,b,c,d,Wn,Bw);


% step 5: Use Bilinear transformation to find discrete equivalent:
if ~analog,
    [a,b,c,d] = bilinear(a,b,c,d,fs);
end

if     Nout == 4     % zero-pole-gain [Msg,Z,P,K]

       [z,b,c] = ss2zp(a,b,c,d,1);
        a      = buttzeros(btype,n,Wn,analog);

elseif Nout <= 3     % polynomial     [Msg,B,A]

        b = poly(a);
        a = buttnum(btype,n,Wn,Bw,analog,b);
        % num = poly(a-b*c)+(d-1)*den;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = buttnum(btype,n,Wn,Bw,analog,den)

% This internal function returns more exact numerator vectors
% for the num/den case.
% Wn input is two element band edge vector

if analog

  switch btype
    case 'lp'  % lowpass
        b = [zeros(1,n) n^(-n)];
        b = real( b*polyval(den,-j*0)/polyval(b,-j*0) );
    case 'bp'  % bandpass
        b = [zeros(1,n) Bw^n zeros(1,n)];
        b = real( b*polyval(den,-j*Wn)/polyval(b,-j*Wn) );
    case 'hp'  % highpass
        b = [1 zeros(1,n)];
        b = real( b*den(1)/b(1) );
    case 'bs'  % bandstop
        r = j*Wn*((-1).^(0:2*n-1)');
        b = poly(r);
        b = real( b*polyval(den,-j*0)/polyval(b,-j*0) );
  end

else

  Wn = 2*atan2(Wn,4);

  switch btype
    case 'lp'  % lowpass
        r = -ones(n,1);
        w = 0;
    case 'bp'  % bandpass
        r = [ones(n,1); -ones(n,1)];
        w = Wn;
    case 'hp'  % highpass
        r = ones(n,1);
        w = pi;
    case 'bs'  % bandstop
        r = exp(j*Wn*( (-1).^(0:2*n-1)' ));
        w = 0;
  end

    b = poly(r);
    % now normalize so |H(w)| == 1:
    kern = exp(-j*w*(0:length(b)-1));
    b = real(b*(kern*den(:))/(kern*b(:)));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = buttzeros(btype,n,Wn,analog)

% This internal function returns more exact zeros.
% Wn input is two element band edge vector

if analog

  % for lowpass and bandpass, don't include zeros at +Inf or -Inf
  switch btype
    case 'lp', z = zeros(0,1);
    case 'bp', z = zeros(n,1);
    case 'hp', z = zeros(n,1);
    case 'bs', z = j*Wn*((-1).^(0:2*n-1)');
  end

else

    Wn = 2*atan2(Wn,4);

  switch btype
    case 'lp', z = -ones(n,1);
    case 'bp', z = [ones(n,1); -ones(n,1)];
    case 'hp', z = ones(n,1);
    case 'bs', z = exp(j*Wn*( (-1).^(0:2*n-1)' ));
  end

end

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [btype,analog,errStr] = iirchk(Wn,varargin)

%IIRCHK  Parameter checking for BUTTER, CHEBY1, CHEBY2, and ELLIP.
%   [btype,analog,errStr] = iirchk(Wn,varargin) returns the 
%   filter type btype (1=lowpass, 2=bandpss, 3=highpass, 4=bandstop)
%   and analog flag analog (0=digital, 1=analog) given the edge
%   frequency Wn (either a one or two element vector) and the
%   optional arguments in varargin.  The variable arguments are 
%   either empty, a one element cell, or a two element cell.
%
%   errStr is empty if no errors are detected; otherwise it contains
%   the error message.  If errStr is not empty, btype and analog
%   are invalid.


Nin = nargin;

errStr = '';

% Define defaults:
analog = 0;    % 0=digital, 1=analog

lw = prod(size(Wn));

if ~( isnumeric(Wn) & any( lw == [ 1  2 ] ) )
    errStr = 'Wn must be a one or two element numeric.';
    return
end

bt = { 'lp' 'bp' };

btype = bt{lw};


if Nin > 3
   errStr = 'Too many input arguments.';
   return
end

% Interpret and strip off trailing 's' or 'z' argument:
if Nin > 1 

   if ~chkcstr(varargin,1)
       errStr = 'FilterType and Option must be Strings.';
       return
   end

  switch lower(varargin{Nin-1})
    case 's'
        analog = 1;
        varargin(Nin-1) = [];
    case 'z'
        analog = 0;
        varargin(Nin-1) = [];
    otherwise
        if Nin > 2
            errStr = 'Analog flag must be either ''z'' or ''s''.';
            return
        end
  end

end

% Check for correct Wn limits
if any( Wn <= 0 )
   errStr = 'The cutoff frequencies must be greater than zero.';
   return
elseif ~analog & any( Wn >= 1 )
   errStr = 'The cutoff frequencies must be within the interval of (0,1).';
   return
end

% At this point, varargin will either be empty, or contain a single
% band type flag.

if length(varargin)==1   % Interpret filter type argument:

  switch lower(varargin{1})
    case { 'lp' 'low'  'lowpass' }
        btype = 'lp';
    case { 'bp' 'pass' 'bandpass' 'band' }
        btype = 'bp';
    case { 'hp' 'high' 'highpass' }
        btype = 'hp';
    case { 'bs' 'stop' 'bandstop' }
        btype = 'bs';
    otherwise
        if Nin == 2
            errStr = 'Option string must be one of "lp" | "bp" | "hp" | "bs" | "z" | "s".';
        else  % nargin == 3
            errStr = 'Filter string must be one of "lp" | "bp" | "hp" | "bs".';
        end
        return
  end

  nr = 1 + any(strcmp(btype,{'bp' 'bs'}));

  if ~( lw == nr )
       errStr = sprintf('For the "%s" filter option, Wn must have %.0f element.',btype,nr);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ok,str] = chkcstr(str,opt)

% CHKCSTR  Checks Input for CellString, contains Strings !
%
%  [ok,str] = chkcstr(str,Option)
%
%  Option ~= 0 ==> CharacterArrays not allowed,
%
%   default: Option == 0   ==>  CharacterArrays --> CellString
%
 
if nargin < 2
   opt = 0;
end

if strcmp(class(str),'char') & isequal(opt,0)
   n = size(str,1);
   if n == 1
      str = strrep(str,char(32),char(1));
   end
   str = cellstr(str);
   if n == 1
      str = strrep(str,char(1),char(32));
   end
end

ok = iscellstr(str);
if ~ok
   return
end

try
  s = cat(2,str{:});
catch
  ok = 0;
  return
end
 
ok = ( strcmp(class(s),'char')  &  ( prod(size(s)) == size(s,2) ) );

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [z,p,k] = buttap(n)

%BUTTAP Butterworth analog lowpass filter prototype.
%   [Z,P,K] = BUTTAP(N) returns the zeros, poles, and gain
%   for an N-th order normalized prototype Butterworth analog
%   lowpass filter.  The resulting filter has N poles around
%   the unit circle in the left half plane, and no zeros.
%
%   See also BUTTER, CHEB1AP, CHEB2AP, ELLIPAP.

%   Author(s): J.N. Little and J.O. Smith, 1-14-87
%   	   L. Shure, 1-13-88, revised
%   (C) 1988-2002 The MathWorks, Inc.
%   $Revision: 1.6 $  $Date: 2002/04/15 01:09:52 $


% Poles are on the unit circle in the left-half plane.
z = [];
p = exp(i*(pi*(1:2:n-1)/(2*n) + pi/2));
p = [p; conj(p)];
p = p(:);
if rem(n,2)==1   % n is odd
    p = [p; -1];
end
k = real(prod(-p));

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [zd, pd, kd, dd] = bilinear(z, p, k, fs, fp, fp1)

%BILINEAR Bilinear transformation with optional frequency prewarping.
%   [Zd,Pd,Kd] = BILINEAR(Z,P,K,Fs) converts the s-domain transfer
%   function specified by Z, P, and K to a z-transform discrete
%   equivalent obtained from the bilinear transformation:
%
%      H(z) = H(s) |
%                  | s = 2*Fs*(z-1)/(z+1)
%
%   where column vectors Z and P specify the zeros and poles, scalar
%   K specifies the gain, and Fs is the sample frequency in Hz.
%
%   [NUMd,DENd] = BILINEAR(NUM,DEN,Fs), where NUM and DEN are 
%   row vectors containing numerator and denominator transfer
%   function coefficients, NUM(s)/DEN(s), in descending powers of
%   s, transforms to z-transform coefficients NUMd(z)/DENd(z).
%
%   [Ad,Bd,Cd,Dd] = BILINEAR(A,B,C,D,Fs) is a state-space version.
%
%   Each of the above three forms of BILINEAR accepts an optional
%   additional input argument that specifies prewarping. 
%
%   For example, [Zd,Pd,Kd] = BILINEAR(Z,P,K,Fs,Fp) applies prewarping 
%   before the bilinear transformation so that the frequency responses
%   before and after mapping match exactly at frequency point Fp
%   (match point Fp is specified in Hz).
%
%   See also IMPINVAR.

%   Author(s): J.N. Little, 4-28-87 
%   	   J.N. Little, 5-5-87, revised
%   (C) 1988-2004 The MathWorks, Inc.
%   $Revision: 1.7.4.2 $  $Date: 2004/12/26 22:15:23 $

%   Gene Franklin, Stanford Univ., motivated the state-space
%   approach to the bilinear transformation.


[mn,nn] = size(z);
[md,nd] = size(p);

if (nd == 1 && nn < 2) && nargout ~= 4	% In zero-pole-gain form
	if mn > md
		error('Numerator cannot be higher order than denominator.')
	end
	if nargin == 5		% Prewarp
		fp = 2*pi*fp;
		fs = fp/tan(fp/fs/2);
	else
		fs = 2*fs;
	end
	z = z(finite(z));	 % Strip infinities from zeros
	pd = (1+p/fs)./(1-p/fs); % Do bilinear transformation
	zd = (1+z/fs)./(1-z/fs);
% real(kd) or just kd?
	kd = (k*prod(fs-z)./prod(fs-p));
	zd = [zd;-ones(length(pd)-length(zd),1)];  % Add extra zeros at -1

elseif (md == 1 && mn == 1) || nargout == 4 %
	if nargout == 4		% State-space case
		a = z; b = p; c = k; d = fs; fs = fp;
		error(abcdchk(a,b,c,d));
		if nargin == 6			% Prewarp
			fp = fp1;		% Decode arguments
			fp = 2*pi*fp;
			fs = fp/tan(fp/fs/2)/2;
		end
	else			% Transfer function case
		if nn > nd
			error('Numerator cannot be higher order than denominator.')
		end
		num = z; den = p;		% Decode arguments
		if nargin == 4			% Prewarp
			fp = fs; fs = k;	% Decode arguments
			fp = 2*pi*fp;
			fs = fp/tan(fp/fs/2)/2;
		else
			fs = k;			% Decode arguments
		end
		% Put num(s)/den(s) in state-space canonical form.  
		[a,b,c,d] = tf2ss(num,den);
	end
	% Now do state-space version of bilinear transformation:
	t = 1/fs;
	r = sqrt(t);
	t1 = eye(size(a)) + a*t/2;
	t2 = eye(size(a)) - a*t/2;
	ad = t2\t1;
	bd = t/r*(t2\b);
	cd = r*c/t2;
	dd = c/t2*b*t/2 + d;
	if nargout == 4
		zd = ad; pd = bd; kd = cd;
	else
		% Convert back to transfer function form:
		p = poly(ad);
		zd = poly(ad-bd*cd)+(dd-1)*p;
		pd = p;
	end
else
	error('First two arguments must have the same orientation.')
end

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [at,bt,ct,dt] = lp_transf(mode,a,b,c,d,wo,bw)

% LP_TRANS Lowpass to low/band/high/stop - pass Transformation
% LP_TRANS(Mode,a,b,c,d,Wo,[Bw])
% LP_TRANS(Mode,NUM,DEN,Wo,[Bw])
%
%LP2LP Lowpass to lowpass analog filter transformation.
%LP2BP Lowpass to bandpass analog filter transformation.
%LP2HP Lowpass to highpass analog filter transformation.
%LP2BS Lowpass to bandstop analog filter transformation.

mode = lower(mode);

chk = ( mode(1) == 'b' );  % BandPass or BandStop

tfc = ( nargin == 3+chk ); % Transfer function case

if tfc

   % convert column vector to rows
   if ( size(a,2) == 1 ), a = a(:).'; end
   if ( size(b,2) == 1 ), b = b(:).'; end

   % Transform to state-space
   wo = c; if chk, bw = d; end 

   [a,b,c,d] = tf2ss(a,b);

end

error(abcdchk(a,b,c,d));

nb = size(b,2);  [mc,ma] = size(c);

switch mode

  case 'lp'

    % Transform lowpass to lowpass
    at = wo*a;
    bt = wo*b;
    ct = c;
    dt = d;

  case 'hp'

    % Transform lowpass to highpass
    at =  wo*inv(a);
    bt = -wo*(a\b);
    ct = c/a;
    dt = d - c/a*b;

  case 'bp'

    % Transform lowpass to bandpass
    q = wo/bw;
    at = wo*[a/q eye(ma); -eye(ma) zeros(ma)];
    bt = wo*[b/q; zeros(ma,nb)];
    ct = [c zeros(mc,ma)];
    dt = d;

  case 'bs'

    % Transform lowpass to bandstop
    q = wo/bw;
    at =  [wo/q*inv(a) wo*eye(ma); -wo*eye(ma) zeros(ma)];
    bt = -[wo/q*(a\b); zeros(ma,nb)];
    ct = [c/a zeros(mc,ma)];
    dt = d - c/a*b;

end

if tfc % Transform back to transfer function
   [z,k] = tzero(at,bt,ct,dt);
   num = k * poly(z);
   den = poly(at);
   at = num;
   bt = den;
end

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,A,B,C,D] = abcdchk(A,B,C,D)

%ABCDCHK Checks dimensional consistency of A,B,C,D matrices.
%   ERROR(ABCDCHK(A,B,C,D)) checks that the dimensions of A,B,C,D
%   are consistent for a linear, time-invariant system model.
%   An error occurs if the nonzero dimensions are not consistent.
%
%   [MSG,A,B,C,D] = ABCDCHK(A,B,C,D) also alters the dimensions
%   any 0-by-0 empty matrices to make them consistent with the others.

%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.22.4.2 $  $Date: 2004/01/24 09:22:39 $

if nargin < 4, D = []; end
if nargin < 3, C = []; end
if nargin < 2, B = []; end
if nargin < 1, A = []; end

[ma,na] = size(A);
[mb,nb] = size(B);
[mc,nc] = size(C);
[md,nd] = size(D);

if mc==0 && nc==0 && (md==0 || na==0)
   mc = md; nc = na; C = zeros(mc,nc);
end
if mb==0 && nb==0 && (ma==0 || nd==0)
   mb = ma; nb = nd; B = zeros(mb,nb);
end
if md==0 && nd==0 && (mc==0 || nb==0)
   md = mc; nd = nb; D = zeros(md,nd);
end
if ma==0 && na==0 && (mb==0 || nc==0)
   ma = mb; na = nc; A = zeros(ma,na);
end

if ma~=na && nargin>=1
    msg = makeMsg('The A matrix must be square', ...
                  'AMustBeSquare');
elseif ma~=mb && nargin>=2
    msg = makeMsg('The A and B matrices must have the same number of rows.', ...
                  'AAndBNumRowsMismatch');
elseif na~=nc && nargin>=3
    msg = makeMsg('The A and C matrices must have the same number of columns.', ...
                  'AAndCNumColumnsMismatch');
elseif md~=mc && nargin>=4
    msg = makeMsg('The C and D matrices must have the same number of rows.', ...
                  'CAndDNumRowsMismatch');
elseif nd~=nb && nargin>=4
    msg = makeMsg('The B and D matrices must have the same number of columns.', ...
                  'BAndDNumColumnsMismatch');
else
    msg.message = '';
    msg.identifier = '';
    msg = msg(zeros(0,1));
end

function msg = makeMsg(message, identifier)
msg.message = message;
msg.identifier = ['MATLAB:abcdchk:' identifier];

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function c = poly(x)

%POLY Convert roots to polynomial.
%   POLY(A), when A is an N by N matrix, is a row vector with
%   N+1 elements which are the coefficients of the
%   characteristic polynomial, DET(lambda*EYE(SIZE(A)) - A) .
%
%   POLY(V), when V is a vector, is a vector whose elements are
%   the coefficients of the polynomial whose roots are the
%   elements of V . For vectors, ROOTS and POLY are inverse
%   functions of each other, up to ordering, scaling, and
%   roundoff error.
%
%   ROOTS(POLY(1:20)) generates Wilkinson's famous example.
%
%   Class support for inputs A,V:
%      float: double, single
%
%   See also ROOTS, CONV, RESIDUE, POLYVAL.

%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 5.14.4.2 $  $Date: 2004/03/02 21:47:54 $

[m,n] = size(x);
if m == n
   % Characteristic polynomial (square x)
   e = eig(x);
elseif (m==1) || (n==1)
   e = x;
else
   error('MATLAB:poly:InputSize','Argument must be a vector or a square matrix.')
end

% Strip out infinities
e = e( isfinite(e) );

% Expand recursion formula
n = length(e);
c = [1 zeros(1,n,class(x))];
for j=1:n
    c(2:(j+1)) = c(2:(j+1)) - e(j).*c(1:j);
end

% The result should be real if the roots are complex conjugates.
if isequal(sort(e(imag(e)>0)),sort(conj(e(imag(e)<0))))
    c = real(c);
end

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [z,gain] = tzero(a,b,c,d)
%TZERO  Transmission zeros of LTI systems.
% 
%    Z = TZERO(SYS) returns the transmission zeros of the LTI 
%    system SYS.
%
%    [Z,GAIN] = TZERO(SYS) also returns the transfer function 
%    gain if the system is SISO.
%   
%    Z = TZERO(A,B,C,D) works directly on the state space matrices
%    and returns the transmission zeros of the state-space system:   
%             .
%             x = Ax + Bu     or   x[n+1] = Ax[n] + Bu[n]
%             y = Cx + Du           y[n]  = Cx[n] + Du[n]
%
%    See also PZMAP, POLE, EIG.

%   Clay M. Thompson  7-23-90
%       Revised: A.Potvin 6-1-94, P.Gahinet 5-15-96
%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 1.35.4.2 $  $Date: 2004/01/24 09:22:49 $

% Extracts from the system matrix of a state-space system [ A B C D ] a regular
% pencil [ Lambda*Bf - Af ] which has the NU Invariant Zeros of the system as
% Generalized Eigenvalues.
%  
%  Reference: Adapted from "Computation of Zeros of Linear Multivariable
%             Systems", A. Emami-Naeini, and P. Van Dooren; Automatica 
%             Vol. 18, No. 4, pp. 415-430, 1982.

% The transmission zero calculation can be tested by checking that the
% matrix: [A B;C D]-lambda*[I 0; 0 0] is rank deficient, where lambda
% is an element of Z.

ni = nargin;
error(nargchk(4,4,ni));
[msg,a,b,c,d]=abcdchk(a,b,c,d); error(msg);
gain = [];

if isempty(b) && isempty(c),
   z = zeros(0,1);
   gain = d; 
   return
end

% Epsilon used for transmission zero calculation.
Zeps = 10*eps*norm(a,'fro');

nn = size(a,1);
pp = size(c,1);
mm = size(b,2);
issiso = (mm==1 & pp==1);

% Construct the Compound Matrix [ B A ] of Dimension (N+P)x(M+N) 
%                               [ D C ]
bf = [b a; d c];

% Reduce this system to one with the same invariant zeros and with
% D(*) full rank MU (The Normal Rank of the original system)
[bf,mu,nu] = tzreduce(bf,mm,nn,pp,Zeps,pp,0);
Rank = mu;
if nu==0, 
   z = zeros(0,1);
else
   % Pretranspose the system
   mnu = mm+nu;
   numu = nu+mu;
   af = zeros(mnu,numu);
   af(mnu:-1:1,numu:-1:1) = bf(1:numu,1:mnu).';

   if mu~=mm,
      pp = mm;
      nn = nu;
      mm = mu;

      % Reduce the system to one with the same invariant zeros and with
      % D(*) square invertible

      [af,mu,nu] = tzreduce(af,mm,nn,pp,Zeps,pp-mm,mm);
      mnu = mm+nu;
   end

   if nu==0,
      z = zeros(0,1);
   else
      % Perform a unitary transformation on the columns of [ sI-A B ]
      %                          [ sBf-Af X ]              [   -C D ]
      % in order to reduce it to [   0    Y ] with Y & Bf square invertible
      bf(1:nu,1:mnu) = [zeros(nu,mm),eye(nu)];
      if Rank~=0,
         nu1 = nu+1;
         i1 = nu+mu;
         i0 = mm;
         for i=1:mm
            i0  = i0-1;
            cols = i0 + (1:nu1);
            [dummy,s,zero] = housh(af(i1,cols)',nu1,Zeps);
%REVISIT: temp. fix
%            af(1:i1,cols) = af(1:i1,cols)*(eye(nu1)-s*dummy*dummy');
%            bf(1:nu,cols) = bf(1:nu,cols)*(eye(nu1)-s*dummy*dummy');
            af(1:i1,cols) = af(1:i1,cols)-s*(af(1:i1,cols)*dummy)*dummy';
            bf(1:nu,cols) = bf(1:nu,cols)-s*(bf(1:nu,cols)*dummy)*dummy';
            i1 = i1-1;
         end % for
      end % if Rank~=0

      % Solve Generalized zeros of sBF - AF
      z = eig(af(1:nu,1:nu),bf(1:nu,1:nu));
   end % if nu==0
end

% Compute transfer function gain if necessary
if (nargout==2) && issiso,  % (mm*pp==1),
   if nu==nn,
      gain=bf(nu+1,1);
   else
      gain=bf(nu+1,1)*prod(diag(bf(nu+2:nn+1,nu+2:nn+1)));
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [abcd,mu,nu] = tzreduce(abcd,m,n,p,Zeps,ro,sigma)
%TZREDUCE Utility function for TZERO.
%   [ABCD.MU,NU] = TZREDUCE(ABCD,M,N,P,Zeps,RO,SIGMA) Extracts
%       the reduced system from the input system matrix 
%   such that the new system D matrix has full row rank.

%   Clay M. Thompson  7-23-90

%  Extracts from the (N+P)x(M+N) system [ B A ],  a (NU+MU)x(M+NU) 'reduced' 
%         [ B' A' ]                     [ D C ]
%  system [ D' C' ] having the same transmission zeros but with D' Full Row 
%  Rank.  The system [ A' B' C' D' ] overwrites the old system.
%
% Reference: Adapted from "Computation of Zeros of Linear Multivariable
%            Systems", A. Emami-Naeini, and P. Van Dooren; Automatica
%            Vol. 18, No. 4, pp. 415-430, 1982.

% Initialize
Sum2 = zeros(1,max(p,m));
mu = p;
nu = n;
while mu~=0,
   ro1 = ro;
   mnu = m+nu;
   numu = nu+mu;

   if m~=0,
      ro1 = ro1+1;
      irow = nu;

      % Compress rows of D(*).  First exploit triangular shape
      for icol=1:sigma-1
         rows = irow + (1:ro1);
         [dummy,s,zero] = housh(abcd(rows,icol),1,Zeps);
% REVISIT: temp. fix
%         abcd(rows,icol:mnu) = (eye(ro1)-s*dummy*dummy')*abcd(rows,icol:mnu);
         abcd(rows,icol:mnu) = abcd(rows,icol:mnu)-s*dummy*(dummy'*abcd(rows,icol:mnu));
         irow = irow+1;
      end

      % Continue householder with pivoting
      if sigma==0,
         sigma = 1;
         ro1 = ro1-1;
      end

      if sigma~=m,
         Sum2(sigma:m) = sum(abcd(irow+1:irow+ro1,sigma:m).*conj(abcd(irow+1:irow+ro1,sigma:m)));
      end

      for icol=sigma:m;
         % Pivot if necessary
         if icol~=m,
            Rows = 1:numu;
            [dum,ibar] = max(Sum2(icol:m));
            ibar = ibar+icol-1;
            if ibar~=icol,
               Sum2(ibar) = Sum2(icol); 
               Sum2(icol) = dum;
               dum = abcd(Rows,icol);
               abcd(Rows,icol)=abcd(Rows,ibar);
               abcd(Rows,ibar)=dum;
            end
         end

         % Perform Householder transformation
         [dummy,s,zero] = housh(abcd(irow+1:irow+ro1,icol),1,Zeps);
         if zero,
            break
         end
         if ro1==1,
            return
         end
% REVISIT: temp. fix
         abcd(irow+1:irow+ro1,icol:mnu) = abcd(irow+1:irow+ro1,icol:mnu) - ...
                             s*dummy*(dummy'*abcd(irow+1:irow+ro1,icol:mnu));
         irow = irow+1;
         ro1 = ro1-1;
         Sum2(icol:m) = Sum2(icol:m) - abcd(irow,icol:m) .* conj(abcd(irow,icol:m));
      end % for
   end % if
   tau = ro1;
   sigma = mu-tau;

   % Compress the columns of C(*)
   if (nu<=0),
      mu = sigma; 
      nu = 0;
      return
   end

   i1 = nu+sigma;
   mm1 = m+1;
   n1 = nu;
   if tau~=1,
      Sum2(1:tau) = sum((abcd(i1+1:i1+tau,mm1:mnu).*conj(abcd(i1+1:i1+tau,mm1:mnu)))');
   end

   for ro1=1:tau;
      ro = ro1-1;
      i = tau-ro;
      i2 = i+i1;

      % Pivot if necessary
      if i~=1,
         [dum,ibar] = max(Sum2(1:i));
         if ibar~=i,
            Sum2(ibar) = Sum2(i); Sum2(i) = dum;
            dum = abcd(i2,mm1:mnu);
            abcd(i2,mm1:mnu) = abcd(ibar+i1,mm1:mnu);
            abcd(ibar+i1,mm1:mnu) = dum;
         end
      end

      % Perform Householder Transformation 
      cols = m + (1:n1);
      [dummy,s,zero] = housh(abcd(i2,cols)',n1,Zeps);
      if zero,
         break
      end
      if n1~=1
% REVISIT: temp. fix
%         abcd(1:i2,cols) = abcd(1:i2,cols)*(eye(n1)-s*dummy*dummy');
         abcd(1:i2,cols) = abcd(1:i2,cols)-s*(abcd(1:i2,cols)*dummy)*dummy';
         mn1 = m+n1;
%         abcd(1:n1,1:mn1) = (eye(n1)-s*dummy*dummy')*abcd(1:n1,1:mn1);
         abcd(1:n1,1:mn1) = abcd(1:n1,1:mn1)-s*dummy*(dummy'*abcd(1:n1,1:mn1));
         Sum2(1:i) = Sum2(1:i)-(abcd(i1+1:i1+i,mn1) .* conj(abcd(i1+1:i1+i,mn1)))';
         mnu = mnu-1;
      end
      n1 = n1-1;
   end % for

   if ~zero,
      ro = tau;
   end
   nu = nu-ro;
   mu = sigma+ro;
  
   if ro==0,
      return
   end
end % while

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u,s,zero] = housh(u,j,heps)
%HOUSH  Construct a householder transformation H=I-s*UU'.  Used in TZERO.
%   [U,S,ZERO] = HOUSH(U,J,Heps)

%   Clay M. Thompson  7-23-90

%  Constructs a Householder transformation H=I-s*UU' that 'mirrors' a 
%  vector u to the Jth unit vector.  If NORM(U)<Eps then Zero=1 [True]
%
% Reference: Adapted from "Computation of Zeros of Linear Multivariable
%            Systems", A. Emami-Naeini, and P. Van Dooren; Automatica
%            Vol. 18, No. 4, pp. 415-430, 1982.

s = sum(u.*conj(u));
alfa = sqrt(s);
if (alfa<=heps), 
   zero=1; 
   return
end
zero=0;

% Transform is I-2vv'/(v'v) where v = u+beta*ej and beta = (uj/|uj|)*norm(u)
beta = (sign(u(j))+(u(j)==0)) * alfa;
s = 1./(s+real(conj(beta)*u(j)));   % u'*u+conj(beta)*uj
u(j) = u(j)+beta;

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [num, den] = ss2tf(a,b,c,d,iu)
%SS2TF  State-space to transfer function conversion.
%   [NUM,DEN] = SS2TF(A,B,C,D,iu)  calculates the transfer function:
%
%               NUM(s)          -1
%       H(s) = -------- = C(sI-A) B + D
%               DEN(s)
%   of the system:
%       .
%       x = Ax + Bu
%       y = Cx + Du
%
%   from the iu'th input.  Vector DEN contains the coefficients of the
%   denominator in descending powers of s.  The numerator coefficients
%   are returned in matrix NUM with as many rows as there are 
%   outputs y.
%
%   See also TF2SS, ZP2TF, ZP2SS.

%   J.N. Little 4-21-85
%   Revised 7-25-90 Clay M. Thompson, 10-11-90 A.Grace
%   Copyright 1984-2003 The MathWorks, Inc. 
%   $Revision: 1.24.4.1 $  $Date: 2003/05/01 20:42:46 $

error(nargchk(4,5,nargin));
[msg,a,b,c,d]=abcdchk(a,b,c,d); error(msg);

[mc,nu] = size(d);
if nargin==4,
  if (nu<=1)
    iu = 1;
  else
    error('MATLAB:ss2tf:NeedIU',...
          'IU must be specified for systems with more than one input.');
  end
end

den = poly(a);
if ~isempty(b), b = b(:,iu); end
if ~isempty(d), d = d(:,iu); end

% System is just a gain or it has only a denominator:
if isempty(b) && isempty(c)
    num = d;
    if isempty(d) && isempty(a)
        den = [];
    end
    return;
end

nc = length(a);
num = ones(mc, nc+1);
for i=1:mc
    num(i,:) = poly(a-b*c(i,:)) + (d(i) - 1) * den;
end

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [z,p,k] = ss2zp(a,b,c,d,iu)
%SS2ZP  State-space to zero-pole conversion.
%   [Z,P,K] = SS2ZP(A,B,C,D,IU)  calculates the transfer function in
%   factored form:
%
%                -1           (s-z1)(s-z2)...(s-zn)
%       H(s) = C(sI-A) B + D =  k ---------------------
%                             (s-p1)(s-p2)...(s-pn)
%   of the system:
%       .
%       x = Ax + Bu
%       y = Cx + Du
%
%   from the single input IU.  The vector P contains the pole 
%   locations of the denominator of the transfer function.  The 
%   numerator zeros are returned in the columns of matrix Z with as 
%   many columns as there are outputs y.  The gains for each numerator
%   transfer function are returned in column vector K.
%
%   See also ZP2SS,PZMAP,TZERO, EIG.

%   J.N. Little 7-17-85
%   Revised 3-12-87 JNL, 8-10-90 CLT, 1-18-91 ACWG, 2-22-94 AFP.
%   Copyright 1984-2003 The MathWorks, Inc.
%   $Revision: 1.33.4.1 $  $Date: 2003/05/01 20:42:47 $

error(nargchk(4,5,nargin));
[msg,a,b,c,d]=abcdchk(a,b,c,d); error(msg);

[nx,ns] = size(a);

if nargin==4,
    if nx>0,
        [nb,nu] = size(b);
    else
        [ny,nu] = size(d);
    end
    if (nu<=1), 
        iu = 1;
    else
        error('MATLAB:ss2zp:NeedIU',...
              'IU must be specified for systems with more than one input.');
    end
end

% Remove relevant input:
if ~isempty(b), b = b(:,iu); end
if ~isempty(d), d = d(:,iu); end

% Trap gain-only models
if nx==0 && ~isempty(d), z = []; p = []; k = d; return, end

% Do poles first, they're easy:
p = eig(a);

% Compute zeros and gains using transmission zero calculation
% Took out check for tzreduce since that now ships with SP
[ny,nu] = size(d);
z = [];
k = zeros(ny,1);
for i=1:ny
   [zi,gi] = tzero(a,b,c(i,:),d(i,:));
   [mz,nz] = size(z);
   nzi = length(zi);
   if i==1,
      z = zi;
   else
      linf = inf;
      z = [[z; linf(ones(max(0,nzi-mz),1),ones(max(nz,1),1))], ...
          [zi;linf(ones(max(0,mz-nzi),1),1)]];
  end
  k(i) = gi;
end

% end ss2zp.m

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [a,b,c,d] = tf2ss(num, den)
%TF2SS  Transfer function to state-space conversion.
%   [A,B,C,D] = TF2SS(NUM,DEN)  calculates the state-space 
%   representation:
%       .
%       x = Ax + Bu
%       y = Cx + Du
%
%   of the system:
%               NUM(s) 
%       H(s) = --------
%               DEN(s)
%
%   from a single input.  Vector DEN must contain the coefficients of
%   the denominator in descending powers of s.  Matrix NUM must 
%   contain the numerator coefficients with as many rows as there are
%   outputs y.  The A,B,C,D matrices are returned in controller 
%   canonical form.  This calculation also works for discrete systems.
%
%   For discrete-time transfer functions, it is highly recommended to
%   make the length of the numerator and denominator equal to ensure 
%   correct results.  You can do this using the function EQTFLENGTH in
%   the Signal Processing Toolbox.  However, this function only handles
%   single-input single-output systems.
%
%   See also TF2ZP, SS2TF, ZP2SS, ZP2TF.

%   J.N. Little 3-24-85
%   Copyright 1984-2003 The MathWorks, Inc.
%   $Revision: 1.22.4.1 $  $Date: 2003/05/01 20:42:48 $
%   Latest revision 4-29-89 JNL, 7-29-96 PG

[mnum,nnum] = size(num);
[mden,n] = size(den);
% Check for null systems
if  (n == 0 && nnum == 0), a=[]; b=[]; c=[]; d=[]; return, end

if min(mden,n)>1,
   % Error out if DEN is an array
   error('MATLAB:tf2ss:NeedRowDenom', 'Denominator must be a row vector.');
elseif mden>1,
   % Transpose DEN when a column vector
   den = den.';
end

% Strip leading zeros from denominator
inz = find(den ~= 0);
den = den(inz(1):end);
[mden,n] = size(den);

% Check for proper numerator
if nnum > n
    % Try to strip leading zeros to make proper
    if (all(all(num(:,1:(nnum-n)) == 0)))
        num = num(:,(nnum-n+1):nnum);
        [mnum,nnum] = size(num);
    else
        error('MATLAB:tf2ss:DenomInvalidOrder',...
              ['Order of denominator must be greater than or equal to ',...
	      'order of numerator.']);
    end
end

% Pad numerator with leading zeros, to make it have the same number of
% columns as the denominator, and normalize it to den(1)
num = [zeros(mnum,n-nnum) num]./den(1);

% Do the D-matrix first
if length(num)
    d = num(:,1);
else
    d = [];
end

% Handle special constant case:
if n == 1
    a = [];
    b = [];
    c = [];
    return
end

% Now do the rest, starting by normalizing den to den(1),
den = den(2:n) ./ den(1);
a = [-den; eye(n-2,n-1)];
b = eye(n-1,1);
if mnum > 0
    c = num(:,2:n) - num(:,1) * den;
else
    c = [];
end

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [a,b,c,d] = zp2ss(z,p,k)
%ZP2SS  Zero-pole to state-space conversion.
%   [A,B,C,D] = ZP2SS(Z,P,K)  calculates a state-space representation:
%       .
%       x = Ax + Bu
%       y = Cx + Du
%
%   for a system given a set of pole locations in column vector P,
%   a matrix Z with the zero locations in as many columns as there are
%   outputs, and the gains for each numerator transfer function in
%   vector K.  The A,B,C,D matrices are returned in block diagonal
%   form.  
%
%   The poles and zeros must correspond to a proper system. If the poles
%   or zeros are complex, they must appear in complex conjugate pairs,
%   i.e., the corresponding transfer function must be real.
%     
%   See also SS2ZP, ZP2TF, TF2ZP, TF2SS, SS2TF.

%   J.N. Little & G.F. Franklin  8-4-87
%   Revised 12-27-88 JNL, 12-8-89, 11-12-90, 3-22-91, A.Grace, 7-29-96 P. Gahinet
%   Copyright 1984-2003 The MathWorks, Inc.
%   $Revision: 1.27.4.1 $  $Date: 2003/05/01 20:42:54 $

[z,p,k,isSIMO,msg] = parse_input(z,p,k);
error(msg);
     
if isSIMO
    % If it's multi-output, we can't use the nice algorithm
    % that follows, so use the numerically unreliable method
    % of going through polynomial form, and then return.
    [num,den] = zp2tf(z,p,k); % Suppress compile-time diagnostics
    [a,b,c,d] = tf2ss(num,den);
    return
end

% Strip infinities and throw away.
p = p(isfinite(p));
z = z(isfinite(z));

% Group into complex pairs
np = length(p);
nz = length(z);
z = cplxpair(z,1e6*nz*norm(z)*eps + eps);
p = cplxpair(p,1e6*np*norm(p)*eps + eps);

% Initialize state-space matrices for running series
a=[]; b=zeros(0,1); c=ones(1,0); d=1;

% If odd number of poles AND zeros, convert the pole and zero
% at the end into state-space.
%   H(s) = (s-z1)/(s-p1) = (s + num(2)) / (s + den(2))
if rem(np,2) & rem(nz,2)
    a = p(np);
    b = 1;
    c = p(np) - z(nz);
    d = 1;
    np = np - 1;
    nz = nz - 1;
end

% If odd number of poles only, convert the pole at the
% end into state-space.
%  H(s) = 1/(s-p1) = 1/(s + den(2)) 
if rem(np,2)
    a = p(np);
    b = 1;
    c = 1;
    d = 0;
    np = np - 1;
end 

% If odd number of zeros only, convert the zero at the
% end, along with a pole-pair into state-space.
%   H(s) = (s+num(2))/(s^2+den(2)s+den(3)) 
if rem(nz,2)
    num = real(poly(z(nz)));
    den = real(poly(p(np-1:np)));
    wn = sqrt(prod(abs(p(np-1:np))));
    if wn == 0, wn = 1; end
    t = diag([1 1/wn]); % Balancing transformation
    a = t\[-den(2) -den(3); 1 0]*t;
    b = t\[1; 0];
    c = [1 num(2)]*t;
    d = 0;
    nz = nz - 1;
    np = np - 2;
end

% Now we have an even number of poles and zeros, although not 
% necessarily the same number - there may be more poles.
%   H(s) = (s^2+num(2)s+num(3))/(s^2+den(2)s+den(3))
% Loop thru rest of pairs, connecting in series to build the model.
i = 1;
while i < nz
    index = i:i+1;
    num = real(poly(z(index)));
    den = real(poly(p(index)));
    wn = sqrt(prod(abs(p(index))));
    if wn == 0, wn = 1; end
    t = diag([1 1/wn]); % Balancing transformation
    a1 = t\[-den(2) -den(3); 1 0]*t;
    b1 = t\[1; 0];
    c1 = [num(2)-den(2) num(3)-den(3)]*t;
    d1 = 1;
%   [a,b,c,d] = series(a,b,c,d,a1,b1,c1,d1); 
% Next lines perform series connection 
    [ma1,na1] = size(a);
    [ma2,na2] = size(a1);
    a = [a zeros(ma1,na2); b1*c a1];
    b = [b; b1*d];
    c = [d1*c c1];
    d = d1*d;

    i = i + 2;
end

% Take care of any left over unmatched pole pairs.
%   H(s) = 1/(s^2+den(2)s+den(3))
while i < np
    den = real(poly(p(i:i+1)));
    wn = sqrt(prod(abs(p(i:i+1))));
    if wn == 0, wn = 1; end
    t = diag([1 1/wn]); % Balancing transformation
    a1 = t\[-den(2) -den(3); 1 0]*t;
    b1 = t\[1; 0];
    c1 = [0 1]*t;
    d1 = 0;
%   [a,b,c,d] = series(a,b,c,d,a1,b1,c1,d1);
% Next lines perform series connection 
    [ma1,na1] = size(a);
    [ma2,na2] = size(a1);
    a = [a zeros(ma1,na2); b1*c a1];
    b = [b; b1*d];
    c = [d1*c c1];
    d = d1*d;

    i = i + 2;
end

% Apply gain k:
c = c*k;
d = d*k;

%----------------------------------------------------------------------------
function [z,p,k,isSIMO,msg] = parse_input(z,p,k)
%PARSE_INPUT   Make sure input args are valid.

% Initially assume it is a SISO system
isSIMO = 0;
msg = '';

% Check that p is a vector
if ~any(size(p)<2),
    msg = 'You must specify a vector of poles.';
    return
end
% Columnize p
p = p(:);

% Check that k is a vector
if ~any(size(k)<2),
    msg = 'The gain must be a scalar or a vector.';
    return
end
% Columnize k
k = k(:);


% Check size of z
if any(size(z)<2),
    % z is a vector or an empty, columnize it
    z = z(:);
else
    % z is a matrix
    isSIMO = 1;
end

% Check for properness
if size(z,1) > length(p),
    % improper
    msg = 'Must be a proper system.';
    return
end

% Check for the appropriate length of k
if length(k) ~= size(z,2) && (~isempty(z))
    msg = 'The length of K must be equal to the number of columns of Z.';
end

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function z = cplxpair(x,tol,dim)
%CPLXPAIR Sort numbers into complex conjugate pairs.
%   Y = CPLXPAIR(X) takes a vector of complex conjugate pairs and/or
%   real numbers.  CPLXPAIR rearranges the elements of X so that
%   complex numbers are collected into matched pairs of complex
%   conjugates.  The pairs are ordered by increasing real part.
%   Any purely real elements are placed after all the complex pairs.
%   Y = CPLXPAIR(X,TOL) uses a relative tolerance of TOL for
%   comparison purposes.  The default is TOL = 100*EPS.
%
%   For X an N-D array, CPLXPAIR(X) and CPLXPAIR(X,TOL) rearranges
%   the elements along the first non-singleton dimension of X.
%   CPLXPAIR(X,[],DIM) and CPLXPAIR(X,TOL,DIM) sorts X along 
%   dimension DIM.

%   L. Shure 1-27-88
%   Revised 4-30-96 D. Orofino
%   Copyright 1984-2003 The MathWorks, Inc. 
%   $Revision: 5.14.4.1 $  $Date: 2003/05/01 20:41:23 $

if isempty(x),
  z = x; return  % Quick exit if empty input
end
if nargin == 3,
  nshifts = 0;
  perm = [dim:max(ndims(x),dim) 1:dim-1];
  x = permute(x,perm);
else
  [x,nshifts] = shiftdim(x);
  perm = [];
end

% Supply defaults for tolerance:
if nargin<2 || isempty(tol), tol=100*eps; end

% Reshape x to a 2-D matrix:
xsiz   = size(x);         % original shape of input
x      = x(:,:);          % reshape to a 2-D matrix
z      = zeros(size(x));  % preallocate temp storage
errmsg = 'Complex numbers can''t be paired.';

for k = 1:size(x,2),
  % Get next column of x:
  xc = x(:,k);

  % Find purely-real entries:
  idx = find(abs(imag(xc)) <= tol*abs(xc));
  nr = length(idx);     % Number of purely-real entries
  if ~isempty(idx),
    % Store sorted real's at end of column:
    z(end-nr+1:end,k)  = sort(real(xc(idx)));
    xc(idx) = [];  % Remove values from current column
  end

  nc = length(xc); % Number of entries remaining in input column
  if nc>0,
    % Complex values in list:
    if rem(nc,2)==1
      % Odd number of entries remaining
      error('MATLAB:cplxpair:ComplexValuesPaired',errmsg);  
    end

    % Sort complex column-vector xc, based on its real part:
    [xtemp,idx] = sort(real(xc));
    xc = xc(idx);  % Sort complex numbers based on real part

    % Check if real parts occur in pairs:
    %   Compare to xc() so imag part is considered (in case real part is nearly 0).
    %   Arbitrary choice of using abs(xc(1:2:nc)) or abs(xc(2:2:nc)) for tolerance
    if any( abs(xtemp(1:2:nc)-xtemp(2:2:nc)) > tol.*abs(xc(1:2:nc)) ),
      error('MATLAB:cplxpair:ComplexValuesPaired', errmsg);
    end

    % Check real part pairs to see if imag parts are conjugates:
    nxt_row = 1;  % next row in z(:,k) for results
    while ~isempty(xc),
      % Find all real parts "identical" to real(xc(1)):
      idx = find( abs(real(xc) - real(xc(1))) <= tol.*abs(xc) );
      nn = length(idx); % # of values with identical real parts
      if nn<=1
        % Only 1 value found - certainly not a pair!
        error('MATLAB:cplxpair:ComplexValuesPaired', errmsg); 
      end

      % There could be multiple pairs with "identical" real parts.
      % Sort the imaginary parts of those values with identical real
      %   parts - these SHOULD be the next N entries, with N even.
      [xtemp,idx] = sort(imag(xc(idx)));
      xq = xc(idx);  % Get complex-values with identical real parts,
      % which are now sorted by imaginary component.
      % Verify conjugate-pairing of imaginary parts:
      if any( abs(xtemp + xtemp(nn:-1:1)) > tol.*abs(xq) ),
        error('MATLAB:cplxpair:ComplexValuesPaired', errmsg);
      end
      % Keep value with pos imag part, and compute conjugate for pair.
      % List value with most-neg imag first, then its conjugate.
      z(nxt_row : nxt_row+nn-1, k) = reshape([conj(xq(end:-1:nn/2+1)) ...
                                                   xq(end:-1:nn/2+1)].',nn,1);
      nxt_row = nxt_row+nn;  % Bump next-row pointer
      xc(idx) = [];          % Remove entries from xc
    end

  end % of complex-values check
end % of column loop

% Reshape Z to appropriate form
z = reshape(z,xsiz);
if ~isempty(perm),
  z = ipermute(z,perm);
end
if nshifts~=0,
  z = shiftdim(z,-nshifts);
end

% end of cplxpair.m

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [b,nshifts] = shiftdim(x,n)
%SHIFTDIM Shift dimensions.
%   B = SHIFTDIM(X,N) shifts the dimensions of X by N.  When N is
%   positive, SHIFTDIM shifts the dimensions to the left and wraps the
%   N leading dimensions to the end.  When N is negative, SHIFTDIM
%   shifts the dimensions to the right and pads with singletons.
%
%   [B,NSHIFTS] = SHIFTDIM(X) returns the array B with the same
%   number of elements as X but with any leading singleton 
%   dimensions removed. NSHIFTS returns the number of dimensions 
%   that are removed. If X is a scalar, SHIFTDIM has no effect.
%
%   SHIFTDIM is handy for creating functions that, like SUM
%   or DIFF, work along the first non-singleton dimension.
%
%   Examples:
%       a = rand(1,1,3,1,2);
%       [b,n]  = shiftdim(a); % b is 3-by-1-by-2 and n is 2.
%       c = shiftdim(b,-n);   % c == a.
%       d = shiftdim(a,3);    % d is 1-by-2-by-1-by-1-by-3.
%
%   See also CIRCSHIFT, RESHAPE, SQUEEZE.

%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 1.15.4.2 $  $Date: 2004/03/02 21:47:08 $


siz = size(x);
if nargin==1
  n = find(siz~=1,1,'first')-1; % Find leading singleton dimensions
end

if n > 0  % Wrapped shift to the left
  n = rem(n,ndims(x));
end 

if  isempty(n) || isequal(n,0)
  b = x;    % Quick exit if no shift required
  nshifts = 0;
elseif (n > 0)  
  if isequal(siz(1:n),ones(1,n))
    s = [siz ones(1,n+2-length(siz))]; % Leading singletons 
    b = reshape(x,s(n+1:end));
  else
    b = permute(x,[n+1:ndims(x) 1:n]);
  end
else  % Shift to the right (padding with singletons).
    b = reshape(x,[ones(1,-n),siz]);
end
nshifts = n;


