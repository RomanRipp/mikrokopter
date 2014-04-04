function x = sampline(x,varargin)

% SAMPLINE   Smooth and lowsample Vectors
%
% Smooth and lowsamples Vectors or 2-dimensional  Matrices, 
%   build from them.
% The Vectors can be separated into Segments by NaN-Lines.
%   Each Segment will smoothed and lowsampled separatly.
% 
% X = SAMPLINE ( X , Window )
%
% Window = Sample + i * Smooth
%
%  defaults: Sample = 3
%            Smooth = 5
%
% If the Value of Smooth is ZERO, it will calculated like:
%
%    Smooth = 2 * ceil( Sample / 2 ) + 1
%
% If the Value of Sample is ZERO, it will calculated like:
%
%    Sample = Smooth - 2;
%
% Note: Smooth should be allways an oddinary Number, 
%       if not, the Value will raised by ONE:
%
%    Smooth =  2 * floor( Smooth / 2 ) + 1
%
%-------------------------------------------------------
%
% Optional Inputs:
%
%  SAMPLINE( ... , Mode , Accuracy , DIM )
%
%    Mode gives Type of SmoothWindow, using MEANIND1.
%
%    Mode = 'linear' | 'binomial' | {'cosine'} | 'gauss'
%
%    DIM  works along Dimension DIM, DIM = 1 | 2
%          (first real integer, following "Window")
%         default: DIM == longest Dimension of X
%
%    Accuracy gives the Value to detect closed Segments,
%       in case of MatriceInput X.
%      Closed Segments will closed for smoothing.
%            
%-------------------------------------------------------
%
% see also: MEANIND1 (required), WINDOW, IND2GRP, GRP2IND
%
%-------------------------------------------------------
%
% Example:
%
%   load coast % Mapping Toolbox
%
%   x = sampline([long lat],5);
%
%   figure, hold on
%
%   plot(long,lat);
%   plot(x(:,1),x(:,2),'r.-');
%

Nin = nargin;

if Nin < 1
   error('Input XY is missing.');
end

si = size(x);
ps = prod(si);

if ~( ps == si(1)*si(2) )
    error('XY must be a 2-dimensional Matrice.');
end

%***************************************************
% Get Inputs

[smp,smt,mode,acc,dim] = checkin(si,varargin{:});
%%% [smp smt], acc
if ~any( [ smp  smt ] > 1 )  % No Smooth and Sampling
    return
end

if si(dim) < 2
   return
end

flip = ( dim == 2 );

if flip 
   x  = permute(x,[2 1]);
   si = si([2 1]);
end

sr = si(2);     % Length of Seperator-Row (NaN)

%********************************************************
% Get Segment, separated by NaN-Rows

ok = double( sum(isnan(x),2) == sr );   % Ok for NaN-Row

if all(ok)               % NaN's only
   x = x(1,:);
   if flip
      x = permute(x,[2 1]);
   end
   return
end

ii = find(~ok);

i0 = cat( 1 , 1 , find( diff(ii,1,1) > 1 )+1 , size(ii,1)+1 );

ls = diff(i0,1,1);     % Length

i0 = ii(i0(1:end-1));  % StartIndex

%********************************************************
% Check Segments

% New Length after LowSampling

ln = ceil(ls/smp);

% Check for Closed Segments

if sr == 1
   cl = zeros(size(i0));
else
   cl = x(i0,:) - x(i0+ls-1,:);
   cl = ( abs(cl)   <= acc );
   cl = ( sum(cl,2) == sr  );
end

% Single Point Segments !!!

sok = ( ( ln == 1 ) & ~cl );

if any(sok)
   ii = find(sok);
   ok(i0(ii))          = -1;
   ok(i0(ii)+ls(ii)-1) = -1;
end

% Find Good Segments !!!

sok = ( ( ln > 3*cl ) & ( ~cl | ( ls >= smt ) ) );

if ~any(sok)
    x = NaN*zeros(any(ok),sr);
    if flip
       x = permute(x,[2 1]);
    end
    return
end

% Remove Bad Segments

if ~all(sok)

    ii = find(~sok);
    ii = ind2grp(i0(ii),ls(ii),1);

    ok(ii) = 0;

    ii = find(sok);

    i0 = i0(ii);
    ls = ls(ii);
    ln = ln(ii);
    cl = cl(ii);

end

%********************************************************
% Smooth Segments

if smt > 1 

   s2 = ( smt - 1 ) / 2;  % Half SmoothWindow

   ns = size(i0,1);

   for ii = 1 : ns

       jj = i0(ii) - 1 + ( 1 : ls(ii) );

       if cl(ii)
          % !!! x(1,:) == x(ls,:) !!!
          kk = ( 1 : s2 );
          kk = cat( 2 , jj(ls(ii)-s2+kk-1) , jj , jj(kk+1) );
           m = meanind1(x(kk,:),smt,mode);
          x(jj,:) = m(s2+(1:ls(ii)),:);    
       else       
          x(jj,:) = meanind1(x(jj,:),smt,mode);
       end

   end

end

%********************************************************
% Get Segments

ok = -ok;   % "-1" for NaN-Row !!!

ok(ind2grp(i0,ln,smp)) = 1;

% ok =  1  Good Data
%      -1  NaN-Row
%       0  Bad Data
%

%----------------------------------
% Use the EndPoint too

ok(i0+ls-1) = 1;


ii = find( ~( ok == 0 ) );

x  = x(ii,:);
ok = ok(ii,:);

% Remove Duplicate NaN-Rows

ok = ( ok == -1 );
if sum(ok) > 1
   ok = find(ok);
   ok = ok( find( diff(ok,1,1) == 1 ) + 1 );
   x(ok,:) = [];
end

% Flip back

if flip 
   x = permute(x,[2 1]);
end

%********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ii = ind2grp(i0,l,s);

n = size(i0,1);

%  l = ceil(l/s);   !!! Allready done !!!

ii = s * ones(sum(l),1);
jj = cumsum( cat(1,1,l) , 1 );

ii(jj(1:n)) = i0;

if n > 1
   ii(jj(2:n)) = ii(jj(2:n))-(i0(1:n-1,1)+s*(l(1:n-1)-1));
end

ii = cumsum(ii,1);

%********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [smp,smt,mode,acc,dim] = checkin(si,varargin);

n = nargin - 1;

def = [ 3  5 ];    % default: [ smp  smt ]

int = [];
acc = [];
dim = [];

mode = 'c';

%****************************
% Get Inputs from varargin

for ii = 1 : n

    v = varargin{ii};

    s = size(v);
    p = prod(s);

    if ( strcmp(class(v),'char') & ( p == s(2) ) & ~( p == 0 ) )

       mode = v;

    else

       ok = ( isnumeric(v) & ( p == 1 ) );
       if ok
          vv = [ real(v)  imag(v) ];
          ok = all( ( vv >= 0 ) & isfinite(vv) );
       end

       if ~ok
           error('Inputs must be single positive finite numerics.')
       end

       rv = ( round(vv) == vv );  % Integer
       ir = ( vv(2) == 0 );       % Real

       if     isempty(int) & all(rv)
          int = vv;
       elseif isempty(dim) & ir & any( v(1) == [ 1  2 ] )
          dim = v;
       elseif isempty(acc) & ir
          acc = v;
       end

    end

end

%****************************

if isempty(dim)
   dim = 1 + ( si(2) > si(1) );
end

if isempty(acc)
   acc = 1e-10;
end

%****************************
% Check Window

if isempty(int)
   int = [ 0  0 ];
end

if all(int==1) | all(int==0)
   int = int + def .* ( int == 0 );
   smp = int(1);
   smt = int(2);
   return
end

smp = int(1);
smt = int(2);

if smt == 0
   smt = 2 * ceil( smp / 2 ) + 1;
else
   smt = 2 * floor( smt / 2 ) + 1;
   if smp == 0
      smp = smt - 2;
   end
end
