function [x,y] = line2arrow(xy,n,l,d,d0,r,varargin)

% ARROWS    Calculates multiple Arrows along a Polygon
%
% [X,Y] = ARROWS( XY , Number , Length , Distance , Offset , Resolution , ...
%                      BandWidth , HeadLength , HeadWidth , HeadFak , Mode );
%
% Number   = Number of Arrows
%
% Offset   = Offset01 | [ OffsetBegin  OffsetEnd ]
%
% Length   = Length   | [ L(1) .. L(Number)   ]
%
% Distance = Distance | [ D(1) .. D(Number-1) ]
%
% The Length of XY is equal to SUM(Offset) + SUM(Length) + SUM(Distance)
%
% BandWidth, HeadLength, HeadWidth, HeadFak are normalized to Length of Arrow
%
% For absolute Units use a nonzeros imaginary Part for:
%     BandWidth, HeadLength, HeadWidth
%
% Mode:  1   Pike at Beginn
%        0   no Pike
%       -1   Head only
%
% see also: ARROWBAND
%

Nin = nargin;

Nout = nargout;

Msg = '';
nl  = char(10);

if Nin < 1
   error('Not enough InputArguments.')
end
if Nin < 2
   n = [];
end
if Nin < 3
   l = [];
end
if Nin < 4
   d = [];
end
if Nin < 5
   d0 = [];
end
if Nin < 6
   r = [];
end

%*******************************************************
% Check XY-Line

s  = size(xy);
ok = ( ( prod(s) == s(1)*s(2) )  &  ( prod(s) >= 4 ) );
if ok
   ok = any( s == 2 );
   if ( ok  &  ( s(1) == 2 )  &  ~( s(2) == 2 ) );
      xy = permute(xy,[2 1]);
   end
end

if ~ok
     Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
             'XY must define a Line, 2-Rows or 2-Columns with min. 2 Points.'];
end

%*******************************************************
% Defaults

if isempty(n)
   n = 1;
else
   n = n(1);
end

if isempty(l)
   l = 1*ones(1,n);
else
   l = l(:)';
  nl = size(l,2);
  ok = any( nl == [ 1  n ] );
  if ~ok
     Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
             'Length must be a Scalar or Vector with length of Number.' ];
  elseif ( nl == 1 )
    l = l(ones(1,n));
  end
end

if isempty(d)
   d = 1/2*ones(1,n-1);
else
   d = d(:)';
  nd = size(d,2);
  ok = any( nd == [ 1  n-1 ] );
  if ~ok
     Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
             'Distance must be a Scalar or Vector with length of Number-1.' ];
  elseif ( nd == 1 )
    d = d(ones(1,n-1));
  end
end

if isempty(d0)
   d0 = [ 0  0 ];
else
   d0 = d0(:)';
   if size(d0,2) == 1
      d0 = [ d0  d0 ];
   else
      d0 = d0([1 2]);
   end
end

if isempty(r)
   r = 100*ones(1,n);
else
   r = r(:)';
  nr = size(r,2);
  ok = any( nr == [ 1  n ] );
  if ~ok
     Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
             'Resolution must be a Scalar or Vector with length of Number.' ];
  elseif nr == 1
     r = r(ones(1,n));
  end
  r = ceil(max(2,abs(r)));
end

%*******************************************************
if ~isempty(Msg)
    error(Msg)
end
%*******************************************************

%  Vector in Length for ArrowSegments

m  = max(r);

ii = ones(m,1) * (l./(r-1));  % [ Start1 Start2 ... StartN ]

ii(1,2:n) = l( 1 : n-1 ) + d;
ii(1,1)   = d0(1);
ii(1,:)   = cumsum(ii(1,:),2);

ii = cumsum(ii,1);

ii = reshape(ii,m*n,1);

%--------------------------------------------

w = sum(l) + sum(d)*(n>1) + sum(d0);

ii = ii / w;   % [ 0 .. 1 ]

ii(1) = ii(1) + 1e-10;
ii(m) = ii(m) - 1e-10;

dxy = sqrt(sum(diff(xy,1,1).^2,2));

jj = ( dxy == 0 );
if any(jj)
   if all(jj)
      error('Singular Point in XY.');
   else
      warning('Remove Duplicate Points.')
   end
       jj = find(jj);
   dxy(jj)     = [];
    xy(jj+1,:) = [];
end

dxy = cat( 1 , 0 , dxy );

dxy = cumsum(dxy);

xy = interp1( dxy , xy , ii*dxy(end) );

xy = reshape(xy,m,n,2);

xy = permute(xy,[1 3 2]);  % [ m by XY by n ]

x = NaN*zeros(1,n);
y = NaN*zeros(1,n);

for ii = 1 : n

    z = arrowband(xy(1:r(ii),:,ii),varargin{:});

    nz = size(z,1);
    nn = size(x,1);

    if nz > nn
       x = cat( 1 , x , NaN*ones(nz-nn,n) );
       y = cat( 1 , y , NaN*ones(nz-nn,n) );
    end

    x(1:nz,ii) = z(:,1);
    y(1:nz,ii) = z(:,2);

end

%--------------------------------------------
% Check for NaN's, fill End

ii = ( isnan(x) | isnan(y) );
if ~any(ii)
    if ( Nout == 1 ) & ( n == 1 )
       x = [ x  y ];
    end
    return
end

  ii  = find(ii);
x(ii) = NaN;
y(ii) = NaN;

nn = size(x,1);

if ~any(isnan(x(nn,:)))
    if ( Nout == 1 ) & ( n == 1 )
       x = [ x  y ];
    end
    return
end

% Last not NaN

jj = nn - sum(cumprod(double(isnan(x(nn:-1:1,:))),1),1);

ok = ( ( 0 < jj ) & ( jj < nn ) );

if any(ok)

   ok = find(ok);

   for ii = ok(:)'

       x((jj(ii)+1):nn,ii) = x(jj(ii),ii);
       y((jj(ii)+1):nn,ii) = y(jj(ii),ii);

   end

end

if ( Nout == 1 ) & ( n == 1 )
   x = [ x  y ];
end

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [xy] = arrowband(xy,bw,hl,hw,hf,pike)

% function XY = arrowband( XY , BandWidth , ...
%                       HeadLength,HeadWidth,HeadFak,Pike)
% 
%   XY = [ X  Y ]  Coordinates of MiddleLine
%
%    BandWidth, HeadLength, HeadWidth       normalized |  absolut if imag
%   
%    Pike = PikeMode:  1   Pike at Beginn
%                      0   no Pike
%                     -1   Head only
%
%   C0 = [Lon0 Lat0]  Coordinates of  LeftLine
%   C1 = [Lon1 Lat1]  Coordinates of RightLine
%  

Nin = nargin;

if Nin < 1
   error('Not enough InputArguments.')
end
if Nin < 2
 bw = [];
end
if Nin < 3
 hl = [];
end
if Nin < 4
 hw = [];
end
if Nin < 5
 hf = [];
end
if Nin < 6
 pike = [];
end

%*******************************************************
% Check XY-Line

s  = size(xy);
ok = ( ( prod(s) == s(1)*s(2) )  &  ( prod(s) >= 4 ) );
if ok
   ok = any( s == 2 );
   flip = ( ok  &  ( s(1) == 2 )  &  ~( s(2) == 2 ) );
   if flip
      xy = permute(xy,[2 1]);
   end
end

if ~ok
     error('XY must define a Line, 2-Rows or 2-Columns with min. 2 Points.');
end
%*************************************
% Defaults

if isempty(bw)
 bw = 0.1;
end

if isempty(hl)
 hl = 1/3;
end

if isempty(hw)
 hw = 3*bw;
end

if isempty(hf)
 hf = 2/3;
end

if isempty(pike)
   pike = 1;
else
   pike = sign(pike);
end

%*************************************

band = ~( bw == 0 );

head = ~( hl == 0 );

%*************************************

n = size(xy,1);

dd = diff(xy,1,1);

ww = cat( 1 , 0 , sqrt(sum(dd.^2,2)) );

ww = cumsum(ww);

ll = ww(n);

ii = ( 1 : n-2 );

dd = cat( 1 , dd(1,:) , (dd(ii,:)+dd(ii+1,:))/2 , dd(n-1,:) ); % Gradient

dd = dd ./ ( sqrt(sum(dd.^2,2)) * [ 1  1 ] ); % normalized

dd = dd(:,[2 1]);        % [  dy  dx ];

dd(:,1) = -dd(:,1);      % [ -dy  dx ];

%*************************************
% Scale

if imag(bw) == 0;
   bw = bw * ll;
else
   bw = real(bw);
end

bw = abs(bw);

if imag(hl) == 0
   hl = hl * ll;
else
   hl = real(hl);
end

hl = abs(hl);


if imag(hw) == 0
   hw = hw * ll;
else
   hw = real(hw);
end

hw = abs(hw);


%--------------------------------

bw = bw / 2;  % !!!
hw = hw / 2;  % !!!

hw = max( hw , bw );
hf = max( hf , bw/hw );

hf = min(hf,1-(1e-10));

%*********************************************************
% Initialize

xy0 = ones(0,2,2);

xyp = xy0;  % Pike
xyh = xy0;  % Head

%*********************************************************
% Pike 

if ( pike == 1 )  &  band  &  head

   xyp = getline( xy , ww , dd , bw*hl/hw , [0 bw] , 1 );

   pike = ~isempty(xyp);

end

%*********************************************************
% Head

if head

   [ xyh , ii ] = getline( xy , ww , dd , -(ll-hl) , [hw 0] , 0 );

   head = ~isempty(xyh);

end

%-------------------------------
% HeadBase

base = head;

if base

   hb = bw * ( 1 - ( pike == -1 ) );

   [xyb,jj] = getline( xy(ii,:) , ww(ii,:) , dd(ii,:) , ll-hf*hl , [hb hw] , 1 );

   base =  ~isempty(xyb);

   if base

      xyh = cat( 1 , xyb , xyh );

       ii = ( 1 : max(ii(jj))-1 );

       xy = xy(ii,:);
       dd = dd(ii,:);

   end

end

%-------------------------------
% Base at Head

if ~base & band & head

   ii = find( ww < ll-bw*hl/hw );

   xy = xy(ii,:);
   dd = dd(ii,:);

end

%*********************************************************
% Band

if band

%   if isempty(xy)

%      xy = xy0;

%   else

      xy = xy(:,:,[1 1]);

      xy(:,:,1) = xy(:,:,1) + bw * dd;
      xy(:,:,2) = xy(:,:,2) - bw * dd;

%   end

end

%--------------------------------
% Pike + Band + Head

if pike == -1
   xy = xyh;
else
   xy = cat( 1 , xyp , xy , xyh );
end

xy = cat( 1 , xy(end:-1:1,:,1) , xy(:,:,2) );

if flip
   xy = permute(xy,[2 1]);
end

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [xy,ii] = getline(xy,ww,dd,p,l,f)

% [XY,II] = GETLINE( XY , Way , D , Parameter , Limit , Flip );
%

ii = find( sign(p)*ww <= p  );

if isempty(ii)
   xy = ones(0,2,2);
   return
end

ni = size(ii,1);

if f
   ii = ii(ni:-1:1);
end

ww = linspace(l(1),l(2),ni)';

xy = xy(ii,:,[1 1]);

xy(:,:,1) = xy(:,:,1) + ww(:,[1 1]) .* dd(ii,:);
xy(:,:,2) = xy(:,:,2) - ww(:,[1 1]) .* dd(ii,:);

