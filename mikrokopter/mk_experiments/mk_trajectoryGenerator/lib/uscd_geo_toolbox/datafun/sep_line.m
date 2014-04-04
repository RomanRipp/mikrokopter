function [x1,y1] = sep_line(p,x0,x,y);

% SEP_LINE  Transforms and Separates Data by NaN
%
%  [ XT , YT ] = SEP_LINE( Period , XOffset , X , Y )
%
% Transform X, using CENTER( Period , XOffset , X )
%
% Fill Gabs in X ( differences larger Period/2 )
%  with NaN's, same for corresponding Y
%
%  X = [ N by 1 ]
%  Y = [ N by M ]
%


Nin  = nargin;
Nout = nargout;


if Nin < 3
 error('Not enough Input Arguments.');
end
 
if Nin < 4
  y1 = [];
else
  y1 = y;
end

if isempty(x)
  return
end


if ~( isnumeric(x)  &  ( prod(size(x)) == size(x,1) ) )
  error('X must be a numeric ColumnVector.')
end


if ~isempty(y);
  if ~( isnumeric(y)  &  ( size(y,1) == size(x,1) ) & ...
        (  prod(size(y)) == size(y,1)*size(y,2) ) )
    error([ 'Y must be a 2-dimensional Matrice '   ...
            'with the same Number of Rows like X.'     ]) 
  end
end


%------------------------------------
% Transform Data

x = center(p,x0,x);


%------------------------------------
% Find breaks

dx = diff(x);

ib = ( abs(dx) > p/2 );

if ~any(ib)
  x1 = x;
  return
end

ib = find(ib);     % Index before Border


%-------------------------------------
% X-Borders

xlim = [ 0   p ] + x0;

% BorderValue

dx = dx(ib);

xb1 = xlim( ( 1 + ( dx < 0 ) ) );
xb2 = xlim( ( 2 - ( dx < 0 ) ) );


%-------------------------------------

sb = size(ib,1);   % Number of BorderHits

s0 = size(x,1);

s1 = s0 + 3 * sb;  % Size of new Vector


%-------------------------------------
% Index for x in new Vector

i1       = ones(s0,1);
i1(ib+1) = i1(ib+1) + 3;
i1       = cumsum(i1,1);

%-------------------------------------
% Index for ib in new Vector

ib = ib + cumsum( 3*ones(sb,1) , 1 ) - 3;

% Index on NaN-Separator in new Vector
ib = ib + 2; 

%-------------------------------------
% New Vector

x1       = NaN * zeros(s1,1);
x1(i1)   = x;

x1(ib-1) = xb1;
x1(ib+1) = xb2;

if isempty(y)  |  ( Nout < 2 ) 
  return
end

%-------------------------------------
% Index, surrounding Border

i01 = ones(sb,1)*[ -2  2 ] + ib(:,[1 1]);

%-------------------------------------
% Center New

x00 = x0 + 180;

x01 = center(p,x00,x1(i01));

if ( sb == 1 )  &  ( size(x01,1) == 2 ) 
  x01 = permute(x01,[2 1]);
end


% BorderValue for new Center

xb = x00+p/2;

sy = size(y,2);

y1        = NaN * zeros(s1,sy);
y1(i1,:)  = y;


% linear Weight at xb

x01 = ( x01(:,2) - xb ) ./ ( x01(:,2) - x01(:,1) );

y1(ib-1,:) = y1(i01(:,2),:) - ...
             ( y1(i01(:,2),:) - y1(i01(:,1),:) ) .* ( x01*ones(1,sy) );

y1(ib+1,:) = y1(ib-1,:);

 


