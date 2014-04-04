function [ii,jj,kk] = grid2tri(si,mode)

% GRID2TRI  Triangulation of Grid
%
% [ I1 , I2 , IND ] = GRID2TRI( SIZ , Mode )
%
%             IND   = GRID2TRI( SIZ , Mode )
%
%  SIZ defines the Size of an 2-dimensional Grid
%  SIZ = [ NRow  NColumn ]
%
% A Grid defined by SIZ has PROD( SIZ - 1 ) 
%  rectangular Faces. Each of them is devided into 2
%   triangular Faces. The Indices of the Vertices
%  of the triangular Faces are returned.
%
%  The Outputs have 3 Rows for the Vertices of the Triangular
%   and 2 * PROD( SIZ - 1 ) Columns for the Faces.
%
%  I1 and I2 are Subscripts for RowIndex and ColumnIndex
%   of the Vertices.
%
%  IND is the linear Index of the Vertices.
%  IND = I1 + ( I2 - 1 ) * SIZ(1)
%
%------------------------------------------------------------
%
% Mode is a single (imaginary) Number to define the Order
%  and the Shape of the Triangles, following the Columns
%  and switching at each Row.
%
% Mode = Swap + Shape * i
%
%  Swap defines the Mode to change the Direction of the Columns:
%
%  Swap = 0  no switch
%         1  switch at the even     Rows
%        -1  switch at the oddinary Rows
%
%         0: [1,1] .. [1,N] , [2,1] .. [2,N] , ... , [M,1] .. [2,N]
%         1: [1,1] .. [1,N] , [2,N] .. [2,1] , ...
%        -1: [1,N] .. [1,1] , [2,1] .. [2,N] , ...
% 
%        where [ M  N ] is the Size of the rectangular Faces:
%              [ M  N ] = SIZ - 1;  M = NRow - 1;  N = NColumn - 1;
%
%  Shape defined the Triangles
%
%  Shape = 0  "V" - Shape
%          1  "\" - Shape
%         -1  "/" - Shape
%
%  default: Mode = 1;  switch at the even Rows, "V"-Shape
%
%*******************************************************
%
% Example:
%
%%*****************************
%% Matlab Logo
%
%  zz = membrane(1,8);
%  si = size(zz);
%
%  x = linspace(-1,1,si(2));
%  y = linspace(-1,1,si(1));
%
%  xx = sign(x) .* (abs(x).^1.5);
%  yy = sign(y) .* (abs(y).^1.5);
%
%  [xx,yy] = meshgrid(xx,yy); % NonRegular Grid
%
%  % Triangulation
%  [i1,i2,ind] = grid2tri(si);
%
%%-----------------------------
%% One Patch per Face
%
%  figure, hold on, view([-37 30]);
%
%  patch( 'xdata' ,  x(i2)  , ...
%         'ydata' ,  y(i1)  , ...
%         'zdata' , zz(ind) , ...
%     'edgecolor' , 'k'    , ...
%     'facecolor' , 'r'          );
%
%%-----------------------------
%% Vertices and Faces
%
% figure, hold on, view([-37 30])
%
% patch( 'vertices'        , [ xx(:) yy(:) zz(:) ] , ...
%        'faces'           ,   ind'   , ...
%        'facevertexcdata' ,  zz(:)   , ...
%        'cdatamapping'    , 'scaled' , ...
%        'edgecolor'       , 'k'      , ...
%        'facecolor'       , 'interp'         );
%
%%*****************************
%% Mode Demonstration
%
% si = [ 3  5 ];
%
% x = ( 1 : si(2) );  x = x(ones(1,si(1)),:);
% y = ( 1 : si(1) )'; y = y(:,ones(1,si(2)));
%
% nc = 2 * prod(si-1); cm = hsv(ceil(1.1*nc));
% cc = ( 1 : nc ); cc = cc(ones(1,3),:);
%
% for swp = [ 0  1  -1 ], for shp = [ 0 1 -1 ]
%
%     ind = grid2tri(si,swp+shp*i);
%
%     figure, colormap(cm(1:nc,:)); hold on
%     axis([ 1 si(2)  1 si(1) ]), box on, axis ij
%     set(gca,'tickdir','out','xtick',1:si(2),'ytick',1:si(1));
%
%     patch( 'xdata' ,  x(ind)  , ...
%            'ydata' ,  y(ind)  , ...
%            'zdata' ,  cc      , ...
%            'cdata' ,  cc      , ...
%     'cdatamapping' , 'direct' , ...
%        'edgecolor' , 'interp' , ...
%        'facecolor' , 'none'   , ...
%       'linestyle'  , '--'     , ...
%       'linewidth'  ,  1       , ...
%       'marker'     , 'o'      , ...
%       'markersize' ,  8       , ...
%  'markeredgecolor' , 'k'      , ...
%  'markerfacecolor' , 'flat'            );
%
%     s1 = ''; if ~( swp == 0 ), s1 = char([44-swp  49]); end
%     s2 = ''; if ~( shp == 0 ), s2 = char([44-shp 105]); end
%     if isempty(s1) & isempty(s2), s1 = '0'; end
%
%     title(sprintf('Mode = %s%s',s1,s2))
%
%     cb = colorbar('horiz');
%     set(cb,'xlim',[1 nc]+0.5*[-1 1]);
%     set(findobj(cb,'type','image'),'xdata',[1 nc]);
%
% end, end
%
%

Nout = nargout;

ok = ( isnumeric(si) & ( prod(size(si)) == 2 ) );
if ok
   ok = all( isfinite(si) & ( si > 1 ) & ( mod(si,1) == 0 ) ); 
end

if ~ok
    error('Size must be a 2-Element Vector with Integers larger 1.');
end

if nargin < 2
   mode = 1;
else
   if ~( isnumeric(mode) & ( prod(size(mode)) == 1 ) )
       error('Mode must be a single Numeric.');
   end
end

mode = [ real(mode)  imag(mode) ];
mode = sign(mode);

swap = mode(1);
mode = mode(2);

%----------------------------------------
% Indices, run over Columns

si = ( si(:)' - 1 ) .* [ 1  2 ];

ii = ( 1 : prod(si) );

jj = ii - si(2) * floor( ii / si(2) );
jj = jj + si(2) * ( jj == 0 );
jj = ceil( jj / 2 );

ii = ceil( ii / si(2) );

%----------------------------------------
% Swap JJ if even/odd II ans NonZero Mode(1)

if swap == 0
   im = ones(size(ii));     % Allways Odd II
else
   im = ii - 2 * floor( ii / 2 );           % True for Odd II
   jm = ( 1 + swap ) / 2 - swap * im;       % IM | ~IM
   jj = jj + ( 1 + si(2)/2 - 2*jj ) .* jm;
end

if mode == 0
   jm = jj - 2 * floor( jj / 2 );  % True for Odd JJ
else 
   % Mode:    1  ==>  JM ==  IM
   %         -1  ==>  JM == ~IM
   jm = ( 1 - mode ) / 2 + mode * im;
end

%****************************************
% Three Points

ii = ii([1 1 1],:);
jj = jj([1 1 1],:);

%----------------------------------------
% 1. Point

ii(1,:) = ii(1,:) + ( 1 - jm );
jj(1,:) = jj(1,:) + ( 1 - im );

%----------------------------------------
% 2. Point

kk = ( 1 : 2 : si(1)*si(2)-1 );

ii(2,kk) = ii(2,kk) + jm(kk);
jj(2,kk) = jj(2,kk) + ( 1 - im(kk) );

ii(2,kk+1) = ii(2,kk+1) + ( 1 - jm(kk+1) );
jj(2,kk+1) = jj(2,kk+1) + im(kk+1);

%----------------------------------------
% 3. Point

ii(3,:) = ii(3,:) + jm;
jj(3,:) = jj(3,:) + im;


%****************************************

if     Nout < 2
       ii = ii + ( jj - 1 ) * (si(1)+1);
elseif Nout > 2
       kk = ii + ( jj - 1 ) * (si(1)+1);
end

