function [ini,xy] = isoarea(x,y,z,v)

% ISOAREA  returns Area of IsoPatches
%
%  [ INI , XY ] = ISOAREA( X , Y , Z , V )
%
%  returns the Area and CenterPoint of the IsoPatches
%   in Z with the IsoValues of the Vector V.
%
%   INI is a 5-Column Matrice for N IsoPatches:
%
%       [ Area XMena YMean Value Lenght ]
%
%   XY is a CellArray which Elements are
%      the X-Y-Coordinates of the IsoPatches
%      in a 2-Column-Matrice: [ X  Y ]
%
%
%  [ INI , XY ] = ISOAREA( Z , V )
%
%  where Z is an [ M x N ] - Matrice use
%   ( 1 : N ) for X and ( 1 : M ) for Y
%
% 
% see also: CONTOURS, POLYAREA
%
%------------------------------------------------------------
% Example: Area's of horizontal ConeSlices
%
%%------ Slices  --------------------------------------------
%
%   x = ( -1 : 0.01 : 1 );
%   y = x';
%
%   n = size(x,2);
%   z = sqrt( x(ones(1,n),:).^2 + y(:,ones(1,n)).^2 );
%
%   v = ( 0 : 0.1 : 1 );
%
%   [ini,xy] = isoarea(x,y,z,v);
%
%%------ 2D -------------------------------------------------
%
%   figure, box on, grid on, hold on
%
%   plot(ini(:,4),ini(:,1),'b-','linewidth',7);
%   plot(ini(:,4),pi*ini(:,4).^2,'r--','linewidth',3);
%
%   xlabel('Radius'); 
%   ylabel('Area:  \Pi \bullet R^2','interpreter','tex');
%
%   legend({'by ISOAREA' 'by Formula'},2)
%
%%------ 3D -------------------------------------------------
%
%   figure('colormap',jet(128),'renderer','zbuffer');
%
%   axes( 'xlim' , [ min(x)    max(x)    ] , ...
%         'ylim' , [ min(y)    max(y)    ] , ...
%         'zlim' , [ min(z(:)) max(z(:)) ] , ...
%         'clim' , [ 0     max(ini(:,1)) ] , ...
%   'dataaspectratio' , [ 1  1 1/sqrt(2) ] , ...
%         'view' , [ 30  15 ], ...
%     'nextplot' , 'add' );
%
%   jj = ( 1 : 20 : n );
%
%   surface( x(jj) ,y(jj) , z(jj,jj)-0.01 , ...
%            'facecolor' , 'none' , ...
%            'edgecolor' , 'k'    , ...
%            'linestyle' , '-');
%
%   for ii = 1 : size(ini,1)
%
%       patch( 'xdata' , xy{ii}(:,1) , ...
%              'ydata' , xy{ii}(:,2) , ...
%              'zdata' , ini(ii,4)*ones(ini(ii,5),1) , ...
%              'cdata' , ini(ii,1) , ...
%              'cdatamapping' , 'scaled' , ...
%              'facecolor' , 'flat'      , ...
%              'edgecolor' , 'k'                 );
%
%    end
%
%    ax = colorbar('vert');
%    set( get(ax,'xlabel') , 'string' , 'Area' );
%

ini = zeros(0,4);
xy  = cell(0,1);

Nin  = nargin;
Nout = nargout;

if ~any( Nin == [ 2  4 ] )
    error('Inputs must be X, Y, Z and V or Z and V.');
end

%--------------------------------------------------------
% Check ContourVector

if Nin == 2
   v = y;
end

if isempty(v)
   return
end

if prod(size(v)) == 1
   v = v([1 1]);
end

%--------------------------------------------------------
% Get ContourLines

try
  if Nin == 2
     cs = contours(x,v);
  else
     cs = contours(x,y,z,v);
  end
catch
  error(sprintf('Error using CONTOURS.\n%s',lasterr));
end

if isempty(cs)
   return
end

%--------------------------------------------------------
% Get StartIndex of Segements
%
% Value  X1 X2 ...
% Length Y1 Y2 ...
%

ns = size(cs,2);
i0 = 1;
nl = 1;

while 1
      i1 = i0(nl) + cs(2,i0(nl)) + 1;
      if i1 > ns
         break
      end
      i0 = cat( 2 , i0 , i1 );
      nl = nl + 1;
end

%--------------------------------------------------------
% Get Segments:  Area, Center etc ...

ini = zeros(nl,4);

ini = zeros(nl,5);

ini(:,[4 5]) = cs([1 2],i0)';  % [ Value Length ]

if Nout == 2
   xy = cell(nl,1);
end

for ii = 1 : nl

    jj = ( 1 : ini(ii,5) ) + i0(ii);

    ini(ii,[2 3]) = mean(cs(:,jj),2)';    % Center

    ini(ii,1) = polyarea(cs(1,jj),cs(2,jj));

    if Nout == 2
       xy{ii} = cs(:,jj)';  % [ X  Y ]
    end

end
