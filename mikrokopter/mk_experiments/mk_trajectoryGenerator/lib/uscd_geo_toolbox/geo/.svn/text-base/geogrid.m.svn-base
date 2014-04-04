function [gx,gy,grd,xyg] = geogrid(area,int,res)

% GEOGRID  Computes a Grid for geodetic Coordinates
%
%    GEOGRID( AREA ,   INT )
%    GEOGRID( AREA , -NINT )
%
% AREA = [ LonMin LonMax LatMin LatMax ]
%
%  INT   Intervall for Grid, 
%        the imaginary Part defines the Startvalue
%
% NINT   Number of Intervalls in Grid, default: 6
%
% The Imaginary part of the second Input defines
%  the Resolution or Spacing:
%
%    GEOGRID( AREA ,   INT ,  RES )
%    GEOGRID( AREA , -NINT ,  RES )
%
%    GEOGRID( AREA ,   INT , -SPC )
%    GEOGRID( AREA , -NINT , -SPC )
%
%  RES   Resolution for Grid, default: 100
%  SPC   Spacing for Grid,
%        the imaginary Part defines the Startvalue
%
%    GEOGRID( ... , [] ,  -SPC ) 
% uses the default value SPC for INT, equal to 
%    GEOGRID( ... , SPC , -SPC )
%
%------------------------------------------------------
% Outputs:
%
% XY = GEOGRID( ... )  XY-Coordinates of Grid
%
% [ GX , GY , GRD , XYG ] = GEOGRID( ... )
%
% GX = XYmeridional (constant Lon)
% GY = XYzonal      (constant Lat)
% 
% GRD  Grid by INT or NINT
% XYG  Grid by RES or SPC
%
% The XY-Cooridates of the meridional and zonal 
% GridLines are seperated by NaN-Rows
%
%------------------------------------------------------
%
% GEOGRID( AXE , ... ) uses the AxesLimits ([XLim YLim]
%  specified by AxesHandle AXE for AREA. 
%  If no Output is requested, the Grid is 
%  plotted into AXE with dotted Lines, 
%  the created LineHandle has the Tag GEOGRID
%

if nargin < 2
   int = [];
end

if nargin < 3
  res = [];
end

Nout = nargout;

%********************************************************
% Check Inputs

msg = {};

axe = [];

s = size(area);
p = prod(s);

ok = ( isnumeric(area) & any(p==[1 4]) );
if ok
   ok - all(isfinite(area));
end

if ~ok
    msg = cat(1,msg,{'First Input must be a finite numeric (4-Element Area or AxesHandle).'});
else
    if ( p == 1 )
       ok = ishandle(area);
       if ok
          ok = strcmp(get(area,'type'),'axes');
       end
       if ~ok
           msg = cat(1,msg,{'Single first Input must be an AxesHandle.'});
       else
           axe = area;
           area = cat(1,get(area,'xlim'),get(area,'ylim'));
       end
    else
       area = area(:)';
       while area(2) <= area(1)
             area(2) = area(2) + 360;
       end
       area([3 4]) = area([3 4]+[1 -1]*(area(3)>area(4)));
       area = cat( 1 , area([1 2]) , area([3 4]) );
    end
end

%--------------------------------------------------------

if ~isempty(int)
    ok = ( isnumeric(int) & ( prod(size(int)) == 1 ) );
    if ok
       ok = isfinite(int);
    end
    if ~ok
        msg = cat(1,msg,{'Second Input must be a single finite numeric.'});
    end
end

if ~isempty(res)
    ok = ( isnumeric(res) & ( prod(size(res)) == 1 ) );
    if ok
       ok = isfinite(res);
    end
    if ~ok
        msg = cat(1,msg,{'Third Input must be a single finite numeric.'});
    end
end

%--------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Invalid Inputs.\n%s',msg))
end


%********************************************************

if isempty(int)
   int = 0;
end

if isempty(res)
   res = 100;
end

idv = imag(int);
int = real(int);

rdv = imag(res);
res = real(res);


if int == 0
   if res < 0
      int = abs(res);
   else
      int = -6;
   end
end

%--------------------------------------------------------

[grd,int] = getgrid(area,int,idv);

%--------------------------------------------------------

if res > 0
   xyg = cell(1,2);
   for ii = 1 : 2
       xyg{ii} = linspace(area(ii,1),area(ii,2),res);
   end
else
    [xyg,res] = getgrid(area,-res,rdv);
end

   for ii = 1 : 2
       xyg{ii} = xyg{ii}(:);  % ColumnVectors
   end

%--------------------------------------------------------

for ii = [ 1  2 ]

    jj = 3 - ii;   % Opposite

    m = prod(size(grd{ii}));
    n = prod(size(xyg{jj}));

    if ( m == 0 ) | ( n == 0 )
       grd{ii} = [ NaN  NaN ];
    else
       n = n + 1;  % Insert NaN
       g = cumsum(ones(m*n,2),1);
       g(:,2) = mod(g(:,2),n);
       kk = find( g(:,2) == 0 ); % Insert NaN here
       g(kk,2) = 1;
       g(:,1) = ceil((g(:,1))/n);
       g(:,1) = grd{ii}(g(:,1));
       g(:,2) = xyg{jj}(g(:,2));
       grd{ii} = g(:,[ii jj]);
       grd{ii}(kk,:) = NaN;
    end

end

%********************************************************

if Nout <= 1
   gx = cat(1,grd{:});
else
   gx = grd{1};
   gy = grd{2};
end

if ( Nout == 0 )
   if ~isempty(axe)
       h = line('xdata',gx(:,1),'ydata',gx(:,2), ...
                'linestyle',':','marker','none', ...
                'parent',axe,'tag','GEOGRID');
       if res < 0
          set( h , 'linestyle' , 'none' , 'marker' , '+' );
       end
       clear gx
   end
end

%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [grd,int] = getgrid(lims,int,dev)


if int < 0

   int = abs(int);

   tick_int = [ 90 ( [ 60 30 20 10  5  2  1 ]/1         ) , ...    % deg
                   ( [    30 20 10  5  2  1 ]/60        ) , ...    % min
                   ( [    30 15 12  6       ]/3600      ) , ...    % min/10
                   ( [    30 15 12  6       ]/3600/10   ) , ...    % min/100
                   ( [    30    12  6       ]/3600/100  )       ]; % min/1000


   nt = size(tick_int,2);

   dl = diff(lims,1,2);

   ii = abs( dl(:,ones(1,nt)) / int - tick_int([1 1],:) );

   [ hilf , ii ] = min(ii,[],2);

   [ hilf , kk ] = max(tick_int(ii));

    ok = find( ( rem( tick_int , tick_int(ii(kk)) ) == 0 )   &  ...
                  (  tick_int >= tick_int(ii(kk))  ) ); 

   [hilf,jj] = min( abs( dl(3-kk) / int - tick_int(ok)) );

   ii(3-kk) = ok(jj);

   int = tick_int(ii(1));

   mdl = 0;  % Don''t check with Modulo !!!

end


  grd = cell(1,2);

  dev = mod(dev,int);

  lims = lims - dev;

  limt = lims/int;
  limr = round(limt);
    jj = ( abs(limt-limr) < 1e-10 );
  if any(jj)
     jj = find(jj);
     limt(jj) = limr(jj);
  end


  for ii = 1 : 2

      grd{ii} = int * ( ceil(limt(ii,1)) : floor(limt(ii,2)) ) + dev;

  end
