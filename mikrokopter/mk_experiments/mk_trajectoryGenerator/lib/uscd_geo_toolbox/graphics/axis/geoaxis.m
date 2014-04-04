function geoaxis(axe,md,varargin)

% GEOAXIS   Show Lat/Lon-AxisLabels
%
% GEOAXIS( Axe , Mode )
%
% Mode defines the Longitude/Latitude for XY-Axis
%
% Mode is a String with Maximum 2 Characters for X- and Y-Axis.
%
% Valid Characters:
%
% '0': No GeoLabels
%
% 'x': Longitude for X-Axis 
% 'X': Latitude  for X-Axis
%
% 'y': Latitude  for Y-Axis
% 'Y': Longitude for Y-Axis 
%
% defaults: Axe = GCA
%          Mode = 'xy'
%
% see also: MERCATOR, TIMEAXIS
%
 
%***************************************************

lab = { 'W' ['00' char(176)]  'E' ; 'S' 'EQ' 'N' };

%***************************************************

app = mfilename;
fsp = ( app == filesep );
if any(fsp)
    app = app(max(find(fsp))+1:end);
end
app = upper(app);

%***************************************************

Nin = nargin;

if Nin < 1
   axe = [];
end

if Nin < 2
   md = [];
end

%***************************************************
% Check Inputs

msg = cell(0,1);

if chkstr(axe,1) & ~ischar(md)
   mode = axe;
   axe  = md;
   md   = mode;
end

%---------------------------------------------------
% Check Axes

if isempty(axe)
   fig = get(0,'currentfigure');
   if ~isempty(fig)
       axe = get(fig,'currentaxes');
   end
else
   ok = ( isnumeric(axe) & ( prod(size(axe)) == 1 ) );
   if ok
      ok = ishandle(axe);
      if ok
         ok = strcmp(get(axe,'type'),'axes');
      end
   end
   if ~ok
       msg = cat( 1 , msg , {'Invalid AxesHandle.'} );
   end
end

%---------------------------------------------------
% Check Mode

if isempty(md)
   md = 'xy';
elseif ~chkstr(md,1) & ( size(md,2) <= 2 ) 
   msg = cat( 1 , msg , {'Mode must be a String with max. 2 Elements.'} );
end

%---------------------------------------------------

if ~isempty(msg)
    msg = cat(1,'%s\n',msg{:});
    msg = sprintf('Invalid inputs.\n%s',msg);
    error(msg);
end

%---------------------------------------------------

if isempty(axe)
   return
end

%---------------------------------------------------

if Nin >= 3
   if isnumeric(varargin{1}) & isequal(size(varargin{1}),[1 4])
      % Max be AXE_ZOOM
      set( axe , 'xtickmode' , 'auto' , ...
                 'ytickmode' , 'auto'         );
   end
end

%***************************************************
% Check Modes: '0xXyY'

wrn = '';

lm = lower(md);

ok = ( ( lm == '0' ) | ( lm == 'x' ) | ( lm == 'y' ) );
if ~all(ok)
    wrn = 'Invalid Modes.';
    ii = find(~ok);
    md(ii) = '0';
    lm(ii) = '0';
end

if ( ( size(md,2) == 2 ) & ( lm(1) == lm(end) ) & ~( lm(1) == '0' ) );
   lm = lm(end); 
   md = md(end);
   wrn = 'Duplicate Modes.';
end

if ~isempty(wrn)
    warn(app,wrn)
end

ok = ~( md == '0' );

if ~any(ok)
   return
elseif ~all(ok)
   ii = find(ok);
   md = md(ok);
end

%***************************************************

[upp,siz] = ppunit(axe);

siz = siz([3 4]);

xy = 'xy';

for m = md

    ii = find( xy == lower(m) );
    kk = ii + ( 3 - 2*ii ) * ( m == upper(m) );

    int = 5 * siz(ii) / 300;

    int = ( 1 + ( ii - 1 ) / 2 ) * int;

    [t,l] = get_tick(axe,m,xy(kk),int,lab,app);

    if ~isempty(t)
        set( axe , [m 'tick'] , t , [m 'ticklabel'] , l );
    end

end


%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [tic,tcl] = get_tick(axe,xy,mode,int,lab,app)

acc = 1e-10;   % Accuracy

tic = get(axe,[ xy 'tick']);
lim = get(axe,[ xy 'lim' ]);

tc = tic;   % Original

isy  = ( lower(xy) == 'y' );

auto = strcmp( get(axe,[xy 'tickmode']) , 'auto' );

%---------------------------------------------------------------
% Check for Mercator

ext = [];

if isy

   wrn = '';

   try
      ext = mercator(axe,'Info');
   catch
      wrn = sprintf('Error call MERCATOR( Info ).\n%s',lasterr);
   end

   if isempty(wrn) & ~isempty(ext)
      ok = ( isnumeric(ext) & ( prod(size(ext)) == 1 ) );
      if ok
         ok = ( abs(ext-0.5) <= 0.5 );
      end
      if ~ok
          wrn = 'Invalid Return of MERCATOR( Info ).';
          ext = [];
      end
   end

   if ~isempty(wrn)
       warn(wrn);
   end

end

%---------------------------------------------------------------

if auto & ~isempty(tic)
   tic = tic( find( ( lim(1)-acc < tic ) & ( tic < lim(2)+acc ) ) );
end

if ~isempty(ext) & isy
    lim = merc2deg(lim,ext);
    tic = merc2deg(tic,ext);
end

%*****************************************************************
if auto & ~isempty(tic)
%*****************************************************************
% Prepare Ticks

  tick_int = [ 90 ( [ 60 30 20 10  5  2  1 ]/1         ) , ...    % deg
                  ( [    30 20 10  5  2  1 ]/60        ) , ...    % min
                  ( [    30 15 12  6       ]/3600      ) , ...    % min/10
                  ( [    30 15 12  6       ]/3600/10   ) , ...    % min/100
                  ( [    30    12  6       ]/3600/100  )       ]; % min/1000

  ii = sum( tick_int > diff(lim,1,2)/int );

  ii = ii + ( ii == 0 );

  tick_int = tick_int(ii);

  tic = tick_int * ( ceil(lim(1)/tick_int) : ...
                    floor(lim(2)/tick_int)       );

  if isempty(tic) 
     tic = tick_int * (floor(lim(1)/tick_int) : ...
                        ceil(lim(2)/tick_int)       );
     if size(tic,2) == 2
        tic = lim;
     end     
  end

%*****************************************************************
end
%*****************************************************************

if isempty(tic)
   tcl = [];
   return
end

%*****************************************************************
% Prepare TickLabels

if any( abs( tic-round(tic) ) > 1e-6 )
   flag = 'sexges';  % sexages
else
   flag = 'deg';     % degree
end

%----------------------------------------------------

tc = tic;

if strcmp(mode,'x')

   %  -->  [ -180 .. 180 )
   tc = tc - 360 * floor( (tc+180) / 360 );

   tc(end) = tc(end) + ( 180 - tc(end) ) * ( tc(end) == -180 );

else

   while 1
         ii = ( abs(tc) > 90 );
         if ~any(ii)
             break
         end
         ii = find(ii);
         tc(ii) = tc(ii) + ( sign(tc(ii)) * 180 - 2*tc(ii) );
   end

end

%----------------------------------------------------

lab = lab( 1 + strcmp(mode,'y') , : );

pre = ( isy & strcmp(get(axe,'yaxislocation'),'right') );
pre = ( ~isy | pre );

pre = char( 32 * ones(1,pre) );

app = ' ';

   grad =   fix(tc);
   minu =   fix(  60*(tc-grad));
   secu = round(3600*(tc-grad-minu/60));

    minu_secu =  fix(secu/60);
   secu = secu - 60*minu_secu;
   minu = minu + minu_secu;

    grad_minu =  fix(minu/60);
   minu = minu - 60*grad_minu;
   grad = grad + grad_minu;

   dec_secu = [];
   sec_form = '';

  if any(diff(1e3*grad+minu)==0)
     dec_secu = round(1e2*( 60*(tic-grad) - minu)) ;
     sec_form = '.%2.2d';
  end

  if ~isempty(dec_secu)
      if any(diff(1e3*grad+minu+1e-3*dec_secu)==0)
         dec_secu = round(1e3*( 60*(tic-grad) - minu)) ;
         sec_form = '.%3.3d';
       end
  end

  is_dec = ~isempty(dec_secu);


  nt = prod(size(tic));

  tcl = cell(nt,1);
  tcl(:) = { '' };

  for ii = 1 : nt

      if abs(tic(ii)) < acc

         tcl{ii} = sprintf(' %s ',lab{2});

      else

         if strcmp(flag,'deg')

            tcl{ii} = sprintf('%.0f%s',abs(grad(ii)),char(176));

         else

            if is_dec
               tcl{ii} = sprintf(sec_form,abs(dec_secu(ii)));
            end

            tcl{ii} = sprintf( '%.0f%s%2.2d%s%s' , ...
                                abs(grad(ii)) , char(176) , ...
                                abs(minu(ii)) , tcl{ii} , char(39)  );
         end
        
         tcl{ii} = sprintf('%s%s%s%s',pre,tcl{ii},lab{sign(tc(ii))+2},app);

      end

  end

%----------------------------------------------------

if ~isempty(ext) & isy
    if auto
       tic = deg2merc(tic,ext);
    else
       tic = tc;
    end
end

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  y = deg2merc(y,ext)

if ext ~= 0

   y = 180/pi * log(abs(tan(pi/4 + y*ext*pi/180/2))) / ext;

end


%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  y = merc2deg(y,ext)
 
if ext ~= 0

   y = 2*180/pi * ( atan(exp(y*ext*pi/180)) - pi/4 ) / ext; 

end


%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function warn(app,wrn)

% Display Warning-Message with Beep

ww = warnstat;

if strcmp(ww,'off') | isempty(wrn)
   return
end

% warning('on');
  
% fprintf(1,'\n%s',char(7));

fprintf(1,'\n%s: %s\n',app,wrn);

% fprintf(1,'\n%s',char(7));

% warning(ww);

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ww = warnstat

% WARNSTAT  Returns global WarningStatus
%
%  WARNSTAT returns the Status of WARNING
%
% Matlab R<13   WARNING
% Matlab R>12   WARNING for Identifier ALL
%

ww = warning;

if isstruct(ww)   % New Matlab R>12 Syntax
   try
      id = strcmp({ww.identifier},'all');
      if any(id)
         id = find(id);
         ww = ww(id(1)).state;
      else
         ww = '';
      end
   catch
      ww = '';
   end
elseif ~chkstr(ww)
   ww = '';
end

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [upp,axesize,x,y] = ppunit(axe);

% PPUNIT  Returns AxesUnitPerPixel and AxesPixelSize
%
%  [ UnitPerPixel, PixelSize, X , Y ] = PPUNIT( AxesHandle )
%
%  UnitPerPixel = [ XUnitsPerPixel  YUnitsPerPixel  ZUnitsPerPixel ] ;
%  PixelSize    = [ PixelLeft  PixelBottom  PixelWidth  PixelHight ];
%
%
%  Units --> Points:   
% 
%    [points] = [unit] / UnitsPerPixel / ScreenPixelPerInch * PointsPerInch
%
%   PointsPerInch = 72
%
%

if nargin == 0
   fig = get(0,'currentfigure');
   axe = get(fig,'currentaxes');
end

if isempty(axe)

   upp     = zeros(0,2);
   axesize = zeros(0,4);

   x = [];
   y = [];

   return

end


ok = ( isnumeric(axe) & ( prod(size(axe)) == 1 ) );
if ~ok
    ok = ishandle(axe);
    if ok
       ok = strcmp( get(axe,'type') , 'axes' );
    end
end

if ~ok
   error('Input must be an AxesHandle.');
end


  
axeuni  = get(axe,'units');
          set(axe,'units','pixels')
axesize = get(axe,'position');
          set(axe,'units',axeuni);


%********************************************
if ~isequal( get(axe,'view') , [ 0  90 ] );

  [upp,axesize,x,y] = ppunit3(axe,axesize);

  return

end
%********************************************


dx = diff(get(axe,'xlim'));
dy = diff(get(axe,'ylim'));

mode   = get(axe,'dataaspectratiomode');

if strcmp(mode,'manual')

   aspect = get(axe,'dataaspectratio');

   %  w/dx*ax == h/dy*ay 
   %
   % w = h * dx/dy * ay/ax;
   % h = w * dy/dx * ax/ay; 
   %

   pos = zeros(2,2);

   pos(1,2) = axesize(4);
   pos(1,1) = axesize(4) * dx/dy * aspect(2)/aspect(1);
   pos(2,1) = axesize(3);
   pos(2,2) = axesize(3) * dy/dx * aspect(1)/aspect(2);

   ii = find( sum( ( pos <= axesize([1 1],[3 4]) ) , 2 ) == 2 );

   pos = pos(ii(1),:);
 
   axesize([1 2]) = axesize([1 2])+axesize([3 4])/2-pos/2;
   axesize([3 4]) = pos;

end

 upp = [ dx  dy  0 ] ./ axesize([3 4 4]) ; % UnitPerPixel

%-----------------------------------
% Index for Cube

  ix = [ 1 2 2 1 1 2 2 1 ];
  iy = [ 1 1 2 2 1 1 2 2 ];

  x = [ 0  axesize(3) ];
  y = [ 0  axesize(4) ];
  
  x = x(ix);
  y = y(iy);

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [upp,axesize,x,y] = ppunit3(axe,axesize);

% PPUNIT  Returns AxesUnitPerPixel and AxesPixelSize for 3D-Axes
%


  xl = get(axe,'xlim');
  yl = get(axe,'ylim');
  zl = get(axe,'zlim');

   v = get(axe,'view');

  dx = diff(xl);
  dy = diff(yl);
  dz = diff(zl);


%-----------------------------------
% Index for Cube

  ix = [ 1 2 2 1 1 2 2 1 ];
  iy = [ 1 1 2 2 1 1 2 2 ];
  iz = [ 1 1 1 1 2 2 2 2 ];


  % Index for Axes in ix, iy, iz
  %  i01(1,:) --> i01(2,:)
  %       x  y  z

  i01 = [ 1  1  1 
          2  4  5 ];

%-----------------------------------

mode = get(axe,'dataaspectratiomode');

if strcmp(mode,'manual')

 x = xl(ix);
 y = yl(iy);
 z = zl(iz);

 aspect = get(axe,'dataaspectratio');

else

 x = ix-1;  % 0 | 1
 y = iy-1;
 z = iz-1;

 aspect = [ 1  1  1 ];

end

[x,y] = xy_proj( x , y , z , v , aspect );

% Normalize
x = x - min(x);
y = y - min(y);

xr = max(x);  % Range
yr = max(y);  % Range


if strcmp(mode,'manual')

   %  w/dx*ax == h/dy*ay 
   %
   % w = h * dx/dy * ay/ax;
   % h = w * dy/dx * ax/ay; 
   %

   pos = zeros(2,2);

   pos(1,2) = axesize(4);
   pos(1,1) = axesize(4) * xr/yr;
   pos(2,1) = axesize(3);
   pos(2,2) = axesize(3) * yr/xr;

   ii = find( sum( ( pos <= axesize([1 1],[3 4]) ) , 2 ) == 2 );

   pos = pos(ii(1),:);

   axesize([1 2]) = axesize([1 2])+axesize([3 4])/2-pos/2;
   axesize([3 4]) = pos;

end

 % Cube in Pixel

 x = x/xr * axesize(3);
 y = y/yr * axesize(4);

 % AxesLength in Pixel

 d01 = sqrt( ( x(i01(1,:)) - x(i01(2,:)) ) .^2 + ...
             ( y(i01(1,:)) - y(i01(2,:)) ) .^2   );
       

 upp = ~( d01 < 1 ) .* [ dx  dy  dz ] ./ ( d01 + ( d01 < 1 ) );

 
%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [x,y,z] = xy_proj(x,y,z,v,aspect)

% function  [x,y,z] = xy_proj(X,Y,Z,View,DataAspectRatio)
%
% View = [ Azimut Elevation ]
%

si = num2cell(size(x));

az = v(1)*pi/180; 
el = v(2)*pi/180;

aspect = aspect(:)' / min(aspect(:));

x = x(:)';
y = y(:)';
z = z(:)';

T = [ cos(az)           sin(az)           0
     -sin(el)*sin(az)   sin(el)*cos(az)  cos(el)
      cos(el)*sin(az)  -cos(el)*cos(az)  sin(el) ];

T = T ./ aspect(ones(3,1),:);

xyz = T * [ x ; y ; z ];

x = reshape(xyz(1,:),si{:});
y = reshape(xyz(2,:),si{:});
z = reshape(xyz(3,:),si{:});
 
