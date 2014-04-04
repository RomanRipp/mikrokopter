function h = gmtlines(cnf,varargin)

% GMTLINES  Plot of GMT-Lines from binned NetCDF-Files
%
% H = GMTLINES( CNF , [Axe] , [Area] , [Sample] )
%
% CNF = { GMT  LineColor  LineStyle  LineWidth }
%
%     GMT = { GMT_File  [Levels]  [+MinSkm] }
%     GMT = { GMT_File  [Levels]  [-MinSkm] }
%
%        NaN for MinSkm use the Value of OneSquarePixel
%        NaN +/- N*i for Min/MaxSkm use the Value of N Pixel
% 
%        See READ_GMT for GMT-Files and Levels
%
%        A first UPPERCASE Letter for GMT_File use SHORE_GMT
%          to extract closed Patches. Valid for ShoreLines only.
%
%     LineColor = [ R  G  B ]  |  ColorSpec
%
%     ColorSpec = Matlab's ColorSpecifier or XRGB-ColorName
%
% Area = [ LonMin LonMax  LatMin LatMax ]
%         default: [ XLim  YLim ] of Axes
%
% Sample = [ Sample  Smooth ]
%
%-------------------------------------------------------------
%
% Msg = GMTLINES( CNF ) Checks Configuration only
%
%-------------------------------------------------------------
%
% see also: READ_GMT, SHORE_GMT, SAMPLINE, COLSPEC, XRGB, 
%           PPUNIT, MERCATOR
%


h = [];

Nin  = nargin;
Nout = nargout;

if Nin == 0
   return
end

%----------------------------------------------------------------
% Check Inputs

msg = cell(0,1);

[m,cnf] = check(cnf);

if Nin == 1
   h = m;
   return
end

if ~isempty(m)
    msg = cat(1,msg,{m});
end

axe  = [];
area = [];
smp  = [];

cfg = {'none'};
cfg = cfg(ones(1,3));

%----------------------------------------------------------------

   if Nin > 4

      msg = cat( 1 , msg , {'To many InputArguments.'} ); 

   else

      for ii = 2 : Nin

          v = varargin{ii-1};

          p = prod(size(v));

          m = '';

          ok = isnumeric(v);

          if ok 
             %-----------------------------------------------------
             if ( p == 0 )
             %-----------------------------------------------------
                jj = strcmp(cfg,'none');
                if any(jj)
                   jj = min(find(jj));
                   cfg{jj} = 'inp';
                else
                   ok = NaN;
                end
             %-----------------------------------------------------
             elseif ( p == 1 )
             %-----------------------------------------------------
                ok = isequal(cfg{1},'none');
                if ok
                   cfg{1} = 'inp';
                   ok = ( isnumeric(v) & ( prod(size(v)) == 1 ) );
                   if ok
                      ok = ishandle(v);
                      if ok
                         ok = strcmp(get(v,'type'),'axes');
                      end
                   end
                   if ok
                      axe = v;
                   else
                      m = 'Invalid ##. Input for AxesHandle.';
                   end
                else
                   ok = NaN;
                end
             %-----------------------------------------------------
             elseif ( p == 2 )
             %-----------------------------------------------------
                ok = isequal(cfg{2},'none');
                if ok
                   cfg{2} = 'inp';
                   v(find(isnan(v))) = 0;
                   ok = ( all(isfinite(v)) & ( v >= 0 ) );
                   if ok
                      smp = v;
                   else
                      m = '##. Input for Sample must be positive finite.';
                   end
                else
                   ok = NaN;
                end
             %-----------------------------------------------------
             elseif ( p == 4 )
             %-----------------------------------------------------
                ok = isequal(cfg{3},'none');
                if ok
                   cfg{3} = 'inp';
                    v = v(:)';
                   ok = all(isfinite(v));
                   if ok
                      area = v;
                   else
                      m = '##. Input for Area must be finite.';
                   end
                else
                   ok = NaN;
                end
             %-----------------------------------------------------
             else
             %-----------------------------------------------------
                m = '##. Input has invalid Length.';
             end
             %-----------------------------------------------------
             if isnan(ok)
                m = 'Unrecognized ##. Input.';
             end
             %-----------------------------------------------------
          else
             m = '##. Input must be Numeric.';
          end

          if ~isempty(m)
              m = strrep(m,'##',sprintf('%.0f',ii));
              msg = cat(1,msg,{m});
          end

      end

   end

%----------------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    if Nout == 0
       error(msg);
    end
    h = msg;
    return
elseif isempty(cnf)
    return
end

%----------------------------------------------------------------
% Check for Axes

new = isempty(axe);

if new

   fig = get(0,'currentfigure');
   if isempty(fig)
      fig = figure;
   end

   axe = get(fig,'currentaxes');
   new = isempty(axe);

   if new
      axe = axes( 'parent' , fig , ...
                    'box'  , 'on'       );
   end

end

set( axe , 'nextplot' , 'add' );

ext = mercator(axe,'info');

if isempty(area) & ~new
   area = cat( 2 , get(axe,'xlim') , get(axe,'ylim') );
   if ~isempty(ext)
       area([3 4]) = mercator(axe,0,{area([3 4])});
   end
end

upp = ppunit(axe);

kpp = upp(1) * 60 * 1.852;  % km per pixel

skm = kpp^2;            % MinSquareKilometer

if ~isempty(smp)
    smp = smp(1) + smp(2)*i;
end

%----------------------------------------------------------------
% Read and Plot GMT-Lines

n =  size(cnf,1);
h =  cell(n,1);
m = zeros(n,1);

h(:) = {zeros(0,1)};

fcn = { 'read_gmt'  'shore_gmt' };
 
for ii = 1 : n
 
    if ~isempty(cnf{ii,1}{3})
        if ~isfinite(real(cnf{ii,1}{3}))
            npx = imag(cnf{ii,1}{3});
            npx = npx + ( npx == 0 );
            cnf{ii,1}{3} = npx * skm;
        end
        cnf{ii,1}{3} = cnf{ii,1}{3} * i;  % Imag for using READ_GMT
    end

    shore = ( cnf{ii,1}{1}(1) == upper(cnf{ii,1}{1}(1)) );

    try
       xy = feval(fcn{1+shore},cnf{ii,1}{1},area,cnf{ii,1}{2},cnf{ii,1}{3});
    catch
       xy = [];
       warn(sprintf('Error using %s.\n%s',upper(fcn{1+shore}),lasterr));
    end

    if ~( isempty(xy) | shore )
        xy = {permute(xy,[2 1])};  % [ X ; Y ]
    end

    m(ii) = prod(size(xy));

    if ~isempty(xy) & ~isempty(smp)
        try
           for jj = 1 : m(ii)
               xy{jj} =  sampline(xy{jj},smp,'cos',1e-6,2);
           end
        catch
           xy = [];
           warn(sprintf('Error using SAMPLINE.\n%s',lasterr));
        end
    end

    if ~isempty(xy)

        fc = { 'none' cnf{ii,2} };
        fc = fc{ 1 + shore };

        h{ii} = zeros(m(ii),1);

        for jj = 1 : m(ii)

            if ~isempty(ext)
                xy{jj}(2,:) = mercator(axe,{xy{jj}(2,:)});
            end

            h{ii}(jj) = patch( 'parent' , axe   , ...
                          'xdata' , xy{jj}(1,:) , ...
                          'ydata' , xy{jj}(2,:) , ...
                      'facecolor' , fc          , ...
                      'edgecolor' , cnf{ii,2}   , ...
                      'linestyle' , cnf{ii,3}   , ...
                      'linewidth' , cnf{ii,4}   , ...
                         'marker' , 'none'      , ...
                            'tag' , 'GMT'                );

            setappdata( h{ii}(jj) , 'GMT' , cnf{ii,1} );

        end

    end

end

jj = ( m > 0 );

if      all(jj)

    h = cat( 1 , h{:} );

elseif  any(jj)

   jj = find(jj);
    h = cat( 1 , h{jj} );

else

    h = [];

end


%***********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,cnf] = check(cnf)

%  { GMT  LineColor LineStyle LineWidth   }
%  LineColor ==  [ R G B ] | ColorSpec | [ -1 .. 1 ]  Scaled Dark/Bright


msg = '';
if isempty(cnf)
   return
end

s = size(cnf);
p = prod(s);

if chkstr(cnf,1)
   cnf = {cnf};
end

if ~( iscell(cnf) & ( ndims(cnf) == 2 ) );
    msg = 'Configuration must be a 2D-CellArray.';
    return
end

nc = size(cnf,1);
mc = size(cnf,2);

dc = { 'sc1'  'k' '-' 0.5 };  % Defaults

nd = size(dc,2);

if    mc > nd

   cnf = cnf(:,1:4);

elseif mc < nd

   cnf = cat( 2 , cnf , dc( ones(nc,1) , (mc+1) : nd ) );

end

msg = cell(0,1);

ls = set(0,'defaultpatchlinestyle');

ok = zeros(nc,4);

fn = { 'GMT'  'Color'  'Style'  'Width' };

for ii = 1 : size(cnf,1)

    %--------------------------------------------------------------------
    % GMT
    
    ok(ii,1) = chkstr(cnf{ii,1});
    if ok(ii,1)
       cnf{ii,1} = { cnf{ii,1} []  [] }; 
    else
       psi = prod(size(cnf{ii,1}));
       ok(ii,1) = ( iscell(cnf{ii,1}) & ~isempty(cnf{ii,1}) & ( psi <= 3 ) );
       if ok(ii,1);
          ok(ii,1) = chkstr(cnf{ii,1}{1},1);
          if ok(ii,1) 
             if psi == 1
                cnf{ii,1} = [ cnf{ii,1}  { []  [] } ];
             elseif psi == 2
                cnf{ii,1} = [ cnf{ii,1}(:)'  {[]} ];
             end
             ok(ii,1) = ( isnumeric(cnf{ii,1}{2}) & isnumeric(cnf{ii,1}{3}) & ...
                          ( isempty(cnf{ii,1}{3}) | ( prod(size(cnf{ii,1}{3})) == 1 ) ) );
             if ok(ii,1) & ~isempty(cnf{ii,1}{2})
                ok(ii,1) = all(isfinite(cnf{ii,1}{2}));
             end
          end
       end
    end
       
    %--------------------------------------------------------------------
    % LineColor

    ok(ii,2) = isnumeric(cnf{ii,2});

    if ok(ii,2)
       s = size(cnf{ii,2});
       ok(ii,2) = isequal(s,[1 3]);
       if ok(ii,2)
          ok(ii,2) = all( abs(cnf{ii,2}-0.5) <= 0.5 );
       end
    else
       ok(ii,2) = chkstr(cnf{ii,2},1);
       if ok(ii,2)
         cnf{ii,2} = colspec(cnf{ii,2},1);
          ok(ii,2) = ~isempty(cnf{ii,2});
       end
    end

    %--------------------------------------------------------------------
    % LineStyle

    ok(ii,3) = chkstr(cnf{ii,3},1); 

    if ok(ii,3)
       ok(ii,3) = any( strcmp( cnf{ii,3} , ls ) );
    end

    %--------------------------------------------------------------------
    % LineWidth

    ok(ii,4) = ( isnumeric(cnf{ii,4}) & ( prod(size(cnf{ii,4})) == 1 ) );
    if ok(ii,4)
       ok(ii,4) = ( cnf{ii,4} > 0 );
    end

    %--------------------------------------------------------------------

    if ~all(ok(ii,:))
        jj = find(~ok(ii,:));
        m = sprintf(' %s ',fn{jj});
        m = sprintf('Invalid %.0f. Configuration: %s',ii,m);
        msg = cat( 1 , msg , {m} );
    end

end

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
else
    msg = '';
end

%***********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [n,c] = colspec(v,x)

% COLSPEC   Returns Matlab's or XRGB color specifier 
%
% Get List of Colors:
%
% [ Name , RGB ] = COLSPEC      List of ColorSpecifier
% 
% [ Name , RGB ] = COLSPEC(1)   List of ColorSpecifier and XRGB-Colors
% 
% Find a Color:
%
%   RGB  = COLSPEC( Name )       Color corresponds with Name, Matlab's only
%   RGB  = COLSPEC( Name , 1 )   Color corresponds with Name, incl. XRGB-ColorNames
%
%   Name = COLSPEC( RGB )        Find Name corresponds with Color, Matlab's only
%
%   Name = COLSPEC( RGB , 1+R*i) Find Name corresponds with Color, incl. XRGB
%
%   R is the ToleranceRadius to search for a XRGB-Color, default: 0.03
%
% see also: XRGB
%

Nin = nargin;

nv = 0;

if Nin == 0
   nv = 1;        % Not V
   v  = [];
   x  = 0;
elseif Nin < 2
   nv = ( isnumeric(v) & ( prod(size(v)) == 1 ) );
   if nv
      x = v;
   else
      x = 0;
   end
end

rad = imag(x);       % ToleranceRadius

x   = isequal(real(x),1);  % Use XRGB

%------------------------------------------------------

n = {   'r'        [ 1  0  0 ]
        'g'        [ 0  1  0 ]
        'b'        [ 0  0  1 ]
        'c'        [ 0  1  1 ]
        'm'        [ 1  0  1 ]
        'y'        [ 1  1  0 ]
        'k'        [ 0  0  0 ]
        'w'        [ 1  1  1 ]
        'red'      [ 1  0  0 ]
        'green'    [ 0  1  0 ]
        'blue'     [ 0  0  1 ]
        'cyan'     [ 0  1  1 ]
        'magenta'  [ 1  0  1 ]
        'yellow'   [ 1  1  0 ]
        'black'    [ 0  0  0 ]
        'white'    [ 1  1  1 ]   };

c = cat(1,n{:,2});
n = n(:,1);

%------------------------------------------------------

if nv
   m = size(n,1) / 2;         % Unique Colors
   c = c(1:m,:);
   n = n(1:m);
   if x
      try
         [nx,cx] = xrgb(1);      % Unique Colors
         n  = cat( 1 , n , nx );
         c  = cat( 1 , c , cx );
      catch
         warn(sprintf('Can''t get ColorList by XRGB.\n%s',lasterr))
      end
   end
   return
end

%------------------------------------------------------

s = size(v);
p = prod(s);

%----------------------------------------------
if isnumeric(v) & isequal(s,[1 3])
%----------------------------------------------

   d = c - v(ones(1,size(c,1)),:);
   d = all(d==0,2);

   if  any(d)

       d = find(d);
       d = d(1);
       n = n{d};
       c = c(d,:);

   elseif x

       try

          rad = rad + 0.03 * ( rad == 0 );

          [n,c] = xrgb(1);  % Unique Colors

           m = ones(size(c,1),1);

           d = ( c - v(m,:) );

           d = sum( d.^2 , 2 );

           ok = ( d <= rad^2 );

           if any(ok)
              nk =  sum(ok);
              ok = find(ok);
              if nk > 1
                 [d,ii] = min(d(ok));
                  ok = ok(ii);
              end
              n = n{ok};
              c = c(ok,:);
           else
              n = '';
              c = [];
           end

        catch

           warn(sprintf('Can''t get Color by XRGB.\n%s',lasterr))
           
           n = '';
           c = [];

        end

   else
       n = '';
       c = [];
   end

%----------------------------------------------
elseif ischar(v) & ~isempty(v) & ( p == s(2) )
%----------------------------------------------

   v = lower(v);

   for ii = 1 : size(n,1)

       nn = n{ii}( 1 : min(p,size(n{ii},2)) );

       ok = strcmp( nn , v );      

       if ok
          break
       end

   end


   if  ok
       n = n{ii};
       c = c(ii,:);
   elseif x
       try
          [c,h,n] = xrgb(v);
       catch
          warn(sprintf('Can''t get ColorList by XRGB.\n%s',lasterr))
          n = '';
          c = [];
       end
   else
       n = '';
       c = [];
   end

   v = c;
   c = n;
   n = v;

%----------------------------------------------
else
%----------------------------------------------

   error('Input must be an RGB-Tripel or String.');

%----------------------------------------------
end
%----------------------------------------------

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function warn(wrn)

% Display Warning-Message with Beep

ww = warnstat;

if strcmp(ww,'off') | isempty(wrn)
   return
end

warning('on');
  
fprintf(1,'\n%s',char(7));

warning(wrn);

fprintf(1,'\n%s',char(7));

warning(ww);

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

ok = ( sum(isnan(x),2) == sr );   % Ok for NaN-Row

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

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
 
