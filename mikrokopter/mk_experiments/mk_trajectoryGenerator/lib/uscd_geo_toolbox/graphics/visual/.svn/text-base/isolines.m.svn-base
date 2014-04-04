function h = isolines(cnf,x,y,z,axe,fcn)

% ISOLINES  Extended Contour-Plot
%
% H = ISOLINES( CNF , X , Y , Z , [Axe] , [Fcn] )
%
% H = Vector of created PatchHandles or ErrorMessage !!!
% 
% CNF = Vector for IsoValues or a 4-Column CellArray with
%       IsoValues and LineStyles:
%
%    { IsoValue(s)  LineColor  LineStyle  LineWidth }
%
%    IsoValue+i ==> 'modulo' for same SIGN:
%                  ( MOD(Z,V) == 0 ) & ( SIGN(Z) == SIGN(V) )
%
%    IsoValue-i ==> 'modulo': ( MOD(Z,V) == 0 )
%
%    LineColor = [ R  G  B ]  |  ColorSpec  |  Scale
%
%    ColorSpec = Matlab's ColorSpecifier or XRGB-ColorName
%
%    Scale    = [ -1 .. 1 ], 
%      get Color from FigureColorMap and AxesColorLimit  
%      or get Color by the Function from last optional Input Fcn.
%
%    Scale < 0  ==>    dark Color (lower Value)
%    Scale > 0  ==>  bright Color (lower Saturation)
%
% X = Vector or Matrice for X-Coordinates
% Y = Vector or Matrice for Y-Coordinates
% Z =           Matrice for Z-Coordinates (IsoValues)
%
% A Vector for X should correspond with the 2. Dimension of Z
% A Vector for Y should correspond with the 1. Dimension of Z
%
% Axe = AxesHandle to plot the IsoLines
%
% Fcn = CellArray with FunctionName to get the Color
%       in case LineColor == Scale by using FEVAL:
%
%       C = FEVAL( Fcn{:} , IsoValue )
%
%       C = [ R  G  B ]  |  ColorIndex 
%
% For 3D-Isolines use a NonZero imaginary Part
%  of the 5. Input:
%
% ISOLINES( CNF , X , Y , Z , +i , [Fcn] ) 
% ISOLINES( CNF , X , Y , Z , Axe+i , [Fcn] ) 
%
%-------------------------------------------------------------
%
% Msg = ISOLINES( CNF ) Checks Configuration only
%
%-------------------------------------------------------------
%
% see also: CONTOURS, MKGRID, COLSPEC, XRGB, RGB2HSV, HSV2RGB
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

%----------------------------------------------------------------

   if Nin < 4

      msg = cat( 1 , msg , {'X, Y and Z required.'} ); 

   else

      ok = ( isnumeric(x) & isnumeric(y) & isnumeric(z) );
      ok = ( ok & all( [ndims(x) ndims(y) ndims(z)] == 2 ) );

      if ~ok
          msg = cat( 1 , msg , {'X, Y and Z must be 2D-Numerics.'} );
      else
         [ok,y,x,z] = mkgrid(y,x,z);
         if ~ok
             msg = cat( 1 , msg , {'Matrix Dimensions must be agree.'} );
         elseif ~( ( ndims(z) == 2 ) | isempty(z) )
             msg = cat( 1 , msg , {'2D-Inputs required.'} );
         end
      end

   end

   
   if Nin < 6
      fcn = {};
   end

   if Nin < 5
      axe = [];
   elseif iscell(axe)
      if ~iscell(fcn)
         ax = fcn;
         fcn = axe;
         axe = ax;
      elseif Nin < 6
         fcn = axe; 
         axe = [];
      end
   end

   v3d = 0;

   if ~isempty(axe)
       ok = ( isnumeric(axe) & ( prod(size(axe)) == 1 ) );
       if ok
          v3d = ~( imag(axe) == 0 );
          axe = real(axe);
          if axe == 0
             axe = [];
          else
             ok = ishandle(axe);
             if ok
                ok = strcmp(get(axe,'type'),'axes');
             end
          end
       end
       if ~ok
           msg = cat( 1 , msg , {'Invalid Input for AxesHandle.'} );
       end
   end

   if ~isempty(fcn)
       ok = iscell(fcn);
       if ok
          ok = chkstr(fcn{1},1);
          if ok
             f = which(fcn{1});
             ok = ~isempty(f);
             if ok
                ok = ( exist(f,'file') == 2 );
                if ok
                   [p,f] = fileparts(f);
                   ok = strcmp(f,fcn{1});
                end
             end
          end
       end
       if ~ok
           msg = cat( 1 , msg , {'Invalid Input for ColorFunction.'} );
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
end

%----------------------------------------------------------------
% Check Grid for finite Elements and NOT singular

ok = ~isempty(z);

if ok

   xl = [ min(x(:)) max(x(:)) ];
   yl = [ min(y(:)) max(y(:)) ];
   zl = [ min(z(:)) max(z(:)) ];
 
   ok = ( all(isfinite(xl)) & ( xl(1) < xl(2) ) & ...
          all(isfinite(yl)) & ( yl(1) < yl(2) ) & ...
          all(isfinite(zl)) & ( zl(1) < zl(2) ) );

end

if isempty(cnf) | ~ok
   return
end

%----------------------------------------------------------------
% Check for Axes

if isempty(axe)

   fig = get(0,'currentfigure');
   if isempty(fig)
      fig = figure;
   end

   axe = get(fig,'currentaxes');
   new = isempty(axe);
   if new
      axe = axes( 'parent' , fig  , ...
                    'xlim' ,  xl  , ...
                    'ylim' ,  yl  , ...
                    'clim' ,  zl  , ...
                    'box'  , 'on'       );
      if v3d
         set(axe,'zlim',zl);
      end
   end

end

if strcmp( get(axe,'climmode') , 'auto' )
   set( axe , 'clim' , zl );
end

set( axe , 'nextplot' , 'add' );

%----------------------------------------------------------------
% Sort IsoValues

%  { Z-Value  LineColor LineStyle LineWidth   }
%  Z-Value+i ==> 'mod'
%  Z-Value+i ==> 'mod' at sign
%  LineColor ==  [ R G B ] | ColorSpec | [ -1 .. 1 ]  Dark or Bright

nn = size(cnf,1);

cz = cat(1,cnf{:,1});

ci = sign(imag(cz));         % [ -1  1   0 ]

cz = real(cz);
sg = sign(cz);

ci = ci .* abs(sg);          % Take care for MOD at ZERO

cf = ( 1 + ( sg - 1 ) .* ( ci == -1 ) );  % Make Absolute if "-i"

cz = cz .* cf;
sg = sg .* cf;

ci = ci - 1.5 * ( ci == 1 ); % [ -1 -.5  0 ]
ci = 2 * ci + 1;             % [ -1  0   1 ]

[cz,si] = sort(abs(cz));        % Sort by absolute Value

ci = ci(si);
sg = sg(si);

cnf = cnf(si,:);

[ci,si] = sort(ci);             % ABS last

cz = cz(si);
sg = sg(si);

cnf = cnf(si,:);

%-------------------------------------------------------
% Get IsoLines

cl =  get(     axe           , 'clim'     );
cm =  get( get(axe,'parent') , 'colormap' );

nc = size(cm,1);

dc = diff(cl) / nc;

cl = linspace( cl(1) , cl(2) , nc ) + dc/2;


hz    = cell(nn,1);

hz(:) = { ones(2,0) };  % [ Handle ; Value ]


for ii = 1 : nn

    %------------------------------------
    % Check for Previous

    ip = [];
    ok = ( ii > 1 );
    if ok
       jj = ( 1 : ii-1 );
       kk = ( ( sg(jj) == sg(ii) ) | ( ci(jj) == -1 ) );
       ok =  any(kk);
       if ok
          kk = jj(find(kk));
          for jj = kk(:)'
              ok = ( mod(cz(ii),cz(jj)) == 0 );
              ok = ( ( ok & ~( ci(jj) == 1 ) ) | ( cz(ii) == cz(jj) ) );
              if ok
                 break
              end
          end
          if ok
             ip = jj;
          end
       end
    end

    %------------------------------------
    % Get / Create Lines

    hh = [];

    if ok & ~isempty(ip)

       if ( ci(ii) == 1 )
          ok = ( hz{ip}(2,:) == cz(ii) );
       else
          ok = ( mod(hz{ip}(2,:),cz(ii)) == 0 );
          if ci(ii) == 0
             ok = ( ok & ( sign(hz{ip}(2,:)) == sg(ii) ) );
          end
       end

       if any(ok)
          ok = find(ok);
          hh = hz{ip}(1,ok);
       end

    else

       if cz(ii) == 0

          if ( zl(1) > 0 ) | ( zl(2) < 0 )
             vc = [];
          else
             vc = [ 1  1 ];
          end

       elseif ci(ii) == -1

          vc = ( ceil(zl(1)/cz(ii)) : floor(zl(2)/cz(ii)) );

       else

          iz = ( sg(ii) + 3 ) / 2;
          zc =   sg(ii) * zl(iz);

          if zc < cz(ii)
             vc = [];
          elseif ( ci(ii) == 1 )  % ABS
             vc = [ 1  1 ];
          else
             vc = ( 1 : floor(zc/cz(ii)) );
          end

       end

       if ~isempty(vc)

           vc = sg(ii) * cz(ii) * vc;

           if prod(size(vc)) == 1
              vc = vc([1 1]);
           end

           cs = contours(x,y,z,vc);

           if ~isempty(cs)

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

               i1 = cs(2,i0);  % Length
               lv = cs(1,i0);  % Value

               n2 = NaN * zeros(2,1);

               hh = zeros(1,nl); % Handles

               for jj = 1 : nl

                   kk = ( 1 : i1(jj) ) + i0(jj);

                   xy = cat( 2 , cs(:,kk) , n2 );

                   hh(jj) = patch( 'parent' , axe      , ...
                                    'xdata' , xy(1,:)' , ...
                                    'ydata' , xy(2,:)' , ...
                                'facecolor' , 'none'   , ...
                                'edgecolor' , 'k'      , ...
                                'clipping'  , 'on'     , ...
                                'userdata'  , lv(jj)         );
                   if v3d
                      set(hh(jj),'zdata',lv(jj)+0*xy(1,:)');
                   end
               end

               hz{ii} = cat( 1 , hh , sg(ii)*lv );

           end

       end

    end

    %------------------------------------
    % Set LineStyle

    if ~isempty(hh)

        set( hh , 'linestyle' , cnf{ii,3} , ...
                  'linewidth' , cnf{ii,4}       );

        cc = cnf{ii,2};

        if isnumeric(cc) & ( prod(size(cc)) == 1 )

           cj = 2 + ( cc < 0 );  % [ Sat | Value ]
           cv = 1 - abs(cc);

           for h = hh

               v = get(h,'userdata');  % Value

               if ~isempty(fcn)

                   [msg,cc] = fcn_color(fcn,v);
                   if ~isempty(msg)
                       if Nout == 0
                          error(msg);
                       end
                       h = msg;
                       return
                   end
                   if prod(size(cc)) == 1
                      cc = min(max(cc,1),nc);
                      cc = cm(cc,:);
                   end

               else

                   ic = sum( cl < v );
                   if any( ic == [ 0  nc ] )
                      cc = cm(ic+(ic==0),:);
                   else
                      ik = ic + 1;
                      cc = ( v - cl(ic) ) * ( cm(ik,:) - cm(ic,:) ) / dc;
                      cc = cc + cm(ic,:);
                   end

               end

               if  cv < 1 
                   cc = rgb2hsv(cc);
                   cc(cj) = cc(cj) * cv;
                   cc(3)  = cc(3) + ( 1 - cc(3) ) * ( 1 - cv ) * ( cj == 2 );
                   cc = hsv2rgb(cc);
               end

               set( h , 'edgecolor' , cc );

           end

        else

           set( hh , 'edgecolor' , cc );

        end

    end
              
end

h = cat( 2 , hz{:} );

if isempty(h)
   h = [];
else
   h = h(1,:);
   n = size(h,2);
   if n > 1

      % Sort by Level
      lv = get(h,'userdata');
      lv = cat( 1 , lv{:} );
     [lv,si] = sort(-lv);

      % Set AxesChildren
      ch = get( axe , 'children' );
      ic = zeros(n,1);
      for ii = 1 : n
          ic(ii) = find( ch == h(ii) );
      end

      h = h(si);

      ch(ic) = h;
      set(axe,'children',ch);

   end
end

%***********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,cnf] = check(cnf)

%  { Z-Value  LineColor LineStyle LineWidth   }
%  Z-Value+i ==> 'mod'
%  LineColor ==  [ R G B ] | ColorSpec | [ -1 .. 1 ]  Scaled Dark/Bright


msg = '';
if isempty(cnf)
   return
end

s = size(cnf);
p = prod(s);

if isnumeric(cnf)
   if ( p == max(s) )
      cnf = num2cell(cnf(:));
        s = [ p  1 ];
   else
      msg = 'Numeric Configuration must be a Vector.';
      return
   end
end

if ~( iscell(cnf) & ( ndims(cnf) == 2 ) );
    msg = 'Configuration must be a 2D-CellArray.';
    return
end

nc = size(cnf,1);
mc = size(cnf,2);

dc = { 0 'k' '-' 0.5 };  % Defaults

nd = size(dc,2);

if    mc > nd

   cnf = cnf(:,1:4);

elseif mc < nd

   cnf = cat( 2 , cnf , dc( ones(nc,1) , (mc+1) : nd ) );

end

msg = cell(0,1);

ls = set(0,'defaultpatchlinestyle');

ok = zeros(nc,4);
nn = zeros(nc,1);   % Number of IsoValues

fn = { 'Value'  'Color'  'Style'  'Width' };

for ii = 1 : nc

    %--------------------------------------------------------------------
    % ZValue

    ok(ii,1) = isnumeric(cnf{ii,1});
    if ok(ii,1)
       cnf{ii,1} = cnf{ii,1}(:);
        ok(ii,1) = all(isfinite(cnf{ii,1}));
        nn(ii)   = size(cnf{ii,1},1);
    end
       
    %--------------------------------------------------------------------
    % LineColor

    ok(ii,2) = isnumeric(cnf{ii,2});

    if ok(ii,2)
       s = size(cnf{ii,2});
       ok(ii,2) = isequal(s,[1 3]);
       if ok(ii,2)
          ok(ii,2) = all( abs(cnf{ii,2}-0.5) <= 0.5 );
       else
          ok(ii,2) = isequal(s,[1 1]);
          if ok(ii,2)
            cnf{ii,2} = sign(cnf{ii,2}) * min(abs(cnf{ii,2}),1);
          end
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

%---------------------------------------------
% Check for empty or multiple IsoValues

jj = ( nn == 0 );
if any(jj)
   jj = find(jj);
   cnf(jj,:) = [];
end

jj = ( nn > 1 );
if ~any(jj)
    return
end

jj = -sort(-find(jj));

for ii = jj(:)'

    cc = cell(nn(ii),4);
    cc(:,1) = num2cell(cnf{ii,1});
    cc(:,[2 3 4]) = cnf(ii*ones(1,nn(ii)),[2 3 4]);

    cnf = cat( 1 , cnf(1:(ii-1),:) , cc , cnf((ii+1):end,:) );

end

%***********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,cc] = fcn_color(fcn,v)


msg = '';
cc  = zeros(1,3);

try
    cc = feval(fcn{:},v);
catch
    msg = sprintf('Error get Color by %s.\n%s',upper(fcn{1}),lasterr);
end

if ~isempty(msg)
    return
end

pc = prod(size(cc));
ok = ( isnumeric(cc) & any( pc == [ 1  3 ] ) );

if ok
   if pc == 1
      ok = ( mod(cc,1) == 0 );
   else
      if strcmp(class(cc),'uint8')
         cc = double(cc) / 255;
      end
      cc = cc(:)';
      ok = all( ( 0 <= cc ) & ( cc <= 1 ) );
   end
end

if ~ok
    msg = sprintf('Invalid Color by %s.',upper(fcn{1}));
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

%***********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,varargout] = mkgrid(varargin);

% MKGRID  Generation of common Grid from Inputs
%
% [ Ok, Y1 , Y2 , Y3 , ... ] = MKGRID( X1 , X2 , X3 , ... )
%
%   Ok =  0   if MatrixDimensions are not agree
%             [ Y1 , Y2 , ... ] == [ X1 , X2 , ... ]
%
%         1   if MatrixDimensions are agree
%             [ Y1 , Y2 , ... ] are Matrices with equal Size
%
%        NaN  if missing Input
%
%
%   MatriceInputs must have an equal Size.
% 
%   VectorInputs are checked successive if they match the Length
%    of the previous Vectors and all MatriceInputs, 
%     or the Length in the DimensionOrder of InputOrder. 
%
%
% MKGRID( '#<MemoryLimit>' , ... )  sets a Limit 
%
%   of the Number of OutputElements to: MemoryLimit * 2^20 / 8
%
%   The Unit of MemoryLimit is [Mbytes] in case of class "double".
%
%   In case of exeeding MemoryLimit, VectorInputs returns reshaped
%    in correct DimensionOrder but not gridded !
%
%   default: MemoryLimit = 256 [Mbytes] ==  2 x [ 4096 x 4096 ] double
%                                       ==  8 x [ 2048 x 2048 ] double
%                                       == 32 x [ 1024 x 1024 ] double
%
%   example: MKGRID( '#512' , ... )   % sets the MemoryLimit to 512 Mbyte
%
%            The Maximum Number of Elements in Output is:  512 * 2^20 / 8
%
% 
% See also: NDGRID, MESHGRID
%

%-------------------------------------------------

mb = 256;     % default MemoryLimit [Mbytes] !!!!!!

bp = char(7); % WarningBeep

%-------------------------------------------------

 in =  nargin;
out = nargout - 1;

out = max(0,out);

varargout = cell(1,out);

%-------------------------------------------------
% Check for 1. MemoryInput

if in > 0

   v = varargin{1};
   s = size(v);
   p = prod(s);

   ok = ( ischar(v) & ( s(2) == p ) & ( s(2) > 1 ) );
   if ok
      ok = ( v(1) == '#' );
   end

   if ok

      try 
         m = eval(v(2:s(2)));
      catch
        ok = 0;
         m = 0;
      end

      wrn = '';

      if ~ok
          wrn = sprintf('Cann''t evaluate Expression for MemoryLimit "%s".',v);
      else
         ok = ( isnumeric(m) & ( prod(size(m)) == 1 ) );
         if ok
            ok = ( m > 0 );
         end
         if ~ok
             wrn = sprintf('MemoryLimit in "%s" must be a single numeric larger ZERO.',v);
         end
      end

      if ~ok
          ww = warnstat;
          if ~strcmp(ww,'off')
              warning('on')
              warning([ bp  wrn ]);
              warning(ww);
          end
      else
          varargin = varargin(2:in);
          in = in - 1;
          mb = m;
      end

   end

end

%-------------------------------------------------
         
if in < 1
   ok = NaN;
   return
end

out = min(in,out);

%-------------------------------------------------
% Check NDim

d  = zeros(in,1);

for ii = 1 : in
    d(ii) = ndims(varargin{ii});
end

nd = max(d);  % NDIM

%-------------------------------------------------
% Get Sizes

sz = ones(in,nd);

for ii = 1 : in
    sz(ii,1:d(ii)) = size(varargin{ii});
end

%-------------------------------------------------
% Check all Equal 

ok = sz - sz(ones(1,in),:);
ok = ( ok == 0 );
ok = all(ok(:));

if ok
   for ii = 1 : out
       varargout{ii} = varargin{ii};
   end
   return
end

%-------------------------------------------------
% Check for Vectors

mv = max(sz,[],2);         % Maximum per Matrice

v = ( mv == prod(sz,2) );  % True for Vector

isv = any(v);

if isv

   vd = double( mv(:,ones(1,nd)) == sz );

   vd = sum(cumprod(1-vd,2),2) + 1;

   vd = vd .* v;            % Orig VectorDimension

   wd = vd;                 % New  VectorDimension

   %-----------------------------------
   % Sort VectorDimension by InputOrder

   sv         = ones(in,max(in,nd));
   sv(:,1:nd) = sz;

      iv    = find(v);
   sv(iv,:) = 1;

   sv(iv+(iv-1)*in) = mv(iv);

else

   iv = [];

   vd = zeros(in,1);
   wd = zeros(in,1);
   sv = sz;

end

%-------------------------------------------------
% Check Length of Dimensions

ok = chkdimsize(sz);

%-------------------------------------------------
% Check Vectors successive

if ~ok & isv

    %-------------------------------------------------
    % Loop over Vectors in InputOrder
    % Check with previous Vectors and all Matrices

    for ii = iv(:)'

        w       = zeros(in,1);
        w(1:ii) = v(1:ii);      % Vector until II
        
        iw = ( ~v | w );        % Vector until II or Matrice
        nw =  sum(iw);
        iw = find(iw);

        %------------------------------------------------------
        % Check Size Including Vector

        ok = chkdimsize(sz(iw,:));

        if ~ok & ( nw > 1 )

            %------------------------------------------------------
            % Check Length of Vector with Size exluding Vector

            iw = iw( 1 : (nw-1) );

            d  = max(sz(iw,:),[],1);  % Maximum per Dimension

            jj = ( d == mv(ii) );     % Length of Vector match Maximum

            ok = any(jj);

            if ok

                     wd(ii)  = sum(cumprod(double(~jj),2),2) + 1;
               sz(ii,  :   ) = ones(1,nd);
               sz(ii,wd(ii)) = mv(ii);

            elseif nd < ii                     % Expand NDIM

               sz = cat(2,sz,ones(in,ii-nd));
               nd = ii;

                     wd(ii)  = ii;
               sz(ii,  :   ) = ones(1,nd);
               sz(ii,wd(ii)) = mv(ii);

            end

        end

    end

    %-------------------------------------------------
    % Check Length of Dimensions again

    ok = chkdimsize(sz);

    %-------------------------------------------------
    % Check Length of Dimensions, 
    %  using Input-Order sorted VectorDimensions

    if ~ok

        ok = chkdimsize(sv);

        if ok
           wd(iv) = iv;
               sz = sv;
               nd = size(sv,2);
        end

    end
    %-------------------------------------------------
 
end

%-------------------------------------------------
% Check for Changed VectorDimensions

if ok & ( out < in ) & isv

   jj = ( (out+1) : in );

   if ~all( vd(jj) == wd(jj) ) 
       ww = warnstat;
       if ~strcmp(ww,'off')
           warning('on')
           warning([ bp  'VectorDimension not agree.' ]);
           warning(ww);
       end
   end

end

%-------------------------------------------------
% Check for Return

if ~ok | ( out == 0 )
    for ii = 1 : out
        varargout{ii} = varargin{ii};
    end
    return
end

%-------------------------------------------------
% Check for VectorPermutation

jj = ( 1 : out );

jj = ~( vd(jj) == wd(jj) );

if any(jj)

   jj = find(jj);

   for ii = jj(:)'

       varargin{ii} = varargin{ii}(:);

       prm = ( 1 : nd ) - wd(ii) + 1;  % 1 at wd

       prm = prm - nd * ( ceil(prm/nd) - 1 );
 
       varargin{ii} = permute(varargin{ii},prm);

   end

end

%-------------------------------------------------
% Check Memory

d = max(sz,[],1);

m = out * prod(d,2) * 8 / 2^20;

if m > mb

   ww = warnstat;
   if ~strcmp(ww,'off')
       warning('on')
       warning(sprintf('%sMemoryLimit of %.3g Mbytes exeeded.',bp,mb));
       warning(ww);
   end

    for ii = 1 : out
        varargout{ii} = varargin{ii};
    end

    return

end

%-------------------------------------------------
% Build Matrices

o = cell(1,nd);
k = cell(1,nd);

for ii = 1 : nd
    o{ii} = ones(1,d(ii));
    k{ii} = ( 1 : d(ii) );
end

for ii = 1 : out

     q = k;

    jj = ( sz(ii,:) == 1 );

    if any(jj)
           jj  = find(jj);
         q(jj) = o(jj);
    end

    varargout{ii} = varargin{ii}(q{:});

end


%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = chkdimsize(s);

% Check if all DimLength's are to equal Maximum or ONE
%

d = max(s,[],1);                   % Maximum per Dimension

ok = s - d(ones(1,size(s,1)),:);

ok = ( ( ok == 0 ) | ( s == 1 ) );

ok = all(ok(:));

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

