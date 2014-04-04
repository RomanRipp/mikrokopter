function [cnf,cfg] = showseq(d,varargin);

% SHOWSEQ  Subplot of Frames
%
% [CNF,CFG] = SHOWSEQ( Data )
%
%  Data = [ Ny by Nx by NFrame ] - Matrice
%
% [CNF,CFG] = SHOWSEQ( Data , Property1 , Value1 , ... )
%
% Properties:
% 
% 'type'   'image' | 'surface'  Type of Grafic,                           default: 'image'
%                                Image works faster
%                                Surface for low-resolution data
%
% 'size'   [ szx  szy  ]          MaxWidth and MaxHeight of Figure [pixel], default: [ 640  480 ]
%
% 'xvec'   [ 1 by Nx ]          Vector for X,                             default: ( 1 : Nx )
% 'yvec'   [ 1 by Ny ]          Vector for Y,                             default: ( 1 : Ny )
%
% 'scale'  [ scx  scy  ]          ScaleFactor for X and Y,                  default: [ 1  1 ]
%
% 'order'  [ nfx  nfy  ]          Number of Frames in X and Y,              default: autodetected
%
% 'space'  [ spx  spy ] | spc            Space between Frames
%                                         0 .. 1  normalized to Nx,Ny
%                                         > 1     absolut in Nx, Ny
%
% 'offset' [offx offy ] | [ Left Right Bott Top ]   Offset of Frames to Boundary
%
%
% ColorProperties
%
% 'cmap'   [ R G B ] | { fcn varargin }  ColorMap or CellArray with ColorMapFcn has first Element
%                                         for using FEVAL, example: { 'jet'  65 }
% 'clim'   [ cmin cmax ]        ColorLimit
% 'csym'    cs                  ColorSymetricValue at Center of ColorMap
%                               
% 'bcol'   [ R G B ] | ColorSpec         BackGroundColor   default: Color at SymetricValue or 'w'
% 'lcol'   [ R G B ] | ColorSpec         Color of SeperatorLines and BoundLine
% 'line'   [ lh lv lb ]                  LineWidth of Lines: [ Horiz Vert Bound ], default: [ 1  0  1 ]
%                                         0 or NaN ==> No Line
%
% 'space'  [ spx  spy ] | spc            Space between Frames
%                                         0 .. 1  normalized to Nx,Ny
%                                         > 1     absolut in Nx, Ny
%
% 'offset' [offx offy ] | [ Left Right Bott Top ]   Offset of Frames to Boundary
%                                         0 .. 1  normalized to Nx,Ny
%                                         > 1     absolut in Nx, Ny
%
% 'title'   '<Title>'  | { '<Title>'  Property Value ... }  Title 
%
% 'author'  '<Author>' | { '<Author>' Property Value ... }  Author 
% 'apos'    [ apx  apy ]   AuthorPosition, normalized to Figure,
%                            >= 0 measured from lower left
%                            <  0 measured from upper right
%
% 'label'   { 'Label1' , ... }     CellStringArray for Labels of Frames
% 'lprop'   {  Property Value ... }  TextPropertyValue-Pairs for Labels
% 'lpos'    [lpx lpy ]     LabelPositions, normalized to Frame
%                            >= 0 measured from lower left
%                            <  0 measured from upper right


%*********************************************************************
% Basic InputCheck

Nin = nargin;

if Nin < 1
   error('Input Data is missing');
end

if ~( isnumeric(d) & ( ndims(d) <= 3 ) )
    error('Data must be a max. 3-dimensional numeric.');
end

if ~chkprop(varargin)
    error('Following Inputs must be Property-Value-Pairs.');
end

%*********************************************************************
% Get Configuration

siz = cat( 2 , size(d) , 1 );
[msg,cnf] = getcnf(siz,varargin);

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Invalid Inputs.\n%s',msg))
end

%---------------------------------------------
% Make monotonic increasing

if cnf.xvec(siz(2)) < cnf.xvec(1)
   d = d(:,(siz(2):-1:1),:);
   cnf.xvec = cnf.xvec(siz(2):-1:1);
end

if cnf.yvec(siz(1)) < cnf.yvec(1)
   d = d((siz(1):-1:1),:,:);
   cnf.yvec = cnf.yvec(siz(1):-1:1);
end

%---------------------------------------------
% Check Configuration

cnf = chkcnf(siz,cnf,[min(d(:)) max(d(:))]);

%*********************************************************************
% Create Figure

[msg,cfg] = newfig(siz,cnf);

if ~isempty(msg)
    try, delete(cfg.fig), end
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Error create Figure.\n%s',msg));
end

%*********************************************************************
% Positioning

xx = cnf.xvec - cnf.xvec(1);
yy = cnf.yvec - cnf.yvec(1);

sx = cnf.order(1) * ( xx(siz(2)) + cnf.space(1) ) - cnf.space(1) + ...
     cnf.offset(1) + cnf.offset(2);

sy = cnf.order(2) * ( yy(siz(1)) + cnf.space(2) ) - cnf.space(2) + ...
     cnf.offset(3) + cnf.offset(4);

iy = ceil( ( 1 : siz(3) ) / cnf.order(1) );
ix = ( 1 : siz(3) ) - (iy-1) * cnf.order(1);

x0 = (ix-1) * ( xx(siz(2)) + cnf.space(1) ) + cnf.offset(1);
y0 = (cnf.order(2)-iy) * ( yy(siz(1)) + cnf.space(2) ) + cnf.offset(3);

if cnf.order(1) > 1
   lx = ( 1 : (cnf.order(1)-1) );
   lx = lx * ( xx(siz(2)) + cnf.space(1) ) - 1/2*cnf.space(1) + cnf.offset(1);
   lx = lx / sx;
end

if cnf.order(2) > 1
   ly = ( 1 : (cnf.order(2)-1) );
   ly = ly * ( yy(siz(1)) + cnf.space(2) ) - 1/2*cnf.space(2) + cnf.offset(2);
   ly = ly / sy;
end

x0 = x0 / sx;
y0 = y0 / sy;

xx = xx / sx;
yy = yy / sy;

%--------------------------------------------------------------
% Label and Author

x1 = x0 + xx(siz(2)) * ( cnf.lpos(1) + ( cnf.lpos(1) < 0 ) );
y1 = y0 + yy(siz(1)) * ( cnf.lpos(2) + ( cnf.lpos(2) < 0 ) );

x2 = ( cnf.apos(1) + ( cnf.apos(1) < 0 ) );
y2 = ( cnf.apos(2) + ( cnf.apos(2) < 0 ) );


ht  = { 'left'  'right' };

htl = ht{ 1 + ( cnf.lpos(1) < 0 ) };
hta = ht{ 1 + ( cnf.apos(1) < 0 ) };

vt  = { 'baseline' 'cap' };

vtl = vt{ 1 + ( cnf.lpos(2) < 0 ) };
vta = vt{ 1 + ( cnf.apos(2) < 0 ) };

%--------------------------------------------------------------
% FigurePositioning, fit to Size

pos = cnf.size([1 1],:);

pos(1,2) = pos(1,1) * sy/sx;
pos(2,1) = pos(2,2) * sx/sy;

ii = find( sum( ( pos <= cnf.size([1 1],:) ) , 2 ) == 2 );

pos = pos(ii(end),:);

figpos        = zeros(1,4);
figpos([3 4]) = pos;

figpos([1 2]) = cfg.scr([3 4])-cfg.off-figpos([3 4]);

%--------------------------------------------------------------
% PaperPosition, Fit to PaperSize

mode = { 'portrait'  'landscape' };

set( cfg.fig , 'paperunits'   , 'inches' , ...
           'paperorientation' , 'portrait' );

pap_si = get(cfg.fig,'papersize');          % [ Width  Height ]  
ppi    = get(0,'screenpixelsperinch');

img_si = pap_si - 2 * cfg.pap;          % Max. [ Width Height ]
fig_si = figpos([3 4]) / ppi;           % Act. [ Width Height ]

sc = min( img_si./fig_si , 1 );

flip = ( sc(1) < sc(2) );

mode = mode{ 1 + flip };

pap_si = pap_si([1 2]+[1 -1]*flip);
img_si = img_si([1 2]+[1 -1]*flip);

sc = min( img_si./fig_si , 1 );

pappos([3 4]) = fig_si .* min(sc);
pappos([1 2]) = ( pap_si - pappos([3 4]) ) / 2;

%*********************************************************************
% CData

is_nan = find(isnan(d));

d = min(max(d,cnf.clim(1)),cnf.clim(2));

d(is_nan) = NaN;

nc = size(cnf.cmap,1);

if isnan(cnf.csym)

   dc = ( cnf.clim(2) - cnf.clim(1) ) / nc;
   d = floor( ( d - cnf.clim(1) ) / dc ) + 1;

else

   nc2 = nc / 2;
   mc2 = mod(nc2,1);

   cl = [ cnf.clim(1)  cnf.csym  cnf.clim(2) ];

   dc = diff(cl,1,2) / nc2;
   
   ok = ~( dc == 0 );

   dc = dc + (~ok);

   d = ok(1) * ( d <  cl(2) ) .* ( floor((d-cl(1))/dc(1)) + 1 ) + ...
       ok(2) * ( d >= cl(2) ) .* ( floor((d-cl(2))/dc(2)+mc2) + 1 + nc2-mc2 );

end

%*********************************************************************
% Set

set(cfg.fig,'position',figpos, ...
            'paperorientation' , mode , ...
            'paperposition',pappos);

set(cfg.axe,'dataaspectratio',[sy sx 1]);

zz = zeros(siz(1),siz(2));

for ii = 1 : siz(3)

    if strcmp(cnf.type(1),'s')

       set( cfg.hs(ii) , 'xdata' , xx + x0(ii) , ...
                         'ydata' , yy + y0(ii) , ...
                         'zdata' , zz          , ...
                         'cdata' , d(:,:,ii)      );

    else

       dx = ( xx(siz(2)) - xx(1) ) / siz(2);
       dy = ( yy(siz(1)) - yy(1) ) / siz(1);

       set( cfg.hs(ii) , 'xdata' , xx([1 siz(2)]) + x0(ii) + [ 1  -1 ]*dx/2, ...
                         'ydata' , yy([1 siz(1)]) + y0(ii) + [ 1  -1 ]'*dy/2, ...
                         'cdata' , d(:,:,ii)      );

    end

    set( cfg.hl(ii) , 'position' , [ x1(ii) y1(ii) 0 ] , ...
                      'horizontalalignment' , htl , ...
                        'verticalalignment' , vtl );
end

if cnf.order(1) > 1
   lx = lx(ones(1,3),:);
   set(cfg.hy,'xdata',lx(:));
end

if cnf.order(2) > 1
   ly = ly(ones(1,3),:);
   set(cfg.hx,'ydata',ly(:));
end

set( cfg.ht , 'position' , [ 0.5  1-cnf.offset(4)/sy/2 0 ] , ...
              'horizontalalignment' , 'center' , ...
                'verticalalignment' , 'middle' );

set( cfg.ha , 'position' , [ x2 y2 0 ] , ...
              'horizontalalignment' , hta , ...
                'verticalalignment' , vta );

set(cfg.fig,'visible','on');

%-----------------------------------------------------------
% Check FigurePositon

drawnow

set(cfg.fig,'position', 2*figpos - get(cfg.fig,'position') );

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function cnf = chkcnf(siz,cnf,cl);

%'xvec'   , { [] } , ...
%'yvec'   , { [] } , ...
%'size'   , { [8 6]*1e2 } , ... % PixelSize    [ Width Height ]
%'scale'  , { [ 1  1 ]  } , ... % SizeScale    [ Wx  Wy ] for Data
%'order'  , { [NaN NaN] } , ... % SubOrder     [ Nx  Ny ]
%'clim'   , { [NaN NaN] } , ... % ColorLimit   [ CMin CMax ]
%'csym'   , { NaN }  , ...  % ColorSymetricValue
%'cmap'   , { [] }   , ...  % ColorMap          RGB | feval(cmap{:})
%'bcol'   , { 'w' }  , ...  % BackGroundColor   RGB | ColorSpec
%'lcol'   , { 'k' }  , ...  % LineColor         RGB | ColorSpec
%'line'   , { [1 0 2]}    , ...  % LineWidth: [ Horiz Vertc Bound ]
%'space'  , { [1 1]*.1 } , ...  % XYSpace        norm | abs
%'offset' , { [NaN NaN] } , ...  % BoundOffset [ X  Y ] | [ Left Right Bott Top ]
%'title'  , { '' }  , ...  % Title             String |  { String Prop Val ... }
%'author' , { '' }  , ...  % Author            String |  { String Prop Val ... }
%'apos'   , { [] }  , ...  % AuthorPosition    [ X  Y ]  norm
%'label'  , { {} }  , ...  % SequenceLabels    { Strings }
%'lprop'  , { {} }  , ...  % LabelProperties   { Prop Val Prop Val ... }
%'lpos'   , { [-1 1]*.01 } ); % LabelPosition(s)  [ X  Y ]  norm

%------------------------------------------
% Order

if all(isnan(cnf.order))
   if siz(3) <= 3
      cnf.order = [ siz(3) 1 ];
   else
      nx = sqrt(siz(3));
      if mod(nx,1) == 0
         cnf.order = [ nx  nx ];
      else
         nx = ceil(nx);
         ny = ceil(siz(3)/nx);
         while ny >= nx
              nx = nx + 1;
              ny = ceil(siz(3)/nx);
         end
         cnf.order = [ nx  ny ];
      end
   end
elseif isnan(cnf.order(1))
   cnf.order(1) = ceil( siz(3) / cnf.order(2) );
elseif isnan(cnf.order(2))
   cnf.order(2) = ceil( siz(3) / cnf.order(1) );
end

%------------------------------------------
% Scale

if ~( cnf.scale(1) == 1 )
     cnf.xvec = cnf.xvec(1) + ( cnf.xvec - cnf.xvec(1) ) * cnf.scale(1);
end

if ~( cnf.scale(2) == 1 )
     cnf.yvec = cnf.yvec(1) + ( cnf.yvec - cnf.yvec(1) ) * cnf.scale(2);
end

%------------------------------------------
% Space


cnf.space = ( abs(cnf.space) <= 1 ) .* cnf.space + ...
            ( abs(cnf.space) >  1 ) .* cnf.space ./ siz([2 1]);

cnf.space(1) = cnf.space(1) * ( cnf.xvec(siz(2)) - cnf.xvec(1) );
cnf.space(2) = cnf.space(2) * ( cnf.yvec(siz(1)) - cnf.yvec(1) );

%------------------------------------------
% Offset

if size(cnf.offset,2) == 2
   cnf.offset = cnf.offset([1 1 2 2]);
end

cnf.offset = ( abs(cnf.offset) <= 1 ) .* cnf.offset + ...
             ( abs(cnf.offset) >  1 ) .* cnf.offset ./ siz([2 2 1 1]);

cnf.offset([1 2]) = cnf.offset([1 2]) * ( cnf.xvec(siz(2)) - cnf.xvec(1) );
cnf.offset([3 4]) = cnf.offset([3 4]) * ( cnf.yvec(siz(1)) - cnf.yvec(1) );

jj = isnan(cnf.offset);
if any(jj)
   jj = find(jj);
   cnf.offset(jj) = cnf.space(ceil(jj/2)) / 2;
end

%------------------------------------------
% ColorLimit

acc = 1e-10;

if isnan(cnf.clim(1))
   cnf.clim(1) = min( cl(1) , cnf.clim(2)-acc );
end
if isnan(cnf.clim(2))
   cnf.clim(2) = max( cl(2) , cnf.clim(1)+acc );
end

%------------------------------------------
% ColorMap

if isempty(cnf.cmap)
   cnf.cmap = jet(256);
end

%------------------------------------------
% BackGround

if isempty(cnf.bcol)
   if isnan(cnf.csym)
      cnf.bcol = 'w';
   elseif cnf.csym <= cnf.clim(1)
      cnf.bcol = cnf.cmap(1,:);
   elseif cnf.csym >= cnf.clim(2)
      cnf.bcol = cnf.cmap(end,:);
   else
      nc = size(cnf.cmap,1);
      if mod(nc,2) == 1
         cnf.bcol = cnf.cmap(ceil(nc/2),:);
      else
         cnf.bcol = ( cnf.cmap(nc/2,:) + cnf.cmap(nc/2+1,:) ) / 2;
      end
   end
end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,cnf] = getcnf(siz,vin)

cnf = defaults(siz);

msg = cell(0,1);

cspc = colorspec;

for ii = 1 : 2 : prod(size(vin))

    p = lower(vin{ii});
    v = vin{ii+1};
    s = size(v);
    n = prod(s);

    m = '';

    is_num  = isnumeric(v);
    is_char = chkstr(v,1);
    is_cell = iscell(v);
 
    if ( n == 0 )
       p = 'none';
    end

    %-----------------------------------------------------------

    switch p

      %---------------------------------------------------------
      case { 'xvec'  'yvec' }

           is_x = strcmp(p,'xvec');
           nv   = siz(1+is_x);

           ok = ( is_num & ( 1+(nv>1) <= n ) & ( n <= nv ) );
           if ok
              v = v(:);
              d = sign(diff(v));
              ok = all( ( d == 1 ) | ( d == -1 ) );
           end

           if ~ok
               m = sprintf(' must be a monotonic Vector with max. Length %.0f',nv);
           elseif is_x
               v = v';
           end

      %---------------------------------------------------------
      case { 'bcol'  'lcol' }

           ok = ( is_num & isequal(s,[1 3]) );
           if ~ok & is_char
               ok = any(strcmp(cspc,v));
           end

           if ~ok
               m = ' must be a RGB-Tripel or valid ColorSpec.';
           end

      %---------------------------------------------------------
      case 'cmap'

           if is_cell
              ok = chkstr(v{1},1);
              if ok
                 try
                    v = feval(v{:});
                 catch
                    m = sprintf('Error evaluate "%s" using FEVAL.\n%s',p,lasterr);
                 end
              else
                 m = ' must be a CellArray with FunctionName in 1. Element.';
              end
           end

           if isempty(m)
              ok = ( ( ndims(v) == 2 ) & isnumeric(v) & ( size(v,2) == 3 ) );
              if ok
                 ok = all( ( 0 <= v(:) ) & ( v(:) <= 1 ) );
              end
              if ~ok
                  m = ' must be a RGB-Matrix with Values betweeen 0 and 1.';
              end
           end 

      %---------------------------------------------------------
      case 'type'

           ok = is_char;
           if ok
              v = lower(v(1));
              ok = any(strcmp(v(1),{'i' 's'}));
           end

           if ~ok
               m = ' must be a String, starting with <i>mage or <s>urface.';
           end

      %---------------------------------------------------------
      case { 'author' 'title' }

          if is_char
             d = getfield(cnf,p);
             d{1} = v;
             v = d;
          else
              ok = is_cell;
              if ok
                 ok = chkstr(v{1});
                 if ok & ( n > 1 )
                    ok = chkprop(v(2:n));
                 end
              end
              if ok
                 d = getfield(cnf,p);
                 d(1) = v(1);
                 if n > 1
                    v = v(:)';
                    d = cat(2,d,v(2:end));
                 end
                 d = v;
              else
                 m = ' must be a String and Property-Value-Pairs.'; 
              end
          end

      %---------------------------------------------------------
      case 'lprop'

            if ~chkprop(v)
                m = ' must be a String and Property-Value-Pairs.'; 
            else
                v = cat(2,getfield(cnf,p),v(:)');
            end

      %---------------------------------------------------------
      case { 'size'  'scale'  'clim' }
 
            is_clim = strcmp(p,'clim');

            ok = ( is_num & isequal(s,[1 2]) );
            if ok
               ok = all( isfinite(v) | ( isnan(v) & is_clim ) );
            end
 
            if ~ok
                m = ' must be a 2-Element finite numeric Vector.';
            elseif strcmp(p,'size') & ~all( v > 1 )
                   m = ' must have finite Elements larger 1.';
            elseif is_clim & ~( any(isnan(v)) | ( v(1) < v(2) ) )
                   m = ' must increasing ([Min Max]).';
             end

      %---------------------------------------------------------
      case 'order'

            ok = ( is_num & ( n <= 2 ) );
            if ok
               ok = ( ( mod(v,1) == 0 ) & ( v > 0 ) );
               ok = all( ok | isnan(v) );
            end
            if ~ok
                m = ' must have max 2 Integers larger 0 or NaN.';
            elseif n == 1
                v = cat( 2 , v , NaN );
            elseif ~any(isnan(v)) & ( prod(v) < siz(3) )
                m = sprintf(['Product of "%s" must be equal or larger then ' ...
                             'Number of Sequences %.0f.'],p,siz(3));

            end

      %---------------------------------------------------------
      case 'csym'

            ok =  ( is_num & ( n == 1 ) );
            if ok
               ok = all( isfinite(v) | isnan(v) );
            end
            if ~ok
                m = ' must be a single finite numeric or NaN.';
            end                

      %---------------------------------------------------------
      case { 'apos'  'lpos' }
 
            ok = ( is_num & isequal(s,[1 2]) );
            if ok
               ok = all(isfinite(v));
            end
            if ~ok
                m = ' must be a 2-Element finite Vector.';
            end

      %---------------------------------------------------------
      case 'label'

           [ok,v] = chkcstr(v);
           if ~ok
               m = ' must be a StringArray.';
           end

      %---------------------------------------------------------
      case 'line'

           ok = ( is_num & ( n <= 3 ) );
           if ok
              ok = ( ( v >= 0 ) & isfinite(v) );
              ok = all( ok | isnan(v) );
           end

           if ~ok
               m = ' must have max 3 finite numerics or NaN.';
           elseif n == 1
               v = v(ones(1,3));
           elseif n == 2
               v = cat( 2 , v(:)' , cnf.line(3) );
           end
 
      %---------------------------------------------------------
      case { 'space'  'offset' }

           is_off = strcmp(p,'offset');

           ok = ( is_num & ( ( n == 2 ) | ( ( n == 4 ) & is_off ) ) );
           if ok
              ok = all( isfinite(v) | isnan(v) );
           end

           if ~ok 
               if is_off
                   m = ' must have 2 or 4 finite numerics or NaN.';
               else
                   m = ' must have 2 finite numerics or NaN.';
               end
           end
   
      %---------------------------------------------------------
      case 'none'

           ok = 1;
            
      %---------------------------------------------------------
      otherwise

            m = sprintf('Invalid Property "%s".',p);

    end

    %-----------------------------------------------------------

    if ~isempty(m)
        if strcmp(m(1),' ')
           m = sprintf('Value for "%s"%s',p,m);
        end
        msg = cat(1,msg,{m});
    elseif ~isempty(v)
        cnf = setfield(cnf,p,v);
    end

end

%*********************************************************************
% Check Length of Vectors

is_img = strcmp(cnf.type(1),'i');

for p = { 'xvec'  'yvec' };

    m = '';

    v = getfield(cnf,p{1});

    is_x = strcmp(p{1}(1),'x');

    n = prod(size(v));
 
    nv = siz(1+is_x);

   ok = ( ( n == nv ) | ( ( n == 2 ) & is_img & ( nv > 1 ) ) );

   if ok & is_img & ( n == nv ) & ( nv > 1 )
      d = diff(v);
      ok = all( d-d(1) == 0 );
      if ~ok
          m = ' must be equal spaced.';
      end
   elseif ~ok
      m = sprintf(' must have Length %.0f',nv);
   end

   if ~isempty(m)
       if strcmp(m(1),' ')
          m = sprintf('Value for "%s"%s',p{1},m);
       end
       msg = cat(1,msg,{m});
   end

end

       
%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function def = defaults(siz)

tp = 'image';

at = { '' 'FontName' 'helvetica' 'FontUnits' 'points' ...
          'FontSize'  10  'Fontweight' 'bold' 'FontAngle' 'normal' };

tt = { '' 'FontName' 'helvetica' 'FontUnits' 'points' ...
          'FontSize'  14  'Fontweight' 'bold' 'FontAngle' 'normal' };

lp = at(2:end);          % LabelProp
ps = [ -1  1 ]  *0.01;   % LabelPos, AuthorPos
sp = [  1  1 ] * 0.1;    % Space

sz = [ 640 480 ];        % Size
lw = [ 1  0  1 ];        % LineWidth [H V B]

n1 = NaN;
n2 = [ NaN  NaN ];
o2 = [  1    1  ];

xv = ( 1 : siz(2) );
yv = ( 1 : siz(1) )';

def = struct( ...
              ...
'xvec'   , { xv  } , ...  % XVector
'yvec'   , { yv  } , ...  % YVector
'size'   , { sz  } , ...  % PixelSize         [ Width Height ]
'scale'  , { o2  } , ...  % SizeScale         [ Wx  Wy ] for Data
'order'  , { n2  } , ...  % SubOrder          [ Nx  Ny ]
'type'   , { tp  } , ...  % Type              'image' | 'surface'
'clim'   , { n2  } , ...  % ColorLimit        [ CMin CMax ]
'csym'   , { n1  } , ...  % ColorSymetricValue
'cmap'   , { []  } , ...  % ColorMap          RGB | feval(cmap{:})
'bcol'   , { []  } , ...  % BackGroundColor   RGB | ColorSpec
'lcol'   , { 'k' } , ...  % LineColor         RGB | ColorSpec
'line'   , { lw  } , ...  % LineWidth:        [ Horiz Vertc Bound ]
'space'  , { sp  } , ...  % XYSpace           norm | abs
'offset' , { n2  } , ...  % BoundOffset       [ X  Y ] | [ Left Right Bott Top ]
'title'  , { tt  } , ...  % Title             String |  { String Prop Val ... }
'author' , { at  } , ...  % Author            String |  { String Prop Val ... }
'apos'   , { ps  } , ...  % AuthorPosition    [ X  Y ]  norm
'label'  , { {}  } , ...  % SequenceLabels    { Strings }
'lprop'  , { lp  } , ...  % LabelProperties   { Prop Val Prop Val ... }
'lpos'   , { ps  } );     % LabelPosition(s)  [ X  Y ]  norm



%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function cspec = colorspec

cspec = { ...
    'black'
    'blue'
    'red'
    'green'
    'cyan'
    'magenta'
    'yellow'
    'white'
    'k'
    'b'
    'r'
    'g'
    'c'
    'm'
    'y'
    'w'  };

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,cfg] = newfig(siz,cnf);

off = [ 10  60 ];   % PixelOffset [ Right Top ]
pap = 1;            % PaperOffset [inch]

scr_uni = get(0,'units');      set(0,'units','pixels')
scr_si  = get(0,'ScreenSize'); set(0,'units',scr_uni);

figpos        = zeros(1,4);
figpos([3 4]) = cnf.size;

figpos([1 2]) = scr_si([3 4])-off-figpos([3 4]);


msg = cell(0,1);

%-----------------------------------------------------------

 fig = figure( 'paperunits'         , 'inches'   , ...
               'paperorientation'   , 'portrait' , ...
               'units'              , 'pixels'   , ...
               'position'           ,  figpos    , ...
               'color'              , cnf.bcol   , ...
               'numbertitle'        , 'on'       , ...
               'menubar'            , 'none'     , ...
               'toolbar'            , 'none'     , ...
               'name'               , ''         , ...
               'colormap'           , cnf.cmap   , ...
               'visible'            , 'off'      , ...
               'handlevisibility'   , 'on'       , ...
               'render'             , 'zbuffer'  , ...
               'createfcn'          , ''         , ...
               'resize'             , 'off'             );

%-----------------------------------------------------------

lp = lineprop(cnf.line(3));

axe = axes( 'parent'   , fig            , ...
            'units'    , 'normalized'   , ...
            'position' , [ 0  0  1  1 ] , ...
            'xlim'     , [ 0  1 ]       , ...
            'ylim'     , [ 0  1 ]       , ...
            'xtick'    , []             , ...
            'ytick'    , []             , ...
            'color'    , 'none'         , ...
            'xcolor'   , cnf.lcol       , ...
            'ycolor'   , cnf.lcol       , ...
            'xgrid'    , 'off'          , ...
            'ygrid'    , 'off'          , ...
            'xdir'     , 'normal'       , ...
            'ydir'     , 'normal'       , ...
            lp{:}                       , ...
            'box'      , 'on'           , ...
            'layer'    , 'top'          , ...
            'view'     , [ 0  90 ]      , ...
         'clipping'    , 'off'          , ...
     'dataaspectratio' , [ 1  1  1 ]    , ...
     'dataaspectratiomode' , 'auto'     , ...
            'visible'  , 'on'           , ...
    'handlevisibility' , 'callback'         );

n = siz(3);

%-----------------------------------------------------------
% Surfaces

hs = zeros(n,1);

n4 = NaN * ones(2,2);

if strcmp(cnf.type(1),'s');

   hs(1) = surface( 'parent' , axe , ...
             'xdata'  , [0 1]  , ...
             'ydata'  , [0 1]' , ...
             'cdata'  , n4     , ...
       'cdatamapping' , 'direct'  , ...
          'facecolor' , 'interp' , ...
          'edgecolor' , 'none' , ...
           'clipping' , 'on'   , ...
           'visible'  , 'on'               );

else

   hs(1) = image( 'parent' , axe , ...
             'xdata'  , [0 1]  , ...
             'ydata'  , [0 1]' , ...
             'cdata'  , n4     , ...
       'cdatamapping' , 'direct' , ...
           'clipping' , 'on'   , ...
           'visible'  , 'on'               );

end

for ii = 2 : n
    hs(ii) = copyobj(hs(1),axe);
end

%-----------------------------------------------------------
% XLine

yy = NaN * zeros(3,cnf.order(2)-1);
xx = [ 0 ; 1 ; NaN ] * ones(1,cnf.order(2)-1);

lp = lineprop(cnf.line(1));

hx = line( 'parent' , axe    , ...
           'xdata'  , xx(:)  , ...
           'ydata'  , yy(:)  , ...
           lp{:}             , ...
           'marker' , 'none' , ...
           'color'  , cnf.lcol   , ...
          'visible' , 'on'   , ...
          'clipping' , 'on'      );

%-----------------------------------------------------------
% YLine

xx = NaN * zeros(3,cnf.order(1)-1);
yy = [ 0 ; 1 ; NaN ] * ones(1,cnf.order(1)-1);

lp = lineprop(cnf.line(2));

hy = line( 'parent' , axe    , ...
           'xdata'  , xx(:)  , ...
           'ydata'  , yy(:)  , ...
           lp{:}             , ...
           'marker' , 'none' , ...
           'color'  , cnf.lcol   , ...
          'visible' , 'on'   , ...
          'clipping' , 'on'      );

%-----------------------------------------------------------
% Labels

hl = zeros(n,1);

nl = prod(size(cnf.label));

m = '';
try
   hl(1) = text( cnf.lprop{:} , 'parent' , axe );
catch
   m = lasterr;
end

if ~isempty(m)
   msg = cat(1,msg,{sprintf('Invalid LabelProperties.\n%s',m)});
else
  for ii = 2 : n
      hl(ii) = copyobj(hl(1),axe);
      if ( nl >= ii )
         set(hl(ii),'string',cnf.label{ii});
      end
  end
  if nl >= 1
     set(hl(1),'string',cnf.label{1});
  end
end
 
%-----------------------------------------------------------
% Title

ht = [];
try
   ht = text( 'string' , cnf.title{:} , 'parent' , axe );
catch
   msg = cat(1,msg,{sprintf('Invalid TitleProperties.\n%s',lasterr)});
end
   
%-----------------------------------------------------------
% Author

ha = [];
try
   ha = text( 'string' , cnf.author{:} , 'parent' , axe );
catch
   msg = cat(1,msg,{sprintf('Invalid AuthorProperties.\n%s',lasterr)});
end

%------------------------------------------------------------

cfg = struct( 'fig' , { fig } , ...
              'axe' , { axe } , ...
              'hs'  , { hs }  , ...
              'hx'  , { hx }  , ...
              'hy'  , { hy }  , ...
              'hl'  , { hl }  , ...
              'ht'  , { ht }  , ...
              'ha'  , { ha }  , ...
              'off' , { off } , ...
              'pap' , { pap } , ...
              'scr' , { scr_si }      );

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function lp = lineprop(lw)

ls = '-';

lp = 'linestyle';

if isnan(lw) | ( lw == 0 )
   lp = 'visible';
   ls = 'off';
   lw = 0.5;
end

lp = { lp , ls , 'linewidth' , lw }; 
 
%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = chkprop(v)

ok = 1;
if isempty(v)
   return
end

ok = iscell(v);
if ok
   n = prod(size(v));
   ok = ( mod(n,2) == 0 );
   if ok
      ok = chkcstr(v(1:2:n-1));
   end
end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
   str = cellstr(str);
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

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

   