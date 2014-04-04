function c = tileimg(axe,siz,off,fcn)

% TILEIMG Extract a tile of an Image from an Axes
%
% C =  TILEIMG( AxesHandle , TileSize , Offset , Fcn )
%
% Size of Tile: [ Width Height ]
% PixelOffset, default: 10px
% Fcn   Function to evaluate before GETFRAME
%
% C     extracted image (tile) [ SIZ(1) x SIZ(2) x 3 ]
%
% see also: GETFRAME
%

Nin = nargin;

if Nin == 0
   axe = [];
end

if Nin < 2
   siz = [];
end

if Nin < 3
   off = [];
end

if Nin < 4
   fcn = {};
end

if isempty(axe)
   fig = get(0,'currentfigure');
   if ~isempty(fig)
       axe = get(fig,'currentaxes');
   end
   if isempty(axe)
       return
   end
end

%--------------------------------------------
% Defaults

if isempty(off)
   off = 10;
end

if isempty(siz)
   siz = [ NaN  NaN ];
end

if isempty(fcn)
   fcn = {};
elseif ~iscell(fcn)
   fcn = {fcn};
end

msg = {};

if ~( isnumeric(siz) & ( prod(size(siz)) == 2 ) )
    msg = cat(1,msg,{'Size must e a 2-element numeric.'});
end

if ~( isnumeric(off) & ( prod(size(off)) == 1 ) )
    msg = cat(1,msg,{'Offset must e a 2-element numeric.'});
end

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Invalid Inputs.\n%s',msg));
end

%--------------------------------------------
% Check Size

fig = get(axe,'parent');

[upp,pos] = ppunit(axe);

pos = pos([3 4]);

xl = get(axe,'xlim'); dx = diff(xl);
yl = get(axe,'ylim'); dy = diff(yl);

dxy = [ dx  dy ];

isz = isnan(siz);

if all(isz)
   siz = [ dx  dy ]; 
   if ~strcmp(get(axe,'dataaspectratiomode'),'auto')
       siz = pos .* max(siz./pos);
   end 
elseif any(isz)
   ii = find(isz);
   jj = 3 - ii;    % Opposite
   siz(ii) = siz(jj) * pos(ii) / pos(jj);
end

siz = ceil(siz);

if any( prod(siz) > 12000*12000 )
   error('Image to large');
end

%----------------------------------------------
% Save defaults

fprp = { 'units' 
         'position'
         'visible'   };
        
aprp = { 'dataaspectratio'
         'dataaspectratiomode'
         'units'
         'position'
         'ydir'
         'xlim' 
         'ylim'  
         'xlimmode'
         'ylimmode'
         'visible'     };


fprp = fprp(:,[1 1]);
for ii = 1 : size(fprp,1)
    fprp{ii,2} = get(fig,fprp{ii,1});
end

aprp = aprp(:,[1 1]);
for ii = 1 : size(aprp,1)
    aprp{ii,2} = get(axe,aprp{ii,1});
end

yd = get(axe,'ydir');

set(axe,'dataaspectratiomode' , 'auto', ...  
        'units'    , 'normalized' , ... 
        'position' , [ 0 0 1 1 ]  , ...
        'xlimmode' , 'manual'     , ...
        'ylimmode' , 'manual'     , ...
        'ydir'     , 'reverse' , ...
        'visible'  , 'off'              );

%----------------------------------------------
% Prepare

uni = get(0,'units'); set(0,'units','pixels');
ssi = get(0,'screensize'); set(0,'units',uni);

mrg = 60; % FigureMargin

int = ssi([3 4]) - 2*off - 2*mrg;

nn = ceil(siz./int);

x = ( 1 : siz(1) );
y = ( 1 : siz(2) );

c = uint8(zeros([siz([2 1]) 3]));

r = [ dx  dy ] ./ siz;

%----------------------------------------------
% Loop
% 1,keyboard
for i1 = 1 : nn(1)

    j1 = [ 1  int(1) ] + (i1-1) * int(1);
    j1 = min(j1,siz(1));
    n1 = j1(2) - j1(1) + 1;

    xlm = x(j1) + [ -1  1 ] * ( off + 0.5 );
    dxl = diff(xlm);
    xlm = xlm * r(1) + xl(1);

    k1 = j1 + off;
   
    j1 = ( j1(1) : j1(2) );
    k1 = ( k1(1) : k1(2) ) - j1(1) + 1;

    set(axe,'xlim',xlm)

    for i2 = 1 : nn(2)

        j2 = [ 1  int(2) ] + (i2-1) * int(2);
        j2 = min(j2,siz(2));
        n2 = j2(2) - j2(1) + 1;

        ylm = y(j2) + [ -1  1 ] * ( off + 0.5 );
        dyl = diff(ylm);
        ylm = ylm * r(2) + yl(1);

        k2 = j2 + off;
   
        j2 = ( j2(1) : j2(2) );
        k2 = ( k2(1) : k2(2) ) - j2(1) + 1;

fprintf(1,'%f ',[xlm ylm]);fprintf(1,'   ');
fprintf(1,'%5d ',[dxl dyl]);fprintf(1,'\n');

if 1
        set(axe,'ylim',ylm)

        set(fig,'pos',[ mrg mrg dxl dyl ] );

        if ~isempty(fcn)
            try
                feval(fcn{:});
            catch
                warning(lasterr);
                break
            end
        end

        figure(fig), drawnow, k = getframe(fig);
        
        c(j2,j1,:) = k.cdata(k2,k1,:);
end

       fprintf(1,'%5d ',j1([1 end]),j2([1 end]));
       fprintf(1,'\n');

    end

end


if strcmp(yd,'normal')
   c = c(siz(2):-1:1,:,:);
end

%----------------------------------------------
% Restore Settings

for ii = 1 : size(aprp,1)
    set(axe,aprp{ii,1},aprp{ii,2});
end

for ii = 1 : size(fprp,1)
    set(fig,fprp{ii,1},fprp{ii,2});
end
