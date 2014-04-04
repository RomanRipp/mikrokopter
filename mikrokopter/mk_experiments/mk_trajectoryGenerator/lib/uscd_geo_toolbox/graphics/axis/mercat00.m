function newhandle = mercat00(axe,ext,ext0,maxint)

% MERCAT00  Set Axes to MercatorProjection 
%
% MERCAT00( Axe , Ext_New , Ext_Old , NTickInt )
%  
% Axe   AxeHandle to set in Mercator
%
% Ext_New   Scale for new Projection
% Ext_Old   Scale for old Projection,
%            []  |  0  if no Mercator before ! 
% NTickInt  Number of TickIntervalls per Axe
%
%   
% MERCAT00 changes the  !!!  ydata  !!! of Axes Children
%
%           and following AxesProperties:
%                            ylim       
%                            ytick
%                            yticklabels
%                            xtick
%                            xticklabels
%                            'dataaspectratio' , [1 1 1]
%
% Manual AxesTicks will not be changed, 
%  if No Input "NTickInt" is given !!!
%
%  all vectors are handled as RowVector
%
% NewHandle = MERCAT00( ... ) gives a vector of Handles of new
%  created graphic objects (axes, patches or lines ... )
%
%
% The Calculation for the Projection of the Y-Data:
%
%%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% function  y = deg2merc(y,ext)
%
% y = 180/pi * log(abs(tan(pi/4 + y*ext*pi/180/2))) / ext;
%
%%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% function  y = merc2deg(y,ext)
%
% y = 2*180/pi * ( atan(exp(y*ext*pi/180)) - pi/4 ) / ext; 
%
%%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%
%
%  WARNING:  MERCAT00 didn't work correctly with TYPE IMAGE
%            Scale transformed YData (deg2merc(y,ext)) equal spaced before. 
%



Nin = nargin;
Nout = nargout;

if Nin < 1
 axe = gca;
end

if Nin < 2
 ext = 1;
end

if Nin < 3
 ext0 = [];
end

is_maxint = 1;
if Nin < 4 
 is_maxint = 0;
 maxint = 5;
end


if isempty(ext0)
  ext0 = 0;
end


newhandle = [];


is_set = 1;
 
is_set = ( ext ~= ext0 );



fig = get(axe,'parent');

xlim = get(axe,'xlim');
ylim = get(axe,'ylim');

ylim = merc2deg(ylim,ext0);




%*****************************************************************
% Prepare Ticks


 tick_int = [ 90 ( [ 60 30 20 10  5  2  1 ]/1         ) , ...
                 ( [    30 20 10  5  2  1 ]/60        ) , ...
                 ( [    30 15 12  6       ]/3600      ) , ...
                 ( [    30 15 12  6       ]/3600/10   ) , ...
                 ( [    30    12  6       ]/3600/100  )       ];

 ii = [ NaN  NaN ];

 lims = [ xlim ; ylim ];

 [hilf,ii(1)] = min( abs( diff(lims(1,:)) / maxint - tick_int) );
 [hilf,ii(2)] = min( abs( diff(lims(2,:)) / maxint - tick_int) );

 [ hilf , kk ] = max(tick_int(ii));

  ok = find( ( rem( tick_int , tick_int(ii(kk)) ) == 0 )   &  ...
               (  tick_int >= tick_int(ii(kk))  ) ); 

 [hilf,jj] = min( abs( diff(lims(3-kk,:)) / maxint - tick_int(ok)) );

 ii(3-kk) = ok(jj);


 xticks = tick_int(ii(1)) * ( ceil(lims(1,1)/tick_int(ii(1))) : ...
                             floor(lims(1,2)/tick_int(ii(1)))       );

 yticks = tick_int(ii(2)) * ( ceil(lims(2,1)/tick_int(ii(1))) : ...
                             floor(lims(2,2)/tick_int(ii(1)))       );

if isempty(xticks) 
 xticks = tick_int(ii(1)) * (floor(lims(1,1)/tick_int(ii(1))) : ...
                              ceil(lims(1,2)/tick_int(ii(1)))       );
 if size(xticks,2) == 2
   xticks = lims(1,:);
 end     
end


if isempty(yticks) 
 yticks = tick_int(ii(1)) * (floor(lims(2,1)/tick_int(ii(1))) : ...
                              ceil(lims(2,2)/tick_int(ii(1)))       );
 if size(yticks,2) == 2
   yticks = lims(2,:);
 end     
end


if ~is_maxint

  if strcmp(get(axe,'xtickmode'),'manual')
    xticks = get(axe,'xtick');
  end

  if strcmp(get(axe,'ytickmode'),'manual')
    yticks = merc2deg(get(axe,'ytick'),ext0);
  end

end


%*****************************************************************
% Prepare TickLabels

if any( abs( xticks-round(xticks) ) > 1e-6 )  |  ...
   any( abs( yticks-round(yticks) ) > 1e-6 )
   not_flag = 'sexges';  % sexages
else
   not_flag = 'deg';  % degree
end


%----------------------------------------------------
% X


if isempty(xticks)

  xtickl = '';

else

   %  -->  [ -180 .. 180 )
   xtick0 = xticks - 360 * floor( (xticks+180) / 360 );

   xtick0(end) = xtick0(end) + ( 180 - xtick0(end) ) * ...
                               ( xtick0(end) == -180 );

   xtickl = geolabel(xtick0,'lon',not_flag);

end
 
%----------------------------------------------------
% Y

if isempty(yticks)

  ytickl = '';

else

   ytickl = geolabel(yticks,'lat',not_flag);

end


%***********************************************************************
% Set Y-Data


if is_set

 shh = get(0,'showhiddenhandles');
       set(0,'showhiddenhandles','on');

 children = get(axe,'children');

       set(0,'showhiddenhandles','off');

 for child = children(:)'

   if ~strcmp(get(child,'type'),'text')

      y = merc2deg(get(child,'ydata'),ext0);

      set( child , 'ydata' ,  deg2merc(y,ext) )

   else

      uni = get(child,'units');
      set(child,'units','data');

      textpos=get(child,'position');
 
      textpos(2) = merc2deg(textpos(2),ext0);

      textpos(2) = deg2merc(textpos(2),ext);

      set(child,'position',textpos, ...
                'units',uni)

   end
   % Text

 end
 % child

end
% is_set



ylimM  = deg2merc( ylim ,ext);
ytickM = deg2merc(yticks,ext);

set(axe,  'tickdir'     , 'out'     , ...
          'dataaspectratio' , [ 1  1  1 ]   , ...
          'box'         , 'on'      , ...
          'xlim'        , xlim      , ...
          'xtick'       , xticks    , ...
          'xticklabel'  , xtickl    , ...
          'ytick'       , ytickM    , ...
          'ylim'        , ylimM     , ...
          'yticklabel'  , ytickl                       )  






ticklength = get(axe,'ticklength');
linewidth  = get(axe,'linewidth') ; 


   jj = find( xticks < xlim(1)  |  xticks > xlim(2) );
   xticks(jj) = [];

   jj = find( ytickM < ylimM(1)  |  ytickM > ylimM(2) );
   yticks(jj) = [];

   dxt = diff(xticks);
   dyt = diff(yticks);

   dok_x = isempty(dxt);
   if ~dok_x
     dok_x = ( max(dxt)-min(dxt) < 1e-3*mean(dxt) );
   end

   dok_y = isempty(dyt);
   if ~dok_y
     dok_y = ( max(dyt)-min(dyt) < 1e-3*mean(dyt) );
   end

 if  dok_x  &  dok_y


  [mdiff,kk]= max( [ diff(xlim), ...
                     diff(ylimM)         ]);

   bwidth = 0.010*mdiff; % BorderWidth
 
   [upp,axepos,figpos] = ppunit(axe);

   papuni = get(fig,'paperunits');
   set(fig,'paperunits','inches');
   pappos = get(fig,'paperposition');
   set(fig,'paperunits',papuni);

   % Minimum 1.8 mm  == 10 Points (150 dpi)
   bmin1 = 1.8/25.4 / (pappos(4)*axepos(4)/figpos(4)) * diff(ylimM); 
   
   % Minimum 5 Pixels
   bmin2 = 4 / axepos(4) * diff(ylimM);        % Minimum 5 Pixels
   
   bmin = max([ bmin1  bmin2 ]);

   bwidth = bwidth + (bmin-bwidth)*( bwidth < bmin );


   mdxt = [];
   if ~isempty(dxt)
    mdxt = mean(dxt);
   end
   mdyt = [];
   if ~isempty(dyt)
    mdyt = mean(dyt);
   end

   intv = [ mdxt ; mdyt ];

   intv = sort([ 4*intv ; 2*intv ; intv ; intv/2 ]);

      jj = find( intv - 3*bwidth  >= 0 ) ;

  
   if ~isempty(jj)

      intv = intv(jj(1));

    % Inner, White Border

    y_for_x_0 = ylimM([1 1 2 2]) + [ -0.5  0    0    0.5 ]*bwidth;
    x_for_y_0 =  xlim([1 1 2 2]) + [ -0.5  0    0    0.5 ]*bwidth;
    y_for_x   = ylimM([1 1 2 2]) + [ -1   -0.5  0.5  1   ]*bwidth;
    x_for_y   =  xlim([1 1 2 2]) + [ -1   -0.5  0.5  1   ]*bwidth;

    xx_0 = xlim([1 2]);
    yy_0 = ylimM([1 2]);

    y_x_0 = [ y_for_x_0([1 1 2 2])'*ones(1,length(xx_0)-1) , ...
              y_for_x_0([3 3 4 4])'*ones(1,length(xx_0)-1) , ...
              y_for_x([1 1 2 2])'*ones(1,length(xx_0)-1) , ...
              y_for_x([3 3 4 4])'*ones(1,length(xx_0)-1)    ];


    xx_0 = xx_0([1;2;2;1]); xx_0 = xx_0(:);
    xx_0 = xx_0 * ones(1,4);

    x_y_0 = [ x_for_y_0([1 1 2 2])'*ones(1,length(yy_0)-1) , ...
              x_for_y_0([3 3 4 4])'*ones(1,length(yy_0)-1) , ... 
              x_for_y([1 1 2 2])'*ones(1,length(yy_0)-1) , ...
              x_for_y([3 3 4 4])'*ones(1,length(yy_0)-1)        ];
              
    yy_0 = yy_0([1;2;2;1]); yy_0 = yy_0(:);
    yy_0 =  yy_0 * ones(1,4);


   % Outher, Black/White Border

    xx = [  xlim(1) ...
     fliplr([min(xticks):-intv:xlim(1)]) ...
             min(xticks):intv:max(xticks)  ...
            [max(xticks):intv:xlim(2)] ...
           xlim(2)       ];

    xx(find(abs(diff(xx))<1e-10))=[];

    
    yy = [  ylim(1)            ...
      fliplr([min(yticks)-intv:-intv:ylim(1)]) ...
           [min(yticks):intv:max(yticks)] ...
           [max(yticks):intv:ylim(2)] ylim(2) ];

    yy(find(abs(diff(yy))<1e-10))=[];

    yy = deg2merc(yy,ext);


    nx = size(xx,2);
    ny = size(yy,2);

    xdif = abs(diff(xx([1 2 nx-1 nx])));
    ydif = abs(diff(yy([1 2 ny-1 ny])));

    % first and last segment to small
    xsm = [xdif([1 3]) < bwidth ];
    ysm = [ydif([1 3]) < bwidth ];


    y_x = [ y_for_x([1 1 2 2])'*ones(1,nx-1) , ...
            y_for_x([3 3 4 4])'*ones(1,nx-1)     ];
    xx = xx([(1:nx-1);(2:nx);(2:nx);(1:nx-1)]);
    xx = permute(xx,([ 1  2 ] + [ 1  -1 ]*(size(xx,2)~=(nx-1))) );
    xx = [ xx  xx ]; 
  
    x_y = [ x_for_y([1 1 2 2])'*ones(1,ny-1) , ...
            x_for_y([3 3 4 4])'*ones(1,ny-1)     ];
    yy = yy([(1:ny-1);(2:ny);(2:ny);(1:ny-1)]);
    yy = permute(yy,([ 1  2 ] + [ 1  -1 ]*(size(yy,2)~=(ny-1))) );
       
    yy = [ yy  yy ]; 

    hold_state = get(axe,'nextplot');
    set(axe,'nextplot','add')

    cc = [   get(axe,'xcolor') ; ...
           1-get(axe,'xcolor')            ];  % Colors

    zz = 10 * max( [ max(get(axe,'zlim'))  max(get(axe,'clim')) ] );

    kk = 2;
    newhandle = ...
     patch( 'xdata' , [xx_0 x_y_0] , ...
            'ydata' , [y_x_0 yy_0] , ...
            'zdata' , zz+0*[xx_0 x_y_0] ,  ...
           'facecolor' , cc(kk,:) , ...
           'edgecolor',cc(1,:), ...
           'linewidth',linewidth, ...
           'clipping','off',...
           'tag','MERCATOR_BORDER' , ...
           'parent',axe);              

    for kk = [ 1 ]
    ff = 3-2*kk;  % Faktor [ 1  -1 ]
    ix = [ [ff*xsm(1)+kk:2:nx-1]   [ff*xsm(1)+kk+nx-1:2:2*(nx-1)] ];
    iy = [ [ff*ysm(1)+kk:2:ny-1]   [ff*ysm(1)+kk+ny-1:2:2*(ny-1)] ];

    newhandle = [ newhandle ; ...
     patch(   'xdata', [xx(:,ix) x_y(:,iy)] , ...
              'ydata',[y_x(:,ix) yy(:,iy)]  , ...
              'zdata' , zz+0*[xx(:,ix) x_y(:,iy)] , ...
           'facecolor' , cc(kk,:) , ...
           'edgecolor',cc(1,:), ...
           'linewidth',linewidth, ...
           'clipping','off',...
	   'tag','MERCATOR_BORDER' , ...
           'parent',axe)           ];
    end
    % kk 

   % 4 Corner's


       y_for_x_1 = ylimM([1 1 2 2]) + [-1 0  0  1]*bwidth;
       x_for_y_1 = xlim([1 1 2 2]) + [-1 0  0  1]*bwidth;

       xx_1 = x_for_y_1( [ 1 2 2 1 ; ...
                           3 4 4 3 ; ...
                           3 4 4 3 ; ...
                           1 2 2 1   ]'    );   % !!! <'> !!!
       yy_1 = y_for_x_1( [ 1 1 2 2 ; ...
                           1 1 2 2 ; ...
                           3 3 4 4 ; ...
                           3 3 4 4   ]'    );   % !!! <'> !!!
       if any(size(xx_1)==1)
       % Matlab4
        xx_1 = reshape(xx_1(:),4,4);
        yy_1 = reshape(yy_1(:),4,4);
       end

       newhandle = [ newhandle ; ...
        patch('xdata', xx_1 , ...
              'ydata', yy_1 , ...
              'zdata' , zz+0*xx_1 , ...
           'facecolor' , cc(2,:) , ...
           'edgecolor',cc(1,:), ...
           'linewidth',linewidth, ...
           'clipping','off',...
           'tag','MERCATOR_BORDER' , ...
           'parent',axe)            ];


    ticklength = [1.5*min(bwidth./[diff(xlim) diff(ylimM)])  , ...
                  ticklength(2)];

    set(axe,'nextplot',hold_state, ...
                   'ticklength',ticklength)
    


    % set TitleYPosition

    tt  = get(axe,'title');
    uni = get(tt,'units');
          set(tt,'units','normalized');

    pp    = get(tt,'position');
    pp(2) = 1 + 2*ticklength(1);
 
    set(tt,'position',pp);

    set(tt,'units',uni)


   end
 end




if Nout==0
 newhandle = [];
end 



%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function  y = deg2merc(y,ext)

if ext ~= 0

  y = 180/pi * log(abs(tan(pi/4 + y*ext*pi/180/2))) / ext;

end


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function  y = merc2deg(y,ext)
 
if ext ~= 0

 y = 2*180/pi * ( atan(exp(y*ext*pi/180)) - pi/4 ) / ext; 

end


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function [upp,axesize,figpos] = ppunit(axe);

% PPUNIT  Returns AxesUnitPerPixel  and AxesPixelSize
%
%  [ UnitPerPixel, AxePixelSize, FigurePixelSize ] = PPUNIT( AxesHandle )
%
%  UnitPerPixel = [ XUnitsPerPixel  YUnitsPerPixel ] ;
%  PixelSize    = [ PixelLeft  PixelBottom  PixelWidth  PixelHight ];
%
 

fig = get(axe,'parent');

 
fig_uni = get(fig,'units');
set(fig,'units','pixels')
figpos = get(fig,'position');
set(fig,'units',fig_uni);

axe_uni = get(axe,'units');
set(axe,'units','normalized')
axepos = get(axe,'position');
set(axe,'units',axe_uni);

dx = diff(get(axe,'xlim'));
dy = diff(get(axe,'ylim'));


axesize = axepos.*figpos([3 4 3 4]);

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

   ii = find( sum( ( pos <= [1;1]*axesize([3 4]) ) , 2 ) == 2 );

   pos = pos(ii(1),:);
 
   axesize([1 2]) = axesize([1 2])+axesize([3 4])/2-pos/2;
   axesize([3 4]) = pos;

end

 upp = [ dx  dy ] ./ axesize([3 4]) ; % UnitPerPixel

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = geolabel(tick,ini,not);

% GEOLABEL  Transform Geographical Coordinates into String
%
% Label = GEOLABEL( X , INI , NotationFlag )
%
%  X             Vector for Coordinates
%  INI           'lon'      |  'lat'
%  NotationFlag  'sexages'  |  'deg'
%
 
if isempty(tick)
  str = '';
  return
end

lab = [ 'WE' ; 'SN' ];

lab = lab( 1 + strcmp(ini,'lat') , : );

if nargin < 3
  if any( abs( tick-round(tick) ) > 1e-6 )
    not = 'deg';
  else
    not = 'sexages';
  end
end


   grad =   fix(tick);
   minu =   fix(  60*(tick-grad));
   secu = round(3600*(tick-grad-minu/60));

    minu_secu =  fix(secu/60);
   secu = secu - 60*minu_secu;
   minu = minu + minu_secu;

    grad_minu =  fix(minu/60);
   minu = minu - 60*grad_minu;
   grad = grad + grad_minu;

   dec_secu = [];
   sec_form = '';

  if any(diff(1e3*grad+minu)==0)
   dec_secu = round(1e2*( 60*(tick-grad) - minu)) ;
   sec_form = '.%2.2d';
  end
  if ~isempty(dec_secu)
   if any(diff(1e3*grad+minu+1e-3*dec_secu)==0)
    dec_secu = round(1e3*( 60*(tick-grad) - minu)) ;
    sec_form = '.%3.3d';
   end
  end

  is_dec = ~isempty(dec_secu);


  nt = prod(size(tick));

  str = cell(nt,1);
  str(:) = { '' };

  for ii = 1 : nt

     if strcmp(not,'deg')

       str{ii} = [ sprintf('%.0f',abs(grad(ii))) char(176)  ];

     else

       str{ii} = [ sprintf('%.0f',abs(grad(ii))) char(176) , ...
                   sprintf('%2.2d',abs(minu(ii)))           , ...
                   sprintf(sec_form,abs(dec_secu(ii*(1:is_dec)))) char(39)  ];

     end

  end


  for ii = [ 1  nt(1:(nt>1)) ]

    i01 = ( ( tick(ii) > 0 )  |  ( ( ii == 1 )  &  ( tick(ii) == 0 ) ) );
   
    str{ii} = [ str{ii}  lab( 1 + i01 ) ];

  end
