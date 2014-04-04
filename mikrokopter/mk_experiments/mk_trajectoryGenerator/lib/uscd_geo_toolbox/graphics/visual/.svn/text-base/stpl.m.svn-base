function [fig,axe,axeT] = stpl(t,u,v,varargin)

% STPL  StickPlot of U-V-TimeSeries
%
% [Fig,Axe] = STPL(T,U,V,[BaseDate],[Unit+Amplitude*i],[Labels],[Layout])
%
% T   TimeVectors, N-Element CellArray or Vector
%      Values in dezimal days since BaseDate
%
% U   U-Current [cm/s], N-Element CellArray or N Columns / N Rows
% V   V-Current [cm/s], N-Element CellArray or N Columns / N Rows
%
% BaseYear or BaseDate  for T
%
% optional:
%
% Unit      Lenght for Unit-Arrow
%
% Amplitude Amplitude for single TimeSeries
%            2*Amplitude*Unit == Hight of Y-Axes
%
% Labels    Labels for TimeSeries, N-Element CellArray
%
% Layout    Mode for Layout, Modes defined at begin of this File
%
% Layout-Parameters can be edited at the Begin of this File!
%
% 
%------------------------------------------------------------------
%
% see also:  STPL_ADCP
% required:  TIMAXIS, WYSIWYG, WYGIWYS 
%
%------------------------------------------------------------------
% Example:
%
% scl = 0.001;
% phi = 1 - scl;
%%  Rrx = rrx * sqrt(1-phi^2);  % Wichtung
%
% u = randn(1000,1) / sqrt(1-phi^2);
% v = randn(1000,1) / sqrt(1-phi^2);
% u0 = u; v0 = v;
% for ii = 1 : ceil(1/scl)
%
%  u = phi * u + randn(1000,1);
%  v = phi * v + randn(1000,1);
%
% end
%
% required m-files:  WYSIWYG
%                    TIMEAXIS
%

Nin = nargin;

if Nin < 3
   error('Not enough InputArguments');
end

[msg,base,unit,ampl,labels,mode] = checkin(varargin,3);

%**********************************************
% Layout-Setup
%**********************************************

% Defaults, see below

modes = { 'large' 'normal' 'wide' 'screen' };

if isempty(mode)
   mode = modes{1};
end

switch mode

   %**********************************************
   case { 'large' 'tall' }
   %**********************************************

     % Borders at Paper [ Left  Right  Bottom Top ]
     %  ignored if NonEmpty FigPos

     PaperBorder = [ 1  0.8  1.2   1.2 ];
  
     FigPos = [];   % if empty use PaperBorder
                   
     % Normalized AxesPosition

     AxePos = [ 0.00  0.00  1.00  1.00 ];

     % Normalized Space of Axes in Figure

     Xoff = [ 0.08   0.15 ];   % [ Left Right  ]
     Yoff = [ 0.08   0.08 ];   % [ Top  Bottom ]

     AxeFontSize   = 10;
     LabelFontSize = 8;

     TickLabelIntervall = 8;   % Intervall for TIMEAXIS

   %**********************************************
   case 'normal'
   %**********************************************

     % Borders at Paper [ Left  Right  Bottom Top ]
     % ignored if NonEmpty FigPos

     PaperBorder = [ 1  1.0  1.5   1.5 ]; 
  
     FigPos = [];   %  if empty use PaperBorder

     % Normalized AxesPosition

     AxePos = [ 0.00  0.00  1.00  1.00 ];

     % Normalized Space of Axes in Figure

     Xoff = [ 0.06   0.12 ];   % [ Left Right  ]
     Yoff = [ 0.08   0.10 ];   % [ Top  Bottom ]

     AxeFontSize   = 10;
     LabelFontSize = 08;

     TickLabelIntervall = 8;   % Intervall for TIMEAXIS

   %**********************************************
   case 'wide'
   %**********************************************

     % Borders at Paper [ Left  Right  Bottom Top ]
     % ignored if NonEmpty FigPos

     PaperBorder = [ 1  1.0  1.2   1.2 ];
  
     FigPos = [];

     % Normalized AxesPosition

     AxePos = [ 0.00  0.00  1.00  1.00 ];

     % Normalized Space of Axes in Figure

     Xoff = [ 0.05   0.12 ];   % [ Left Right  ]
     Yoff = [ 0.12   0.08 ];   % [ Top  Bottom ]

     AxeFontSize   = 12;
     LabelFontSize = 10;

     TickLabelIntervall = 12;   % Intervall for TIMEAXIS

   %**********************************************
   case 'screen'
   %**********************************************

     % Borders at Paper [ Left  Right  Bottom Top ]
     % ignored if NonEmpty FigPos

     PaperBorder = [ 1  1.0  1.2   1.2 ];
  
     FigPos = [ 32 100 1024-2*32 560 ];   % Ignore PaperBorder

     % Normalized AxesPosition

     AxePos = [ 0.00  0.00  1.00  1.00 ];

     % Normalized Space of Axes in Figure

     Xoff = [ 0.05   0.12 ];   % [ Left Right  ]
     Yoff = [ 0.12   0.08 ];   % [ Top  Bottom ]

     AxeFontSize   = 14;
     LabelFontSize = 12;

     TickLabelIntervall = 12;   % Intervall for TIMEAXIS

    %**********************************************
   case 'hb800'
   %**********************************************

     % Borders at Paper [ Left  Right  Bottom Top ]
     % ignored if NonEmpty FigPos

     PaperBorder = [ 1  1.0  1.2   1.2 ];
  
     FigPos = [ 32 100 800 600 ];   % Ignore PaperBorder

     % Normalized AxesPosition

     AxePos = [ 0.00  0.02  1.00  0.95 ];

     % Normalized Space of Axes in Figure

     Xoff = [ 0.04   0.10 ];   % [ Left Right  ]
     Yoff = [ 0.10   0.10 ];   % [ Top  Bottom ]

     AxeFontSize   = 10;
     LabelFontSize = 10;

     TickLabelIntervall = 12;   % Intervall for TIMEAXIS

    %**********************************************
   case 'dm960'
   %**********************************************

     % Borders at Paper [ Left  Right  Bottom Top ]
     % ignored if NonEmpty FigPos

     PaperBorder = [ .8  0.8  0.8 0.8 ];
  
     FigPos = [ 32 100 960 720 ];   % Ignore PaperBorder

     % Normalized AxesPosition

     AxePos = [ 0.00  0.02  1.00  0.95 ];

     % Normalized Space of Axes in Figure

     Xoff = [ 0.03   0.08 ];   % [ Left Right  ]
     Yoff = [ 0.08   0.08 ];   % [ Top  Bottom ]

     AxeFontSize   = 10;
     LabelFontSize = 10;

     TickLabelIntervall = 12;   % Intervall for TIMEAXIS

  %**********************************************
   otherwise
   %**********************************************

     str = sprintf('%s, ',modes{:});
     str = str( 1 : end-2 );

     msg = cat(1,msg,{sprintf('Invalid Mode, use any of %s.',str)});
   
%**********************************************
end
%**********************************************

%**********************************************
% Check CellArrays for Time, U and V
 
if iscell(t) & iscell(u) & iscell(v)
   t = t(:);
   u = u(:);
   v = v(:);
   if ~isequal(size(t),size(u),size(v))
       msg = cat(1,{'CellArrays for Time, U and V must have the same Size.'},msg)
   else
       for ii = 1 : prod(size(t))
           ok = isequal(size(t{ii}),size(u{ii}),size(v{ii}));
           if ~ok
               break
           end
       end
       if ~ok
           msg = cat(1,{'Elements of CellArrays for Time, U and V must have the same Size.'},msg)
       end
   end
elseif ~( isnumeric(t) & isnumeric(u) & isnumeric(v) )
   msg = cat(1,{'Inputs for Time, U and V must be all CellArrays or all Numerics.'},msg);
else
   if ~isequal(size(u),size(v))
       msg = cat(1,{'Matrices for U and V must have the same Size.'},msg);
   else
       si = size(t);
       if ~( max(si) == prod(si) )
           msg = cat(1,{'Numeric Time must be a Vector.'},msg);
       elseif ~any( max(si) == size(u) )
           msg = cat(1,{'Length of Time must match Size of U and V.'});
       else
           if si(2) == max(si)
              t = t(:);
              if size(u,2) == si(2)
                 u = permute(u,[2 1]);
                 v = permute(v,[2 1]);
              end
           elseif ~( size(u,1) == si(1) )
              u = permute(u,[2 1]);
              v = permute(v,[2 1]);
           end              
           t = {t};
           t = t(ones(1,size(u,2)));
           u = num2cell(u,1);
           v = num2cell(v,1);
           t = t(:);
           u = u(:);
           v = v(:);
       end
   end       
end

%---------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Invalid Inputs.\n%s',msg));
end

%**********************************************

n = size(u,1);

vmax = zeros(n,1);
vmin = zeros(n,1);
t01  = zeros(n,2);
rat  = zeros(n,1);

for ii = 1 : n

  vmax(ii) = max(abs(v{ii}));
  vmin(ii) = min(min(v{ii}),0);

  rat(ii) = ( max(v{ii})*(max(v{ii})>0) - min(v{ii})*(min(v{ii})<0) );

  if all(isnan(u{ii}))
      rat(ii) =  2 * vmax(ii);
     vmin(ii) = -1 * vmax(ii);
  end

  t01(ii,1) = min(t{ii});
  t01(ii,2) = max(t{ii});

end

ratsum=max(cumsum(rat));

%--------------------------------------------------------------

if isempty(unit)
   vmax = max(vmax);
   pt10 = floor(log10(vmax));
   unit = ( floor(vmax/10^pt10) + 1 ) * 10^pt10;
end

if isempty(ampl)
   ampl = 1.5 * ratsum/(2*unit)
end

if isempty(labels)
   labels    = cell(n,1);
   labels(:) = { '' };
else
   labels = labels(:);
   s1 = size(labels,1);
   if s1 < n
      labels = labels(cat(1,(1:s1)',ones(n-s1,1)));
      labels( si+1 : n ) = {''};
   end
end

%--------------------------------------------------------------

ymax   = 2*unit*ampl;

yrange = rat./ratsum * ymax;

yminl  = cumsum(yrange);

y0     = yminl + (yrange./rat) .* vmin;

ymin   = 0;

%--------------------------------------------------------------
% AxesLimits

nl = 0;
for ll = labels(:)'
    nl = max( nl , size(char(ll{1}),2) );
end

Xoff(2) = ( Xoff(2) - Xoff(1) ) / 10 * nl + Xoff(1);

xlim0 = [ min(t01(:,1)) max(t01(:,2)) ];
xlim  = xlim0 + ([-1 1].*Xoff) * diff(xlim0)/(1-sum(Xoff));

ylim0 = [ -ymax  ymin ];

ylim  = ylim0 + ([-1 1].*Yoff) * diff(ylim0)/(1-sum(Yoff));


dx = diff(xlim);
dy = diff(ylim);

%--------------------------------------------------------------
% TimeAxis-Position

xt01 = Xoff / AxePos(3);

AxePosT    = AxePos;
AxePosT(1) = ( xlim0(1)-xlim(1) ) / dx * AxePos(3) + AxePos(1);
AxePosT(3) = diff(xlim0) / dx * AxePos(3);

AxePosT(4) = 2/3 * AxePosT(4);
AxePosT(2) = AxePosT(2) - AxePosT(4);

%*********************************************************
% Figure

fig = figure( 'paperunits'      , 'inches'   , ...
              'paperorientation', 'portrait' , ...
              'units'           , 'pixels'   , ...
              'position'        , [ 100 50 200 200 ] , ...
              'color'           , [ 1  1  1 ] , ...
              'nextplot'        , 'add'              );


if isempty(FigPos)

   pap_si = get(fig,'papersize');

   pappos = [ PaperBorder([1 3])  ...
           pap_si(1)-sum(PaperBorder([1 2]))  pap_si(2)-sum(PaperBorder([3 4])) ];

   set( fig , 'paperposition' , pappos );
  
   wysiwyg

else

   set( fig , 'position' , FigPos );

   wygiwys

   pappos = get(fig,'position');

end

drawnow

%*********************************************************
% TimeAxes

axeT = axes( 'parent' , fig , ...
            'units'  , 'normalized' , ...
            'position' , AxePosT , ...
            'fontsize' , AxeFontSize , ...
            'fontweight' , 'bold' , ... 
            'box'      , 'off'  , ...
            'tickdir'  , 'out'   , ...
            'ytick'    , []     , ...
            'ycolor'   , [ 1 1 1 ] , ...
            'xcolor'   , [ 0 0 0 ] , ...
            'linewidth' , 0.5 , ...
            'xaxislocation' , 'top' , ...
            'xlim'     , xlim0  , ...
            'color'    , 'none' , ...
            'clipping' , 'on'   , ...
            'nextplot' , 'add'          );
   
timeaxis(axeT,'x',TickLabelIntervall,base);

%*********************************************************
% StickAxes

axe = axes( 'parent' , fig , ...
            'units'  , 'normalized' , ...
            'position' , AxePos , ...
            'fontsize' , AxeFontSize , ...
            'box'      , 'on'   , ...
            'tickdir'  , 'in'   , ...
            'xlim'     , xlim   , ...
            'ylim'     , ylim   , ...
            'xtick'    ,  []    , ...
            'ytick'    ,  []    , ...
            'color'    , 'none' , ...
            'clipping' , 'on'   , ...
            'nextplot' , 'add'        );
  

%-------------------------------------------------------------

axepos = get(axe,'position');

k=(pappos(4)*axepos(4))/(pappos(3)*axepos(3))*dx/dy;
% die Laenge von einer y-Einheit im Druckbild ist die von k x-Einheiten 
% das ist der Faktor fuer die u-Laenge der Sticks in x-Richtung

%-------------------------------------------------------------
% Unit-Arrow

% po = ylim0(2) + 0.05*dy

po = -y0(1) + max(max(v{1}),0) + 0.03*dy;

pf15 = xlim0(1);


trieast =[pf15+k*unit        po
          pf15+k*2/3*unit   po-0.05*unit
          pf15+k*2/3*unit   po+0.05*unit ];

patch( 'parent' , axe , ...
       'xdata'  , trieast(:,1) , ...
       'ydata'  , trieast(:,2) , ...
       'facecolor' , 'k'       , ...
       'edgecolor' , 'k'       , ...
       'tag'       , 'ArrowHead'    ); 


line( 'parent' , axe  , ...
      'xdata'  , pf15+[0 k*unit] , ...
      'ydata'  , [ po po ] , ...
      'color'  , 'k' , ...
      'linestyle' , '-' , ...
      'marker'    , 'none' , ...
      'tag'       , 'ArrowLine'  )

text( 'parent' , axe , ...
      'position' , [ pf15+k*0.50*unit  po+0.04*unit  0 ] , ...
      'string'   , sprintf( '%.3g cm/s' , unit ), ...
      'Horizontalalignment' , 'center'  , ...
        'Verticalalignment' , 'bottom'  , ...
        'FontSize'          , AxeFontSize   , ...
        'FontWeight'        , 'bold' , ...
        'tag'               , 'ArrowLabel'          );

%---------------------------------------------
% New Position of Title

tt  = get(axe,'title');
uni = get(tt,'units');

set(tt,'units','data');

pos = [ xlim(1)+0.5*dx ylim(2)-0.02*dy  0 ];

set( tt , 'position' , pos , ...
          'horizontalalignment' , 'center' , ...
          'verticalalignment'   , 'top'    , ...
          'fontweight'          , 'bold'   , ...
          'fontsize'            , 12              );

set( tt , 'units' , uni );

%---------------------------------------------
% Plot Sticks

for ii = 1 : n 

 % si(2) Zeitserien

 jj=find( ~isnan(v{ii}) );

 anf  = min(jj);
 ende = max(jj);
 

 % Basislinie

 line( 'parent' , axe , ...
       'xdata'  ,  t{ii}([anf ende]) , ...
       'ydata'  , -y0([ii ii]) , ...
       'color'  , 'k' , ...
       'linestyle' , '-' , ...
       'marker'    , 'none' , ...
       'tag'       , sprintf('BaseLine%2.2d',ii) , ...
       'userdata'  , labels{ii} );


 x = t{ii} + k.*u{ii};
 y = v{ii} -   y0(ii);

 
 xx = NaN*ones(3,prod(size(x)));

 yy = NaN*ones(3,prod(size(x)));
 
 xx(1,:) = t{ii}(:)';
 xx(2,:) = x(:)';

 yy(1,:) = -y0(ii)+0*y(:)';
 yy(2,:) = y(:)';


 patch( 'parent' , axe , ...
        'xdata'  , xx(:) , ...
        'ydata'  , yy(:) , ...
        'edgecolor'  , 'k'   , ...
        'linestyle'  , '-'   , ...
        'clipping','off' , ...
        'tag'     ,sprintf('StickPatch%2.2d',ii) , ...
        'userdata' , labels{ii}        )

 % Labels
 
 text(  'parent'            , axe , ...
        'position'          , [ xlim(2)-0.02*dx , -y0(ii) 0 ] , ... 
        'String'            , labels{ii} , ...
        'Color'             , 'k'      , ...
      'Horizontalalignment' , 'right'  , ...
        'Verticalalignment' , 'middle' , ...
                 'Fontsize' , LabelFontSize , ...
                     'tag'  , sprintf('Label%2.2d',ii) );
 

end


%*******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,base,unit,ampl,labels,mode] = checkin(vin,n0);

base    = [];
unit    = [];
ampl    = [];
labels  = {};
mode    = '';

msg = cell(0,1);

for ii = 1 : prod(size(vin))

    v = vin{ii};

    if isnumeric(v)
       if isequal(size(v),[1 3])
          base = v;
       elseif ~( prod(size(v)) == 1 )
          msg = cat(1,msg,{'Numeric Input must be single.'});
       else 
          vr = real(v); vi = imag(v);
          if ( vi == 0 ) & ( 1000 <= vr ) & ( vr <= 3000 )
             base = [ vr 01 01 ];
          else
             unit = vr; 
             ampl = vi;
          end
       end
    elseif chkstr(v,0)
       mode = v;
    elseif chkcstr(v,0)
       [ok,labels] = chkcstr(v);
    else
       msg = cat(1,msg,{sprintf('Unrecognized %.0f. Input.',ii+n0)});
    end

end

if ~isempty(unit)
    if unit == 0
       unit = [];
    elseif ~isfinite(unit)
       msg = cat(1,msg,{'Unit must be finite.'});
    end
end

if ~isempty(ampl)
    if ampl == 0
       ampl = [];
    elseif ~isfinite(ampl)
       msg = cat(1,msg,{'Amplitude must be finite.'});
    end
end
