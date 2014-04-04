function [axe1,axepos1,lim1,lim_fak] = multiaxes(N,xy,lim,axepos,axe_min)

% MULTIAXES    creates multiple Axes in One Direction
%
% MULTIAXES( N , XY , Lim , AxePos , AxeMin )
%
%       N       Number of the Axes to Create
%
%      XY       Specifies the Dimension of multiplying:  'x'  or 'y'
%                default: 'x'
%
%      Lim      AxesLimits of the Common Axis (default: Y-Axes)
%                 [ Min  Max ]  for Direction 'NORMAL'
%                 [ Max  Min ]  for Direction 'REVERSE'  !!!
%                default:  [ 0   1 ]
% 
%      AxePos   AxesPosition of the Inner Axes
%                [ 1 by 4 ] Vector: [ AxeLeft AxeBott AxeWidth AxeHigh ] 
%                                     Normalized to Figure !!!!
%
%      AxeMin   Minimum Distance of the Outher Axis to Border
%                default: 0.06  (normalized to Figure)
%
%
%  AxeHandles = MULTIAXES( ... )  returns the Handles, [ N+1 by 1 ]
%
%     AxeHandles(1) is the Inner Axe
%
%     The last, additional AxeHandle represents the Box  (Top and Right)
%         surrounding the Inner Axe, so are the AxeLine and the YTicks 
%         represented by   the first and last AxeHandle, i.e. to change
%          Color and Tick - Properties of the AxeLine, you have to set
%          both Axes !
%
%
% For following plots, start with the Outher Axes, like
%
%
% for ii = N : -1 : 1
%
%   axes(AxeHandles(ii))
%
%      set( gca , ... )
%     plot( ...       )
%   xlabel( ...       )
%
% end
%
%
%
% Other OutPuts:
%
% [ AxeHandles , AxePos , Lims , LimFak ] = MULTIAXES( ... )
%
%    returns 
%       the Position of the Axes     in        AxePos ( [ N+1 by 4 ] ),
%       the Limits of the Other Axis in        Lims   ( [ N+1 by 2 ] ),
%       the Factor to calculate the Limits in  LimFak ( [ N+1 by 2 ] )
%           ( Lims =  Lim(1) - LimFak*diff(Lim); )
%
%
% Example:   type   >> help multiaxes_example 
%
% 







XY = 'xy';
YX = 'yx';



% Check InputArguments

% N,xy,lim,axepos,axe_min


Nin = nargin;

if Nin == 0
 error('Not enough InputArguments.')
end

if length(N) ~= 1 
 error(' Number of Axes must be  a single, finite, nonzero  IntegerValue.')
end

N = round(N);
if ~finite(N)  |  N <= 0
 error(' Number of Axes must be  a single finite, nonzero  IntegerValue.')
end



if Nin < 2
 xy = 'x';
end 

if ~isstr(xy)  |  length(xy) > 1
 error(' XY must be ''x''  or ''y'' .' )
end

ini = find( lower(xy) == 'xy' );  %  1 for 'x',  2 for 'y' !!!!!!

if isempty(ini)
 error(' XY must be ''x''  or ''y'' .')
end

% PositionIndex for AxePosition = [ AxeLeft AxeBott AxeWidth AxeHigh ]

pini = 3-ini;     %  2 for 'x' ,  1 for 'y'  !!!!!!!!!!!!!!!!!!!!!!!!

 

if Nin < 3   
 lim = [ 0  1 ];
end
if size(lim,1) ~= 1   |   size(lim,2)~=2
 error(' Lim must be a 2-Element Vector:  [ Left Right ]  or  [ Bott  Top ] ')
end

axe_default = 0; 
if Nin < 4
 axepos = get(0,'defaultaxesposition');
 axe_default = 1;
end

if size(axepos,1) ~= 1   |   size(axepos,2)~=4
 error(' AxePos must be a 4-Element Vector:  [ AxeLeft AxeBott AxeWidth AxeHigh ] ')
end


if Nin < 5
 axe_min = 0.06;
end
if length(axe_min) ~= 1
 error(' AxeMin must be a Value between 0 and 1 ')
end

if axe_min < 0  |   1 <  axe_min
 error(' AxeMin must be a Value between 0 and 1 ')
end

if axe_default

 axe_end = sum( axepos([pini pini+2]) );
 
 axe_min + (N-1)*0.06;

 axepos( pini   ) = axe_min + (N-1)*0.06;
 axepos( pini+2 ) = axe_end - axepos(pini); 

end

fig = gcf;



% Ceck for Color of Other Axe

  yx_color = get(fig,'color');
  if strcmp( yx_color , 'none' )
    yx_color = get(0,'defaultaxescolor');
    if strcmp( yx_color , 'none' )
      yx_color = [1 1 1];
      set(fig,'color',yx_color)
    end
  end
   

% Check for Direction of Other Axes
 
    yx_dir = 'normal';
  if diff(lim) < 0
    yx_dir = 'reverse';
  end


XY = XY(ini);
YX = YX(ini);




   axe1 = ones(N+1,1) ;               % Handles   of Axes
axepos1 = ones(N+1,1) * axepos;       % Position  of Axes
   lim1 = ones(N+1,1) * lim;          % YX-Limits of Axes
lim_fak = ones(N+1,1);



% Intervall between Axes, normalized to figure !

int = ( axepos(pini) - axe_min ) / ( N - 1 + [ N == 1 ] ) * ...
      ( 1 -  [ N == 1 ] );



%******************************************************************
% die ausserste Achse wird zuerst gezeichnet

Ind = ( N : -1 : 1 )';   % 1 Column !!!


 %******************************************************************
 % Calculate Positions relativ to Inner Axes
 %----------------------------------------------5--------------------

  axepos1(Ind,pini  ) = axepos(pini  ) - (Ind-1)*int;  % AchsenAbstand vom FigureRand
  axepos1(Ind,pini+2) = axepos(pini+2) + (Ind-1)*int;  % AchsenDehnung

  axepos1(N+1,:)  =  axepos1(1,:);  % Last DummyAxes


 %******************************************************************
 % Calculate Other Limit 
 %----------------------------------------------5--------------------

  lim_fak = axepos1(:,pini+2) / axepos(pini+2) - 1 ; % aktuelle / kleinste - AchsenDehnung

  lim1(:,1)  = lim(1) - lim_fak*diff(lim);


for ii = Ind(:)'

 %********************************************************************
 % Create new Axes and do settings
 %--------------------------------------------------------------------


  axe1(ii) = axes('units'    , 'normalized'       , ...
                  'position' , axepos1(ii,:)      , ...
                  'color'    , 'none'             , ...
             [YX 'lim']      , sort(lim1(ii,:))   , ... 
             [YX 'dir']      , yx_dir             , ...
             [YX 'color']    , yx_color           , ...
             [YX 'tick']     , []                 , ...
             [YX 'ticklabel'] , []               , ...
                'nextplot'    , 'add'             , ...
                'box'         , 'off'                      );

end




%***********************************************************************
% Die oberste (innerste) Achse wird sichtbar gemacht
%----------------------------------------------------------------------



% Achse 1 ist die innerste, y-sichtbare Achse !!!
 
set( axe1(1) ,  ...
       [YX 'color']         , get(0,['defaultaxes' YX 'color' ])  , ...
       [YX 'tickmode']      , 'auto'                   , ...
       [YX 'ticklabelmode'] , 'auto'                       )


% drawnow




%***********************************************************************
% Setze eine dummyAchse fuer den oberen, rechten Rand
%----------------------------------------------------------------------

 axe1(N+1) = axes( ...
             'units'           , get(axe1(1),'units')       , ...
             'position'        , get(axe1(1),'position')    , ...
             'box'             , 'on'                       , ...
             'color'           , 'none'                     , ...
             [YX 'lim']        , get(axe1(1),[YX 'lim'])    , ... 
             [YX 'dir']        , yx_dir                     , ...
             [YX 'color']      , get(axe1(1),[YX 'color'])  , ...
             [YX 'tick']       , get(axe1(1),[YX 'tick'])   , ...
             [YX 'ticklabel']  , []                         , ...
             [XY 'color']      , get(axe1(1),[XY 'color'])  , ...
             [XY 'tick']       , []                         , ...
             [XY 'ticklabel']  , []                              );
             




%***********************************************************************
% Make First Axes ACTIVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%----------------------------------------------------------------------

 axes(axe1(1))   % set(fig,'currentaxes',axe1(1))     

