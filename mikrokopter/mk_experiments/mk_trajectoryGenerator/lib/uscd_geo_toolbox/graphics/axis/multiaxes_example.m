% MULTIAXES    creates multiple Axes in One Direction
%
% AxeHandle = MULTIAXES( N , XY , Lim , AxePos , AxeMin )
%
% Example:
%
%% CTD = [ P   T   S ]
%  CTD = [ 0  16   7
%          2  15   8
%          4  15   8
%          6  14   9
%          8  13  10
%         10  11  12     
%         12   9  15
%         14   8  16
%         16   8  16  ];
%
% N = 2;  % 2 Variables:  Temp and Sal
%
% xlims  = [  5  20 ; ...
%             5  20       ];
%
% colors = [ 1 0 0 ; ...     % Temp Red
%            0 0 1       ];  % Sal  Blue
%
% labels = [ ['Temparature   [' setstr(176) 'C] ' ]
%            ['Salinity      [PSU]'               ]    ];
%
% figure('paperunits'      , 'centimeter'       , ...
%        'paperorientation', 'portrait'         , ...
%        'paperposition'   , [ 2  4  16  20 ]   , ...
%        'position'        , [ 100 100 200 200 ]      );
% wysiwyg
% drawnow
%
% axe = multiaxes(2,'x',[18 0]);
%
% for ii = ( N : -1 :  1 )
%   axes(axe(ii))
%   set(axe(ii),'xcolor', colors(ii,:), ...
%               'xlim'  , xlims(ii,:) , ...
%               'box'   , 'off'              );
%   plot(CTD(:,ii+1),CTD(:,1),'color',colors(ii,:));
%   xlabel(labels(ii,:))
% end
%
% axes(axe(N+1))
%  set(axe(N+1),'xgrid','on', ...
%               'ygrid','on', ...
%               'xlim' , get(axe(1),'xlim') , ...
%               'xtick', get(axe(1),'xtick')      )
% 
% axes(axe(1))
% ylabel('Pressure [dbar]')
%


