function wysiwyg(fig)

% WYSIWYG  WhatYouSeeIsWhatYouGet
%
% WYGIWYS( FigureHandle )
%
%WYSIWYG -- changes the size of the figure on the screen to equal
%       the size of the figure that would be printed, 
%       according to the papersize attribute.  Use this function
%       to give a more accurate picture of what will be 
%       printed.
%       Dan(K) Braithwaite, Dept. of Hydrology U.of.A  11/93
 
if nargin < 1
  fig = get(0,'currentfigure');
end

if isempty(fig)
 return
end

ok = ( isnumeric(fig)  &  ( prod(size(fig)) == 1 ) );
if ok
   ok = ishandle(fig);
   if ok
      ok = strcmp( get(fig,'type') , 'figure' );
   end
end

if ~ok
   error('Input must be a FigureHandle.');
end

unis = get(fig,'units');
ppos = get(fig,'paperposition');

org = get(fig,'position');

set(fig,'units',get(fig,'paperunits'));

pos = get(fig,'position');
pos(3:4) = ppos(3:4);

set(fig,'position',pos);
set(fig,'units',unis);
  
new = get(fig,'position');

if ~all( new([1 2]) == org([1 2]) )
    new([1 2]) = org([1 2]);
    set(fig,'position',new);
end
