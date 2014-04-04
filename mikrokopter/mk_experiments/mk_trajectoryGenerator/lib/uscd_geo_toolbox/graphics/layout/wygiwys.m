function wygiwys(fig,ornt)

% WYGIWYS  WhatYouGetIsWhatYouSee
%
% Switch the PaperPosition to the same View like on the Screen
% The PaperOrientation will adapted to the Extension of the Figure.
% The Figure will positioned in the Center of the Paper
%
% WYGIWYS( FigureHandle )
%
% WYGIWYS( FigureHandle , Orientation )
%
%   use the defined PaperOrientation: 'portrait' | 'landscape'
%


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

mode = { 'portrait'  'landscape' };
if nargin == 2
   if ~( ischar(ornt) & ~isempty(ornt) & ...
         ( prod(size(ornt)) == size(ornt,2) ) )
       error('Orientation must be a String.');
   end
   ornt = lower(ornt(1));
   mode = mode{ 1 + strcmp(ornt,'l') };
   set( fig , 'paperorientation' , mode );
end
   
figuni = get(fig,'units');
papuni = get(fig,'paperunits');

set(fig,     'units' , 'pixels' , ...
        'paperunits' , 'inches'       );

figpos = get(fig,'position');

ppi    = get(0,'screenpixelsperinch');

pappos = zeros(1,4);

pappos([3 4]) = figpos([3 4]) / ppi;

if ~ischar(mode)

    set( fig , 'paperorientation' , mode{1} );

   pap_si = get(fig,'papersize');

   sc     = pap_si ./ pappos([3 4]);
   flip   = ( sc(1) < sc(2) );
   flip   = ( flip & ( ( sc(1) < 1 ) | ( pappos(4) < pappos(3) ) ) );
   mode   = mode{ 1 + flip };

   pap_si = pap_si([1 2]+[1 -1]*flip);

   set( fig , 'paperorientation' , mode );

else

   pap_si = get(fig,'papersize');

end

pappos([1 2]) = (pap_si-pappos([3 4])) / 2;

set(fig,'paperposition',pappos);

set(fig,     'units' , figuni   , ...
        'paperunits' , papuni         );
