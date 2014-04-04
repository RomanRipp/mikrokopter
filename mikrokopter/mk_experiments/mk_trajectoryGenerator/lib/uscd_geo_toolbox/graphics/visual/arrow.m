function  [xx,yy]= arrow(xx0,yy0,tip,head,bas,width);

% ARROW    Calculate an Arrow for Lines
%
%   [ XX , YY ] = ARROW( XX0 , YY0 );
% 
%   Gives the 3 Coordinates of the Head of the Arrow
%
%    [XX0(1,:),YY0(1,:)] ----> [XX0(2,:),YY0(2,:)]
%
%    ( One Arrow, Head per Column:
%        XX0, YY0 ...  [ 2 by N ] ([ Start ; End ])
%        XX , YY  ...  [ 7 by N ] ([ Tail1 ; Pike(End) ; Tail2 ]) )
%
%    Use:
%
% >> plot( [ XX0 ; NaN*ones(1,size(XX0,2)) ; XX ] , ...
%          [ YY0 ; NaN*ones(1,size(XX0,2)) ; YY ]       )
%
%    to plot the Arrows.
%
%
%   [ XX , YY ] = ARROW( ... , TipAngle , HeadLength , BaseAngle , Width );
%
%   defines the TipAngle (degree) and 
%           the Length of the ArrowHead (normalized to the ArrowLength).
%
%   defaults:  TipAngle = 15;  HeadLength = 3/12;
%


Nin = nargin;

if Nin < 2
   error('Not enough InputArguments.');
end

if ~isequal(size(xx0),size(yy0))
   error('XX0 and YY0 must have the same Size.')
end

if Nin < 3
   % TipAngle
   tip = 15;
end

if Nin < 4
   % Normalized Length of Head
   head = 3/12;
end

if Nin < 5
   % BaseAngle
   bas = 2*tip;
end

if Nin < 6
   % Width
   width = NaN;
end

tip   =  tip(1)*pi/180;
bas   =  bas(1)*pi/180;
head  = head(1);


if size(xx0,1) ~= 2  &  size(xx0,2) == 2
   xx0 = xx0';
   yy0 = yy0';
end




  % Normalized ArrowHead in angle 0 to X-Axis
   
   yy = [     head*tan(tip)  ; 0 ;    -head*tan(tip)  ];
   xx = [    ( 1 - head )    ; 1 ;  ( 1 - head )      ];

  if isnan(width)
     width = 1/3 * yy(1);
  end

  yb = width;
  xb = xx(1) + (yy(1)-yb)/tan(bas);

  yy = cat( 1 , yb , yb , yy , -yb , -yb );
  xx = cat( 1 , 0  , xb , xx ,  xb ,  0  );
 
  % InputData
  dx  = diff(xx0,1,1) ;
  dy  = diff(yy0,1,1) ; 
 
  len = sqrt( dx.^2 + dy.^2 );
  ang = atan2(dy,dx);
  
  
  % Scale Normalized Head
  xx = xx*len;
  yy = yy*len;

  N = size(xx,2);
  
  for jj = 1:N

    T = [ [  cos(ang(jj)) sin(ang(jj)) ] ; ...
          [ -sin(ang(jj)) cos(ang(jj)) ]        ];
   xy = [xx(:,jj) yy(:,jj)] * T;
   
   xx(:,jj) = xy(:,1) + xx0(1,jj); 
   yy(:,jj) = xy(:,2) + yy0(1,jj);

  end 

