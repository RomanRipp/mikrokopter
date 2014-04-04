function [ll,yd,xx] = way(wd,wb,x,y,a,b);

% WAY  Calculates way of a Trajectory trough Nozzle.
%
% [ Way , YD , XX ] = WAY( WD , WB , X , Y , A , B )
%
%--------------------------------------------------------
% Input:
%
% The unrotated Nozzle has the Extensions:
%
%   WD in Y-Direction
%   WB in X-Direction  (ZERO if RundStrahl)
%
% The Rotation of the Nozzle is defined by
%  the Angle B, measured anticlockwise positive.
%  (valid for ~( WB == 0 ) only)
%
% The Location of the Trajectory in the Nozzle is defined by 
%   the TargetPoint [ X , Y ] measured relative to NozzleCenter.
%
% The Direction of the Trajectory is defined by
%  the Angle A, measured anticlockwise positive.
%
% The Inputs can be Scalars or Matrices of same Size.
%
%--------------------------------------------------------
% Output:
% 
%   Way  Length of Way of Trajectory trough Nozzle
%
%   YD   Max. possible Length of Way of a Trajectory
%          with Direction A trough a Nozzle with Rotation B
%
%   XX   Distance of Trajectory to NozzleCenter
% 

Nin  = nargin;
Nout = nargout;

%************************************************************
% Check Inputs

if Nin < 3
   x = 0;
end
if Nin < 4
   y = 0;
end
if Nin < 5
   a = 90;
end
if Nin < 6
   b = 0;
end


%------------------------------------------------------------
% Compare Size of Inputs

si = { size(wd) size(wb) size(x) size(y) size(a) size(b) };
for ii = 2 : size(si,2)
    s0 = si{ii-1}; p0 = prod(s0);
    s1 = si{ii-0}; p1 = prod(s1);
    if ~( isequal(s0,s1) | any( [p0 p1] == 1 ) )
        error('Inputs must be scalars of Matrices of same Size.');
    end
    if ( p1 == 1 ) & ~( p0 == 1 )
       si{ii} = s0;
    end
end

si = si{end}; one = ( prod(si) == 1 );

%*************************************************************
% Project: TrajectoryDirection  ==  Y-Axes

%-------------------------------------------------------------
% Angle: Nozzle-Achses WD  / Trajectory
  
p = 90 - a + b;

%-------------------------------------------------------------
% X-Distance of [dx,dy] to Trajectory

a = a * pi/180;

x = y.*cos(a) - x.*sin(a);

if Nout > 2 
   xx = x;
   if ~one & ( prod(size(xx)) == 1 )
       xx = xx(ones(si));
   end
end

x = abs(x);
  
%***************************************
% Check for RundStrahl

jj = ( wb(:) == 0 );

ok = any(jj) + all(jj);

jj = find(jj);

if ~( ok == 0 )

    r = wd / 2;    % Radius

    d = 2 * sqrt( r.^2 - x.^2 );
    d = d .* ( x < r );

end

%***************************************
if ok == 2
%***************************************

   yd = 2 * r;
   ll = d;
   if ~one & ( prod(size(ll)) == 1 )
       ll = ll(ones(si));
   end

   if Nout > 1
      yd = r;
      if ~one & ( prod(size(yd)) == 1 )
          yd = yd(ones(si));
      end
   end

   return

%***************************************
end
%***************************************

%-------------------------------------------------------------
% Angle: Nozzle-Achses WD  / Trajectory

p = p - 360 * floor(p/360);                 % [  000 .. 360 )
p = p - 360 * ( p >= 180 );                 % [ -180 .. 180 )
p = p - 180 * sign(p) .* ( abs(p) > 90 );   % [ -090 .. 090 ]

p = abs(p) * pi / 180;

% Maximum vertical Extension at Center (X==0)

cp = cos(p) / 2;
sp = sin(p) / 2;

x0 = -wb .* cp + wd .* sp;    % Max. Extension of constant Y
x1 =  wb .* cp + wd .* sp;    % Max. Extension in X

ww = warnstat;
     warning('off');

yd = min( wd./cp , wb./sp );  % Max. Extension at Center (X==0)

     warning(ww);

yd = yd / 2;                  % Scale back cp,sp

% [x0 x1 yd]

x0 = ( x1 - abs(x0) );
x1 = ( x1 - abs(x ) );

x1 = x1 .* ( x1 > 0 );                               % ZERO if outside

x1 = x1 + ( 1 - x1 ) .* ( x0 == 0 ) .* ( x1 ~= 0 );  % ONE  if X1 == X0 
x0 = x0 + ( x0 == 0 );

ll = ( yd ./ x0 ) .* x1;

ll = max(min(ll,yd),0);

if Nout > 1
   if ~one & ( prod(size(yd)) == 1 )
       yd = yd(ones(si));
   end
end

%***************************************
% Check for RundStrahl

if ok == 0
   return
end

if ( prod(size(d)) == 1 ) 
   ll(jj) = d;    % Distance
else
   ll(jj) = d(jj);
end

if Nout > 1
   if ( prod(size(r)) == 1 ) 
      yd(jj) = 2 * r;    % 2 * Radius
   else
      yd(jj) = 2 * r(jj);
   end
end


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



