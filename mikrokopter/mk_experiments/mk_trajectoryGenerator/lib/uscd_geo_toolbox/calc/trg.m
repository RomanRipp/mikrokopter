function x = trg(x,varargin);

% TRG   Generates a V-shape Triangle-Function
%
% Y = TRG(X)  , Period = 2*pi
%
% TRG( 0.0*pi ) =  1
% TRG( 0.5*pi ) =  0
% TRG( 1.0*pi ) = -1
% TRG( 1.5*pi ) =  0
% TRG( 2.0*pi ) =  1
%
%--------------------------------------------------------
%
% Y = TRG( X , P + i*V )  Modified Triangle with Period P
% Y = TRG( X , i*V )        use P = 1
%
%   V: Location of Extrema, normalized to Period, ZERO at P/2
%
%   V > 0 moves the lower Extrema from P/2 to right (P)  \|
%   V < 0 moves the lower Extrema from P/2 to left  (0)  |/
%
% lower Extrema: TRG( N*P , P+i*V ) = -1;    N = V+0.5
%
%         Zeros: TRG(    N  * P/2 , P+i*V ) = 0
%                TRG( (1+N) * P/2 , P+i*V ) = 0
% 
% defaults: P = 2*pi; V = 0.0;  normal Triangle
%
%--------------------------------------------------------
%
% Y = TRG(  X , P + i*V , W )  Hold Value arround Extrema
%
%   W:  Half Width of Window arround Extrema, normalized to Period
%
% TRG( X , P + i*V , W ) =  1;  X <= P*W  |  P*(1-W) <= X
%
% TRG( X , P + i*V , W ) = -1;  P*(N-W) <= X <= P*(N+W)
%
% example: TRG(X,1,0.25) returns SIGN(COS(2*PI*X))
%
%   W = [ W1  W2 ] use W1 for   0 <= X <= P/2
%                      W2 for P/2 <= X <= P
% 
% TRG( X , P + i*V , [ W1  W2 ] ) =  1;  X <= P*W1  |  P*(1-W2) <= X
%
% TRG( X , P + i*V , [ W1  W2 ] ) = -1;  P*(N-W1) <= X <= P*(N+W2)
%
%--------------------------------------------------------
%
% Y = TRG(  X , ... , Mode , ... ) use the Modes:
%
%  'l'   Linear  (default)
%  'c'   Cosine  example: TRG(X,'c')   returns COS(X)
%                example: TRG(X,'c',1) returns COS(2*PI*X)
%
%  'cl'  Cosine: 0   <= X <= N*P
%        Linear: N*P <  X <  P
%
%  'lc'  Linear: 0   <= X <= N*P
%        Cosine: N*P <  X <  P
%
%
%--------------------------------------------------------
% Example:
%
%  x = ( -1 : 0.001 : 3 );
%
%  figure, hold on, box on, grid on
%  
%  plot(x,trg(x,1),'k-');
%  plot(x,trg(x,2.0),'k--');
%
%  plot(x,trg(x,1-0.3*i),'r-');
%  plot(x,trg(x,1+0.2*i),'g-');
%
%  figure, hold on, box on, grid on
%
%  x = linspace( -pi , 6*pi , 10000 );
%  plot(x,cos(x/2),'k-');
%
%  plot(x,trg(x,4*pi),'b-');
%  plot(x,trg(x,4*pi,0.1),'g-');
%  plot(x,trg(x,4*pi-0.3*i),'r-');
%  plot(x,trg(x,4*pi-0.3*i,[0.05 0.1]),'m-');
%  
%  figure, hold on, box on, grid on
%
%  x = ( 0 : 0.01 : 10 );
%
%  plot(x,trg(x,'c' ,5-0.1*i),'b-');
%  plot(x,trg(x,'cl',5-0.1*i,[0 0.1]),'r--');
%  plot(x,trg(x,'cl',5+0.5*i,[0.2 0]),'g-');
%  

%****************************************

Nin = nargin;

if Nin == 0
   x = [];
elseif ~isnumeric(x)
   error('X must be numeric.');
elseif ~strcmp(class(x),'double')
   x = double(x);
end

if isempty(x)
   return
end

v = [];
w = [];
m = 'l';

mm = 'lc';

for ii = 1 : Nin-1
    y = varargin{ii};
    s = size(y);
    p = prod(s);
    if ~( p == 0 )
        if ischar(y)
           y = lower(y(1:min(2,p)));
           y = cat(2,y(:)',' ');
           if any( y(1) == mm )
              i2 = ( any( y(2) == mm )  &  ~( y(1) == y(2) ) );
              m  = y( 1 : 1+i2 );
           end
        elseif isnumeric(y) 
           if ( p == 1 ) & isempty(v)
               v = y;
           elseif isempty(w)
               w = y;
           end
        end
    end
end

%----------------------------------------
% Defaults

if isempty(v)
   v = 2*pi;
end

if isempty(w)
   w = [ 0  0 ];
elseif ( prod(size(w)) == 1 )
   w = [ w  w ];
else
   w = w(1:2);
end

%****************************************

p = real(v);    % Period
v = imag(v);    % Deviation of lower Extrema

v = v + 0.5;
p = p + 2*pi * ( p == 0 );

if ~( p == 1 )
      x = x / p;
end

x = x - floor(x);    % [ 0 .. 1 ) Periodic

if ( v == 0.5 ) & all( w == 0 ) & ( size(m,2) == 1 )
   switch m
     case 'l'
        x = -2 * sign(x-0.5) .* ( 1 - 2*x ) - 1; 
     case 'c'
        x = cos( 2*pi * x );
   end  
   return
end


v = v - 1e2*eps * any( v == [ 0  1 ] );

ww = warnstat;
     warning('off');

if size(m,2) == 2

   i1 = ( x <= v );
   i2 = ( x >  v );

   x =  i1 .* min(max( ( 1/(v-0-2*w(1)) * (x-(0+w(1)))     ) , 0 ) , 1 )  +  ...
        i2 .* min(max( ( 1/(1-v-2*w(2)) * (x-(v+w(2))) + 1 ) , 1 ) , 2 );

else

   x =  ( x <= v ) .* min(max( ( 1/(v-0-2*w(1)) * (x-(0+w(1)))     ) , 0 ) , 1 )  +  ...
        ( x >  v ) .* min(max( ( 1/(1-v-2*w(2)) * (x-(v+w(2))) + 1 ) , 1 ) , 2 );

end

warning(ww);

switch m

   case 'c'

      x = cos( 2*pi * min(max(0,x/2),1) );

   case 'l'

      x = 2 * abs( 1 - x ) - 1; 

   case 'cl'

      x = i1 .* cos( 2*pi * min(max(0,x/2),1) ) + ...
          i2 .* ( 2 * abs( 1 - x ) - 1 ); 

   case 'lc'

      x = i2 .* cos( 2*pi * min(max(0,x/2),1) ) + ...
          i1 .* ( 2 * abs( 1 - x ) - 1 ); 

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



