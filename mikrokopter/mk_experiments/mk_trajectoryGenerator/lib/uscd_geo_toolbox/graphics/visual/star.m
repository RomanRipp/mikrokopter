function [x,r,d] = star(n,m)

% STAR  generate a 1. Order Star
%
% [ XY , R , D ] = START( N , P )
%
% N  = Number of Edges
%
% P  = Rotation [deg]  of Star
%      
% One outher Corner of the Star shows 
%  by default to Top, use a NonZero 
%  imaginary part of P to start 
%  with a Corner at right.
%  
% XY = Coordinates of Points
%       [ 2*N+1 by 2 ]
%
% R  = inner Radius
%
% D  = Rotation [deg] of inner Polyeder
%       relative to outher
%
% STAR without an OutputArgument 
%  opens a Figure which shows the result.
%
% see also: KOCHISLE
%
%---------------------------------------------
% Examples for Flowers:
%
% f =  2; p = -1; % Factor and Potenz for IterationRadius
% c = [7.3e-04 0.31 -1.9 6.0]; % Max. Number of Iteration
%
% for n = 5 : 13;
%
%     nc = ceil(polyval(c,n)); cl = linspace(0,1,nc);
%     cl = interp1(linspace(0,1,512),hot(512),cl.^0.5);
%     cl = min(max(cl,0),1);
%
%     rd = (f*n)^p;  % Min. Radius to iterate
%
%     figure 
%     axis(1.1*[ -1  1  -1  1 ]), axis manual, axis equal
%     hold on, box on, set(gca,'tickdir','out');
%
%     rad = 1; dev = 0; ic = 0;
%
%     while rad > rd
%           [xy,r,d] = star(n,dev);
%           ic = ic + 1 - ic * ( ic == size(cl,1) );
%           plot(rad*xy(:,1),rad*xy(:,2),'color',cl(ic,:));
%         rad = rad * r;, dev = dev + d;
%     end
% end
%

Nin = nargin;

if Nin < 1
   n = 5;
end

if Nin < 2
   m = 0;
end

if isempty(n)
   x = []; r = []; d = [];
   return
end


n = 2*n;

d = 2*pi / n;

p = ( 0 : n )';

%%% sin(pi/2-2*d)/cos(d)
%%% cos(2*d)/cos(d)

r = cos(d);
r = [ 1 ; 2*r - 1/r ];
r = max(r,-1);

r = r(mod(p,2)+1);

m = real(m) + 90 * ( imag(m) == 0 );

p = p * d + m*pi/180 ;

x = r(:,[1 1]) .* [ cos(p)  sin(p) ];

r = r(2);
d = d * 180/pi;

if nargout == 0

   figure, 
   axes( 'xlim' , [ -1.1  1.1 ] , ...
         'ylim' , [ -1.1  1.1 ] , ...
         'box'  , 'on' , ...
         'tickdir' , 'out' , ...
         'dataaspectratio' , [ 1 1 1 ] , ...
         'nextplot' , 'add' );

   plot(x(:,1),x(:,2),'k-');

   ii = [ ( 1 : 2 : n )  1 ];

   plot(x(ii,1),x(ii,2),'k--');

   ii = ii + 1;
   ii = ii - n * ( ii > n );

   plot(x(ii,1),x(ii,2),'k:');

   p = ( 0 : 360 ) * pi / 180;

   plot(sin(p),cos(p),'k:')

   clear x

end
