function x = corpus(m,x);

% CORPUS  Returns a unit (rotational) Corpus
%
%   Z = CORPUS( M+i*F , X ); with N = LENGTH(X)
%   Z = CORPUS( M+i*F , N ); with X = LINSPACE( -1 , 1 , N )
%
% returns the Z-Values of an [ N by N ] - Corpus 
%
% M gives the Potenz to modify the shape of the Inclination,
%    use NaN for a Half-Sphere
%
% F Number of Flat Faces
%
%   The fractional Part of F gives the Rotation
%    in Range:  [ 0 .. 360/fix(F) )
%
%   CORPUS(0)        Half UNIT-Sphere
%   CORPUS(1)        Cone
%   CORPUS(1+4.5*i)  Pyramid
%
% X Vector for Corpus, N Number of Rows/Columns
%
% defaults: M = 1    (Cone)
%           F = 0    (rotational)
%           N = 100
%
%
%------------------------------------------------------------
% Examples to demonstrate the Influence of the Potenz
%
%------------------------------------------------------------
% Animated Image
%------------------------------------------------------------
%
%  n = 201;
%  l = [ -1  1 ];
%  c = l + diff(l)/(n-1)/2 * [ 1 -1 ];
%
%  p = ( 0.5 : 0.05 : 2 );
%
%
%  fig = figure( 'position' , [300 200 n n] , 'colormap' , jet(128) );
%
%  axe = axes( 'position' , [0 0 1 1] , 'visible' , 'off' , ...
%      'xlim' , l , 'ylim' , l , 'clim' , [0 1] , 'nextplot' , 'add' );
%
%  img = image( 'xdata' , c , 'ydata' , c' ,'cdata' , []  , ...
%               'cdatamapping' , 'scaled' , 'erasemode' , 'none' );
%
%  while 1
%      for ii = cat( 2 , p , p(end:-1:1) );
%          set(img,'cdata',corpus(ii,n));  drawnow
%      end
%  end
%
%------------------------------------------------------------
% Animated Surface
%------------------------------------------------------------
%
%  l = [ -1  1 ]; x = ( l(1) : 0.005 : l(2) );
%
%  p = ( 0.1 : 0.1 : 3 ); n = size(p,2);
%
%  msk = sqrt( x(ones(size(x)),:).^2 + (x(ones(size(x)),:)').^2 );
%  msk = find( msk > 1 );
%
%  z = corpus(p(1),x); z(msk) = NaN;
%
%  fig = figure( 'position' , [300 300 300 250] , ...
%                'color' , 'w' , 'colormap' , copper(128) );
%
%  axe = axes( 'position' , [0 0.02 1 1] , 'visible' , 'off' , ...
%      'xlim' , 0.8*l , 'ylim' , 0.8*l , 'zlim' , [0 0.8] , 'clim' , l , ... 
%      'dataaspectratio' , [1 1 1] , 'view'  , [20 30] , 'nextplot' , 'add' );
% 
%  srf = surface( 'xdata' , x , 'ydata' , x' , ...
%                 'zdata' , z , 'cdata' , z  , ...
%             'facecolor' , 'interp' , 'edgecolor' , 'none' , ...
%          'cdatamapping' , 'scaled' , 'clipping' , 'off'      );
%
%  material([0.4 0.6 0.0 10 1]);
%
%  lgt = light('color',[1 1 1],'style','infinite'); lightangle(lgt,90,60);
%
%  M = getframe(fig); M = M(1,ones(1,n));
%
%  for ii = 2 : n
%      z = corpus(p(ii),x); z(msk) = NaN;
%      set(srf,'zdata',z,'cdata',z);
%      M(ii) = getframe(fig);
%  end
%
%  delete([lgt srf])
%
%  movie(axe,M,-100,24);
%


Nin = nargin;

if Nin < 1
   m = 1;
end

if Nin < 2
   x = 100;
end

if isempty(x)
   return
end

if prod(size(x)) == 1

   if x <= 0
      x = [];
      return
   end

   n = max(floor(x),1);

   % LINSPACE( -1 , 1 , N )
   x = cat( 2 , -1+(0:n-2)*2/(n-1) , 1 );

else

   s = size(x);
   n = max(s);

   if ~( prod(s) == n )
       error('N must be a Scalar or Vector.');
   end

   if ~( s(2) == n ) 
       x = n(:)';
   end

end

x = x(ones(1,n),:);
   
k = imag(m);
m = real(m);

l = k - fix(k);
k = round(k-l);

if k == 0

   x = sqrt( x.^2 + (x').^2 );

else

   p = 2*pi;

   k = p / k;
   l = l * k;
 
   w = atan2(x',x) - l;
   w = w - p * floor(w/p);  % [ 0 .. p )

   w = k * floor( w / k ) + k/2;

   if  l == 0
       x = x .* cos(w) + (x') .* sin(w);
   else
       x = (  x*cos(l) + (x')*sin(l) ) .* cos(w) + ...
           ( -x*sin(l) + (x')*cos(l) ) .* sin(w);
   end

end

if m == 0

   x = sqrt( 1 - min(x,1).^2 );

else

   x = 1 - x.^m;

end

