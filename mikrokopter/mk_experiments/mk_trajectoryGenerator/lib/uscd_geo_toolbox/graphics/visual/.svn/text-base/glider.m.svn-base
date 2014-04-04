function [v,f,d,m,c] = glider(n,w)

% GLIDER   Returns Faces and Vertices of a Glider
%
% [V,F,D,M,C] = GLIDER(GEO,NRM)
%
%--------------------------------------------------------
%
%  GEO    Geometry: [ Resolution  Length  Radius ]
%
%            Resolution for Cylindrical Corpus
%            Length of Glider
%            Radius of GliderCorpus
%
%         default: [ 36  100  5  ];
%
%  NRM    XYZ-Vector(s) for Direction and Tilt of Glider
%
%         The 1. Vector defines the Direction of the Glider,
%         the 2. Vector defines the Z-Axis relative to Direction.
%
%         default:  [ 1  0  0 ; 0  0  1 ] Head along X-Axis, no Tilt
%
%--------------------------------------------------------
%
%  V    XYZ-Coordinates of Vertices
%
%  F    VerticeIndex for Faces
%
%  D    Index: 1 Corpus
%              2 Stripes
%              3 Head
%              4 Tail
%              5 Wings
%              6 Fane
%
%  M    ColorMap corresponding to D
%
%  C    TrueColor CData for Faces
%
%--------------------------------------------------------
% Example:
%
% [v,f,d,m] = glider;
%
% figure('ColorMap',m);
%
% axes( 'dataaspectratio' , [ 1  1  1 ] , ...
%       'view'            , [ -30  30 ] , ...
%       'nextplot'        , 'add'  );
%
% patch(     'vertices' , v        , ...
%               'faces' , f        , ...
%     'facevertexcdata' , d        , ...
%        'cdatamapping' , 'direct' , ...
%           'facecolor' , 'flat'   , ...
%           'edgecolor' , 'none'          )
%
% h = light('style','infinite','color','w');
%
% lightangle(h,60,60)
%
%--------------------------------------------------------
% Or try this:
%
% glider([100  100 20])
% glider([100  100  2])
%
%--------------------------------------------------------
% GliderPumping
%
%  r = ( 2 : .2 : 20 );
%
%  [v,f,d,m] = glider([100 100 r(1)]);
% 
%  fig = figure('units'    , 'pixels' , ...
%               'position' , [ 500 300 480 360 ] , ...
%               'MenuBar'  , 'none' , ...
%               'ToolBar'  , 'none' , ...
%               'ColorMap' ,  m           );
% 
%  axe = axes( 'units'      , 'normalized' , ...
%         'position'        , [ 0 0 1 1 ] , ...
%         'dataaspectratio' , [ 1  1  1 ] , ...
%         'view'            , [ -30  30 ] , ...
%         'xlim'            , [ -40  20 ] , ...
%         'ylim'            , [ -80  80 ] , ...
%         'zlim'            , [ -10  10 ] , ...
%         'visible'         , 'off'       , ...
%         'nextplot'        , 'add'  );
% 
%  gld = patch( 'vertices' , v        , ...
%                  'faces' , f        , ...
%        'facevertexcdata' , d        , ...
%           'cdatamapping' , 'direct' , ...
%              'facecolor' , 'flat'   , ...
%              'edgecolor' , 'none'   , ...
%              'clipping'  , 'off'           );
% 
%  lgt = light('style','infinite','color','w');
% 
%  lightangle(lgt,60,60)
%
%  n = size(r,2);
%  
%  M = getframe(fig); M = M(1,ones(1,n+2));
%
%  for ii = 2 : n
%      v = glider([100 100 r(ii)]);
%      set(gld,'vertices',v);
%      M(ii+1) = getframe(fig);
%  end
%
%  M(n+2) = M(n+1);
%
%  delete([lgt gld]);
%
%  movie(axe,M,-100,12);
%


%********************************************************
% Geometry

r =   5;  % Radius
l = 100;  % Length

hf =  2.0; % NormLen of Head versus Radius
tf =  4.0; % NormLen of Tail versus Radius

lf = [ 1  1.5 ];  % Minimum NormLen for Head/Tail CylinderLength Changed

z  = l - r * ( hf + tf );      % Length  of Cylinder

zl = 2;                        % Minimum Length versus Radius

sp = [ 0.00 0.45 0.65 1.00 ];  % NormPos of Stripes versus CylinderLength
sw =   0.02;                   % NormWdt of Stripes versus CylinderLength

wl = 3.0;  % NormLen of WingEdge versus Diameter
wp = 0.64; % NormPos of WingEnd  versus CylinderLength
wb = 0.12; % NormWdt of WingBase versus CylinderLength
we = 0.08; % NormWdt of WingFane versus CylinderLength
wa = 60;   % WingAng of BackEdge

lb = [ 1  3 ];    % Limits for NormWdt of WingBase WB, 
                  % relative to Radius/CylinderLength

ll = [ 1/5 2/3 ]; % Limits of WingLength, relative to Length/Diameter

fh = 1.3;  % NormHgt of FaneEdge  versus Diameter from Center
fp = 0.5;  % NormPos of FaneFront versus TailLength
fb = 0.35; % NormWdt of FaneBase  versus TailLength
fe = 0.25; % NormWdt of FaneTop   versus TailLength
fa = 70;   % FaneAng of FrontEdge

dg = 10;   % Resolution [deg]

% ColorMap

m = [ 1.0  1.0  0.0     % Corpus
      0.0  0.0  0.5     % Stripes
      0.0  0.0  0.0     % Head
      1.0  0.8  0.0     % Tail
      1.0  0.7  0.0     % Wings
      1.0  1.0  0.5 ];  % Fane

%********************************************************

v = zeros(0,3);  % Vertices
f = zeros(0,4);  % Faces
d = zeros(0,1);  % FaceColorIndex
c = zeros(0,3);  % FaceColorRGB

%********************************************************
% Check Input for Resolution/Length

ok = ( nargin >= 1 );
if ok
   nn = prod(size(n));
   ok = ( isnumeric(n) & ( nn <= 3 ) );
   if ~ok
       warning('Invalid Input for Resolution/Length.');
   else
       ok = ~isempty(n);
   end
end

if ok
   dg = ceil(360/n(1));
   if     nn == 2
          r  = r * n(2) / l;
          l  = n(2);
          z  = l - r * ( hf + tf );
   elseif nn == 3
          l  = n(2);
          r  = n(3);
          z  = l - r * ( hf + tf );
   end
end

%********************************************************
% Check GliderDimensions

%------------------------------------------
% Check Length of Cylinder

zm = zl * r;      % ZMin
dz = zm - z;

if dz > 0

   ff = ( l - zm ) / r; 
   sf = hf + tf;
   hf = hf * ff/sf;
   tf = tf * ff/sf;

   if any( [ hf  tf ] < lf )
      hf = max(hf,lf(1));
      tf = max(tf,lf(2));
      r = l / ( zl + hf + tf );  % r(z==zl*r)
      warning('Reduce large Radius.');
   end

   z  = l - r * ( hf + tf );

end

%------------------------------------------
% Check Extension of Wings 

lb = lb * r/z;  % Limits for NormWdt of WingBase WB

%---------------------------------------
if     wb < lb(1)         % wider
%---------------------------------------

       % Check surrounding Stripes
       if isempty(sp)
          ok = [];
       else
          ok = ( ( 0 < sp ) & ( sp < 1 ) );
          if any(ok)
             ok = find(ok);
             sb = sp(ok) - sw/2;
             se = sp(ok) + sw/2;
             se = se + ( sw - se ) .* ( sb < 0 );
             sb = max(sb,0);
             sb = sb + ( 1-sw - sb ) .* ( se > 1 );
             se = min(se,1);
             kk = 1 * ( sb >= wp ) + 2 * ( se <= wp-wb );
             jj = ~( kk == 0 );
             if any(jj)
                jj = find(jj);
                kk =   kk(jj);
                ok =   ok(jj);
             else
                ok = [];
             end   
          end
       end

       sw = sw * lb(1)/wb;
       we = we * lb(1)/wb;
       wp = wp - wb/2 + lb(1)/2;
       wb = lb(1);
       wp = min(wp,1-3*sw);

       if ~isempty(ok)
           jj = ( kk == 1 );
           if any(jj)
              jj = find(jj);
              sp(ok(jj)) = max(sp(ok(jj)),wp+sw/2);
           end 
           jj = ( kk == 2 );
           if any(jj)
              jj = find(jj);
              sp(ok(jj)) = min(sp(ok(jj)),wp-wb-sw/2);
           end
       end

%---------------------------------------
elseif wb > lb(2)         % narrow
%---------------------------------------

       we = we * lb(2)/wb;
       wp = wp - wb + lb(2);
       wb = lb(2);

%---------------------------------------
end
%---------------------------------------

%---------------------------------------
% Length of Wings
 
ll = ll * l / (2*r);

wl = min( max( wl , ll(1) ) , ll(2) );


%********************************************************
% Check Input for Direction

T = [];       % Transformation-Matrice

ok = ( nargin >= 2 );
if ok
   s  = size(w);
   p  = prod(s);
   ok = ( isnumeric(w) & any( p == [ 3  6 ] ) );
   if ok
      if p == 3
         w = w(:)';
      else
         ok = ( isequal(s,[2 3]) | isequal(s,[3 2]) );
         if ok & ( s(1) == 3 )
            w = permute(w,[2 1]);
         end
       end
   end
   if ~ok
       warning('Invalid Input for Direction.');
   end
end

if ok

   n  = sqrt(sum(w.^2,2));
   ok = ~( n == 0 );
   if any(~ok)
      warning('ZERO-Input for Direction.');
   end

   if ~any(ok)
       ok = 0;
   elseif ~all(ok)
       def     = [ 1 0 0 ; 0 0 1 ];
         jj    = find(~ok);
       w(jj,:) = def(jj,:);
       n(jj)   = 1;
   end

   w = w ./ n(:,[1 1 1]);

   ex = w(1,:);
   ey = ex([2 1 3]) .* [ -1  1  0 ];

   if all(ey == 0)
      ey = [ 0  1  0 ];
   else
      ey = ey ./ sqrt(sum(ey.^2));
   end

   ez = cross(ex,ey);

   T = cat( 1 , ex , ey , ez );

   if size(w,1) == 2

      ex = [ 1  0  0 ];
      ez = w(2,:);

      if ez(3) == 0
         ex = ez([2 1 3]) .* [ 1  -1  0 ];
      else
         ex(3) = -ez(1)/ez(3);
         ex    = ex ./ sqrt( 1 + ex(3)^2 );
      end

      ey = cross(ez,ex);

       T = cat(1,ex,ey,ez) * T;

   end

end

%********************************************************

nf = ceil(360/dg);
nv = nf + 1;

pp = pi/180 * linspace(0,360,nv)';

yp = r * cos(pp);
zp = r * sin(pp);

fv = ( 1 : nf )';
fv = cat( 2 , fv , fv+1 );

fv = cat( 2 , fv , fv(:,[2 1])+nv );

%********************************************************
% Zylinder with Stripes

if isempty(sp)

   xp = [ 0  1 ];

else

   sp = sp(:)';

   ns = size(sp,2);

   sb = sp - sw/2;
   se = sp + sw/2;

   se = se + ( sw - se ) .* ( sb < 0 );
   sb = max(sb,0);

   sb = sb + ( 1-sw - sb ) .* ( se > 1 );
   se = min(se,1);

   if any( sb(2:ns) < se(1:(ns-1)) )
      warning('Overlapping Stripes');
      ok = ones(1,ns);
      zz  = 1;
      while zz < ns
            jj = ( sb((zz+1):ns) < se(zz) );
            if any(jj)
               jj = find(jj) + zz;
               ok(jj) = 0;
               zz = max(jj);
            end
            zz = zz + 1;
      end
      ns = sum(ok);
      ok = find(ok);
      sp = sp(ok);
      sb = sb(ok);
      se = se(ok);
   end

   xp = cat( 1 , sb , se );

   xp = xp(:)';

   if ~( xp(1) == 0 )
       xp = cat( 2 , 0 , xp );
       sp = [];      % Start with Zylinder !!!
   end

   if ~( xp(ns) == 1 )
       xp = cat( 2 , xp , 1 );
   end

end

xp = z * xp;

np = size(xp,2);

n = size(v,1);

v = cat( 1 , v , cat( 2 , reshape(xp(ones(1,nv),:),nv*np,1) , ...
                          reshape(yp(:,ones(1,np)),nv*np,1) , ...
                          reshape(zp(:,ones(1,np)),nv*np,1)       ) );

np = ( np - 1 );

ff = nv * ( cumsum(ones(nf,4,np),3) - 1 );

ff = ff + fv(:,:,ones(1,np)) + n;

ff = permute(ff,[1 3 2]);

f = cat( 1 , f , reshape(ff,nf*np,4) );

cc = cumsum(ones(nf,np),2) + isempty(sp);
cc = mod(cc,2) + 1;

cc = cc(:);

d = cat( 1 , d , cc );

%********************************************************
% Elliptic Head and Tail
%
%   rad(phi) = rx * ry / sqrt( ry² + (rx²-ry²) * sin²(phi-rot) )
%

np = ceil(90/dg) + 1;

pr = pi/180 * linspace(0,90,np);

%------------------------------------------------------------------
% Head

rr = r * hf;

rr = r * rr ./ sqrt( rr^2 + (r^2-rr^2) * sin(pr).^2 );

yp = cos(pp) * ( rr .* cos(pr) );
zp = sin(pp) * ( rr .* cos(pr) );

xp = -( rr .* sin(pr) );

np = size(xp,2);

n = size(v,1);

v = cat( 1 , v , cat( 2 , reshape(xp(ones(1,nv),:),nv*np,1) , ...
                          reshape(yp,nv*np,1) , ...
                          reshape(zp,nv*np,1)       ) );

np = ( np - 1 );

ff = nv * ( cumsum(ones(nf,4,np),3) - 1 );

ff = ff + fv(:,:,ones(1,np)) + n;

ff = permute(ff,[1 3 2]);

f = cat( 1 , f , reshape(ff,nf*np,4) );

d = cat( 1 , d , 3*ones(nf*np,1) );

%------------------------------------------------------------------
% Tail

rr = r * tf;

rr = r * rr ./ sqrt( rr^2 + (r^2-rr^2) * sin(pr).^2 );

yp = cos(pp) * ( rr .* cos(pr) );
zp = sin(pp) * ( rr .* cos(pr) );

xp = z + ( rr .* sin(pr) );

np = size(xp,2);

n = size(v,1);

v = cat( 1 , v , cat( 2 , reshape(xp(ones(1,nv),:),nv*np,1) , ...
                          reshape(yp,nv*np,1) , ...
                          reshape(zp,nv*np,1)       ) );

np = ( np - 1 );

ff = nv * ( cumsum(ones(nf,4,np),3) - 1 );

ff = ff + fv(:,:,ones(1,np)) + n;

ff = permute(ff,[1 3 2]);

f = cat( 1 , f , reshape(ff,nf*np,4) );

d = cat( 1 , d , 4*ones(nf*np,1) );

%********************************************************
% Wings
%
% wp  NormPos of WingEnd  versus Length
% wb  NormWdt of WingBase versus Length
% we  NormWdt of WingFane versus Length
% wl  NormLen of WingEdge versus Length
% wa  WingAng of BackEdge
%

wa = wa * pi/180;

wl = wl * 2*r;

x0 = wl * [ 0  cos(wa) ];
x1 = x0 - z * [ wb  we ];

xx = cat( 1 , x0 , x1 );
xx = xx(:);
xx = z * wp + xx([1 2 4 3]);

yy = [ 0   sin(wa) ];
yy = yy([1 1 2 2]);
yy = r + wl * yy(:);

%------------------------------------------------------------------
% Right

n = size(v,1);

v = cat( 1 , v , cat( 2 , xx , yy , zeros(4,1) ) );

f = cat( 1 , f , (1:4)+n );

d = cat( 1 , d , 5 );

%------------------------------------------------------------------
% Left

n = size(v,1);

v = cat( 1 , v , cat( 2 , xx , -yy , zeros(4,1) ) );

f = cat( 1 , f , (1:4)+n );

d = cat( 1 , d , 5 );

%********************************************************
% Fane
%
% fp  NormPos of FaneFront versus TailLength
% fb  NormWdt of FaneFront versus TailLength
% fe  NormWdt of FaneFront versus TailLength
% fh  NormHgt of FaneEdge  versus Diameter from Center
% fa  FaneAng of FrontEdge
%

fa = fa * pi/180;

x0 = 2*r * fh * [ 0  cos(fa) ];
x0 =  x0 + tf * r * fp;
x1 =  x0 + tf * r * [ fb  fe ];

xx = cat( 1 , x0 , x1 );
xx = xx(:);
xx = z + xx([1 2 4 3]);

zz = [ 0   sin(fa) ];
zz = zz([1 1 2 2]);
zz = 2*r * fh * zz(:);

%------------------------------------------------------------------

n = size(v,1);

v = cat( 1 , v , cat( 2 , xx , zeros(4,1) , zz ) );

f = cat( 1 , f , (1:4)+n );

d = cat( 1 , d , 6 );

%********************************************************
% TrueColor CData

nc = size(m,1);

c = m(cat(2,d,d+nc,d+2*nc));

%********************************************************
% Center in X

v(:,1) = -( v(:,1) - z*(wp-wb/2) );

%********************************************************
% Transformation

if ~isempty(T)
    v = v * T;
end

%********************************************************
% Plot

if nargout == 0
 
  figure('ColorMap',m);
 
  axes( 'dataaspectratio' , [ 1  1  1 ] , ...
                  'view'  , [ -30  30 ] , ...
        'nextplot'        , 'add'  );
 
  patch(     'vertices' , v        , ...
                'faces' , f        , ...
      'facevertexcdata' , c        , ...
         'cdatamapping' , 'direct' , ...
            'facecolor' , 'flat'  , ...
            'edgecolor' , 'none')
 
  h = light('style','infinite','color','w');
 
  lightangle(h,60,30)

  clear v f d m c

end
