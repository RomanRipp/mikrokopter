function [x,y,cc,yl] = conegrid(c,xl,yl,n);

% CONEGRID  creates a conic Grid on a Sphere
%
% [Lon,Lat] = CONEGRID( Center , Angle , Distance , n );
%
% Center = [ CenterLon CenterLat [Direction] ]
%
% Direction: clockwise to North
%
% Direction + i : Start with Center at DistLim(1)
%                 negative DistLim(1) : shift Start back
%
% Resolution = Nad | [ Na  Nd ]
%              used for AngleLimit / DistanceLimit
%
% Angle     Vector or Matrice for Angles, Zero in Line of Direction
%
%           two-element Vector: AngleLimit [ AngleMin AngleMax ]
%
%           single Value: AngleLimit = [ -Angle/2  +Angle/2 ]
%
% Distance  Vector or Matrice for Distances to Center, foreward in Direction
%           Units: [deg] == [60*nm]
%
%           two-element Vector: DistanceLimit [ DistMin DistMax ]
%
%           single Value: DistanceLimit = [ 0  Distance ]
%
%           a negative Value for DistMin and Non-zero imaginary Part of Direction
%             shifts the Center by 2*DistMin (back in Direction)
% 
%----------------------------------------------------------------
%
% see also: SPH_PROJ, DST_PROJ, OBS_PROJ
%
%----------------------------------------------------------------
% Example: ConeGrids centered on Europe
%----------------------------------------------------------------
%
% area = [ -15  25  35  62 ]; cc = [ 5  48 ];
%
%%----------------------------------------------------------------
%
% load coast
% 
% ii = find(~inpolygon(long,lat,area([1 2 2 1]),area([3 3 4 4])));
%
% long(ii) = NaN; lat(ii) = NaN;         % NaN if outside Area
%
% ii = find( isnan(long) | isnan(lat) );
% jj = find( diff(ii) == 1 ) + 1;        % Duplicate NaN's
%
% long(ii(jj)) = []; lat(ii(jj)) = [];
%
% xf = interp1(1:5,area([1 2 2 1 1]),linspace(1,5,400));
% yf = interp1(1:5,area([3 3 4 4 3]),linspace(1,5,400));
% 
%%----------------------------------------------------------------
%
% figure('position',[200 200 400 400]), hold on
%
% plot(long,lat,'k-'); plot(xf,yf,'k-','linewidth',3);
%
% 
% [x,y] = conegrid([ cc  90 ],[-30 30],(-10:2:10),13);
% plot(x',y','b.-');
%
% [x,y] = conegrid([ cc  30 ],[-20 20],[4 12],9);
% plot(x,y,'g-');
% 
% [x,y,c] = conegrid([ cc 170+i ],[-10 10],[-3.5 10.5],[5 15]);
% surface(x,y,0*y,'edgecolor','b','facecolor','none','linestyle','--');
% plot(c(1),c(2),'b*');
%
% [x,y,c] = conegrid([ cc 170+i ],[-10 10],[2 10],[5 9]);
% surface(x,y,0*y,'edgecolor','r','facecolor','none','linewidth',2);
% plot(c(1),c(2),'r*');
%
%%----------------------------------------------------------------
%% Projection on Sphere for correct Aspect
%
% lim = NaN * ones(1,4); op = [ -1  1  -1  1 ];
%
% for h = get(gca,'children')'
% 
%     x = get(h,'xdata'); y = get(h,'ydata');
%     ii = find( isnan(x) | isnan(y) );
%
%     [x,y,z] = sph_proj([0 90 0],x,y);
%
%     x(ii) = NaN; y(ii) = NaN;
%
%     lm = [ min(x(:)) max(x(:)) min(y(:)) max(y(:)) ];
%     lim = op .* max( lm.*op , lim.*op );
%
%     set( h , 'xdata' , x , 'ydata' , y , 'zdata' , z );
%
% end
%
% set( gca , 'xlim' , lim([1 2]) , 'ylim' , lim([3 4]) , ...
%            'view' , cc , 'dataaspectratio' , [ 1 1 1 ] , ...
%            'position' , [ .02 .02 .96 .96] , 'visible' , 'off' );
%
%%----------------------------------------------------------------
%

if ( prod(size(c)) == 2 )

   c = cat( 1 , c(:) , 180 * ( sign(c(2)) == 1 ) );

end

if prod(size(n)) == 1
   n = [ n  n ];
end

%******************************************************

%-------------------------------------
% Check for Y

if prod(size(yl)) == 1
   yl = [ 0  yl ];
end

bd = ~( imag(c(3)) == 0 );

c(3) = real(c(3));

if prod(size(yl)) == 2

   if bd % & ~( yl(1) == 0 )

      %%% !!! Change Direction if Center changes !!! %%%
      
      p180 = pi/180;

      isn = ( yl(1) < 0 );             % Check for Negative

      dev = ( 1 + isn ) * abs(yl(1));  % Shift of Center

       yl = yl - 2 * yl(1) * isn;


      [x,y,z] = sph_proj([0 90-dev 180],0,90);     

           w =   pi - c(3)   * p180;
           a = ( 90 - c(2) ) * p180;   % Pole --> Old

      [c(1),c(2)] = sph_proj(c,x,y,z);           

           b = ( 90 - c(2) ) * p180;   % Pole --> New

           g =  dev * p180;            % New  --> Old

           cw = ( cos(a) - cos(b)* cos(g) ) / ( sin(b) * sin(g) );
           sw =   sin(w) * sin(a) / sin(b);

          c(3) = atan2(sw,cw) / p180;  % New Direction !!!

   end

   y = 90 - linspace(yl(1),yl(2),n(2))';

else

   y = 90 - yl;

end

%-------------------------------------
% Check for X

if prod(size(xl)) == 1
   xl = xl/2 * [ -1  1 ];
end

if prod(size(xl)) == 2
   x = -linspace(xl(1),xl(2),n(1));
else
   x = -xl;
end

%******************************************************

cc = c;

c(3) = c(3) + 180;

[x,y,z] = sph_proj([0 90],x,y);

[x,y] = sph_proj(c,x,y,z);
