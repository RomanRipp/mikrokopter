function [y,z] = splinen(x,n,p,d);

% SPLINEN   Calculates Spline
%
% Y = SPLINEN( X , N , P , D )
%
% X   Vector or Matrice 
% 
% N   Number of new Points
%        default: N = 100
%
% P   Potenz along Distance between points
%     The Distance is  SQRT(SUM(X.^2,K)),
%      where K is the opposite Dimension of D
%        default: P = 1
%
% D   Dimension to work along
%        default: longest Dimension of X
%
% Use a negative Number of Points to iterate
%  the DistanceWeight 3 Times with SQRT(SUM(Y.^2,K))
%  where at 1. Iteration Y is set to X
%
%--------------------------------------------------
% Example for CurveSplines:
%--------------------------------------------------
% Random Points, different Values for Potenz
%--------------------------------------------------
%
%  x = rand(5,2); 
%  a = splinen(x,100,1.0,1);
%  b = splinen(x,100,1.5,1);
%  c = splinen(x,100,0.5,1);
%
%  figure, hold on, box on
%
%  plot(x(:,1),x(:,2),'k*-');
%  plot(a(:,1),a(:,2),'r-');
%  plot(b(:,1),b(:,2),'g-');
%  plot(c(:,1),c(:,2),'b-');
%
%  title('CurveSpline')
%  legend('Original','P = 0.5', 'P = 1.0' , 'P = 1.5')
%
%--------------------------------------------------
% Regular Triangle
%--------------------------------------------------
%
% m = 3;   % Corners (Triangle)
% n = 100; % Number of SplinePoints between Corners
% u = 3;   % Number of Turns to smooth the Spline at start and end
%
% p = 360/m * ( 0 : u*m ) - 90 - 360/m/2;     % Angles of Corners
%
% xy = [ cos(p*pi/180) ; sin(p*pi/180) ]; % Regular Triangles
%
% a = splinen( xy , u*m*n+1 , 1 , 2 );    % All Turns
%
% xy = xy( : , ( 0 : m ) + 1 );           % 1. Turn only
% b = splinen( xy ,   m*n+1 , 1 , 2 );
%
% figure, hold on, box on, axis equal
%
% plot(xy(1,:),xy(2,:),'k*-');
% plot(xy(1,1),xy(2,1),'r*','linestyle','none');
%
% leg = { 'Original' 'StartPoint' };
%
% c = 'rgb';
% for ii = 1 : u
%     jj = ( 0 : m*n ) + (ii-1)*m*n + 1;
%     plot(a(1,jj),a(2,jj),'-','color',c(ii))
%     leg = cat( 2 , leg , {sprintf('%.0f. Turn',ii)} );
% end
%
% plot(b(1,:),b(2,:),'k--')
% leg = cat( 2 , leg , {'single Turn'} );
%
% title(sprintf('Spline of %.0f anticlockwise Turns',u));
% legend(leg{:});
%
%--------------------------------------------------
% Regular Star, requires STAR
%--------------------------------------------------
%
% m = 7;   % Corners
% n = 100; % Number of SplinePoints between Corners
% u = 3;   % Number of Turns to smooth the Spline at start and end
%
% x = star(m);  % Star, 2*m+1 points (closed polygon)
%
% m  = 2 * m;
% ii = ( 0 : u*m ) + 1;       % 3 Turns
% ii = ii - m * floor(ii/m);
% ii = ii + m * ( ii == 0 );
%
% y = splinen( x(ii,:) , u*m*n+1 , 1 , 1 );
%
% ii = ( 0 : m*n ) + m*n + 1; % 2. Turn
%
% y = y(ii,:);
%
% figure, hold on, box on, axis equal
%
% plot(x(:,1),x(:,2),'k*-');
% plot(y(:,1),y(:,2),'r');
%
% title(sprintf('Spline on %.0f-Star',m));
%

y = [];
z = [];

%************************************************

Nin = nargin;

if Nin < 1
   return
end

if Nin < 2
   n = [];
end

if Nin < 3
   p = [];
end

if Nin < 3
   d = [];
end

%------------------------------------------------

s = size(x);

if isempty(p)
   p = 1;
end

if isempty(d)
   [h,d] = max(s([1 2]));
end

%------------------------------------------------

msg = cell(0,1);

ret = isempty(x);
flp = 0;

if ~ret
    if ~( isnumeric(x) & ( ndims(x) == 2 ) )
        msg = cat( 1 , msg , {'XY must be 2D-Numeric.'} );
    else
        flp = ~( d == 2 );
        if flp
           x = permute(x,[2 1]);
        end
    end
end

ok = ( isnumeric(n) & ( prod(size(n)) == 1 ) );
if ok
   ok = (  mod(n,1) == 0 );
end

if ~ok
    msg = cat( 1 , msg , {'N must be a single positive Integer.'} );
end

if ~( isnumeric(p) & ( prod(size(p)) == 1 ) )
    msg = cat( 1 , msg , {'P must be a single Numeric.'} );
end

%------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Invalid Inputs.\n%s',msg));
end

%************************************************

if ~ret

    dd = sqrt(sum(diff(x,1,2).^2,1));

    jj = ( dd <= 1e3*eps );

    ret = all(jj);

    if ~ret 
        if any(jj)
           % Remove Duplicate Points
             jj    = find(jj);
          dd(jj)   = [];
         x(:,jj+1) = [];
        end
        jj  = isfinite(dd);
        ret = ~any(jj);
    end

end

if ret
   m = size(x,1);
   y = zeros(m*(1-flp),m*flp);
   z = zeros(1*(1-flp),1*flp);
   return
end

%---------------------------------------------------
% Use only finite Intervalls

if all(jj)

   [y,z] = clcspl(x,n,p);

else

   dd(find(~jj)) = 0;

   jj = find(jj);

   i0 = cat( 2 , 1 , find( diff(jj,1,2) > 1 )+1 , size(jj,2)+1 );

    m = size(i0,2) - 1;   % GroupNumber

   ll = diff(i0,1,2);     % GroupLenght
   i0 = jj(i0(1:m));      % GroupStart
   i1 = i0 + ll - 1;      % GroupEnd 

   ll = cumsum(dd,2);
   ll = ll(i1);
   ll(2:m) = ll(2:m) - ll(1:(m-1));  % GroupSum

   sg = sign(n);
    n =  abs(n);

   nn = ceil( n * ll / sum(ll) );    % GroupResolution
   nn =  max( nn , 2 );

   off = sum(nn,2) - n;

    jj = ( nn > 2 );

   while ( off > 0 ) & any(jj)

         jj = find(jj);
         dv = off * nn(jj) / sum(nn(jj));

         dd = floor(dv);

         if all( dd == 0 )
            dd = round(dv);
            if all( dd == 0 )
               dd = ceil(dv);
            end
            if sum(dd) > off
               dv = sum(dd) - off;
               [h,si] = sort(nn(jj));
               dd(si(1:dv)) = dd(si(1:dv)) - 1;
            end
         end

         nn(jj) = nn(jj) - dd;
         off    = sum(nn,2) - n;
            jj  = ( nn > 2 );

   end

   mm = sum(nn) + (m-1);

    y = NaN * zeros( size(x,1) , mm );
    z = NaN * zeros( 1 , mm );
    
   j0 = cumsum(cat(2,1,nn+1),2);

   i1 = i1 + 1;

   jj = ( 2 : m );

   z0 = cat( 2 , 0 , i0(jj) - 1 );

   z(j0(jj)-1) = z0(jj);

   for ii = 1 : m

       jj = j0(ii) + ( 1 : nn(ii) ) - 1;

       [y(:,jj),z(jj)] = clcspl(x(:,i0(ii):i1(ii)),sg*nn(ii),p);

       z(jj) = z(jj) + z0(ii);

   end

end

%---------------------------------------------------

if flp
   y = permute(y,[2 1]);
   z = permute(z,[2 1]);
end

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [y,z] = clcspl(x,n,p);

m = size(x,2);

nn = 2 * ( n < 0 );
n  = abs(n);

ind = ( 0 : (n-1) ) / (n-1);

if     isempty(x)
   y = zeros(size(x))
   z = zeros(1,0);
   return
elseif m == 1
   y = x(:,ones(1,n));
   z = ones(1,n);
   return
elseif m == 2
   y = x(:,ones(1,n),:) + ( x(:,2) - x(:,1) ) * ind;
   z = ind+1;
   return
end

d = sqrt(sum(diff(x,1,2).^2,1)) .^ p;

d = cumsum(cat(2,0,d),2);

%-----------------------------------------

k = size(x,1);

y = zeros(k,n);

g = ind * d(m);

for ii = 1 : k
    y(ii,:) = spline(d,x(ii,:),g);
end

z = interp1(d,(1:m),g);

if nn == 0
   return
end

for ii = 1 : nn

    l = cat(2,zeros(k,1),diff(y,1,2));
    l = cumsum(sqrt(sum(l.^2,1)),2);

    f = l(n) * ind;

    g = interp1(l,g,f);

    for ii = 1 : k
        y(ii,:) = spline(d,x(ii,:),g);
    end

    z = interp1(d,(1:m),g);

end
