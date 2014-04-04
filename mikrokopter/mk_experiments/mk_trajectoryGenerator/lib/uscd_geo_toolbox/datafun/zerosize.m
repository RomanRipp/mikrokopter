function c = zerosize(c,n,p)

% ZEROSIZE  Returns Depth of ZeroFields in Matrice 
%
%   Z = ZEROSIZE( C , N , Period )
%
%
%  C        2-dimensional Matrize
%  N        RecursionDepth to find closed ZERO-Fields in C
%
%  Period = PeriodX + i*PeriodY
%
%  PeriodX/Y = 0 | 1 
%
%  Z       Matrize with same size of C, where the closed ZERO-Fields
%           of C have Non-Zero-Values with there Depth
%
% see also: REPL_NAN
%
% Example:
%
% load('topo.mat','topo');
%
% n = 10;
%
% c = ( topo > 0 );  % Land
% z = zerosize(cat(2,c(:,end-n+1:end),c,c(:,1:n)),(n-1));
% z = z( : , n+1 : end-n );
%
% z(find(z==0)) = n+1;  % Land
% 
% cm = cat( 1 , jet(n) , [0.5 0.5 0.5]);
%
% figure('colormap',cm)
% set(gca,'ydir','normal','box','on','nextplot','add')
%
% image(z,'cdatamapping','direct');
%
% colorbar horiz
%

if isempty(c)
   return
end

if nargin < 2
   n = 10;
end

if nargin < 3
   p = 0;
end

if ~any( c(:) == 0 )
   c = 0*c;
   return
end

%*****************************************************
% Append n Elements at Begin and End if Period

%-----------------------------------------------------
% Period in X

if ~( real(p) == 0 )

   nx0 = size(c,2);

   ix = ( 1-n : nx0+n );
   ix = ix - nx0 * floor( ix / nx0 );  % [ 0 .. nx0-1 ]
   ix = ix + nx0 * ( ix == 0 );

    c = c(:,ix);

end

%-----------------------------------------------------
% Period in Y

if ~( imag(p) == 0 )

   ny0 = size(c,1);

   iy = ( 1-n : ny0+n );
   iy = iy - ny0 * floor( iy / ny0 );  % [ 0 .. ny0-1 ]
   iy = iy + ny0 * ( iy == 0 );

    c = c(iy,:);

end

%*****************************************************
%-----------------------------------------------------

nx = size(c,2);
ny = size(c,1);

n1 = min(n,nx);
n2 = min(n,ny);

c  = double( ~( c == 0 ) );

ix = ( 1 : nx-1 );
iy = ( 1 : ny-1 );

for ii = 1 : max(n1,n2)

    if n1 >= ii 

      c(:,ix+0) = c(:,ix+0) + (ii+1) * ( ( c(:,ix+0) == 0 ) &  ~( c(:,ix+1) == 0 ) );
      c(:,ix+1) = c(:,ix+1) + (ii+1) * ( ( c(:,ix+1) == 0 ) &  ~( c(:,ix+0) == 0 ) );

    end


    if n2 >= ii

      c(iy+0,:) = c(iy+0,:) + (ii+1) * ( ( c(iy+0,:) == 0 ) &  ~( c(iy+1,:) == 0 ) );
      c(iy+1,:) = c(iy+1,:) + (ii+1) * ( ( c(iy+1,:) == 0 ) &  ~( c(iy+0,:) == 0 ) );

    end
 
    if ~any( c(:) == 0 )
       break
    end

end
  
  
c = c + (n+2) * ( c == 0 );

c = c - 1;

%*****************************************************
% Remove n Elements from Begin and End if Period

%-----------------------------------------------------
% Period in X

if ~( real(p) == 0 )

   c  = c(:,(n+(1:nx0)));

end

%-----------------------------------------------------
% Period in Y

if ~( imag(p) == 0 )

   c  = c((n+(1:ny0)),:);

end
