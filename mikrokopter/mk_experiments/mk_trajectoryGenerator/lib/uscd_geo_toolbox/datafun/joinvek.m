function [y,l] = joinvek(x,acc,mode)

% JOINVEK  Returns continuous Liness of NaN-Row separated Segments
%
% [ Y , L ] = JOINVEK( X , Accuracy , Mode ) 
%
% Y = { N by 1 }  CellArray with continous Segments
% L = [ N by 1 ]  Length of Segments
%
% Mode = 1  ==> Segments, which are not closed, 
%                 are appended with a NaN-Row.
%
% Older, faster Version of JOINVEC
%

Nin = nargin;

if Nin < 1
   y = {};
   l = [];
   return
end

if Nin < 2
   acc = [];
end

if Nin < 3
   mode = 0;
end

if ~( isnumeric(x) & ( ndims(x) == 2 ) )
    error('X must be a 2-Dimensional numeric.');
end

if isempty(acc)
   acc = 1e-6;
elseif ~( isnumeric(acc) & ( prod(size(acc)) == 1 ) )
   error('Accuracy must be a 2-Dimensional numeric.');
else
   acc  = abs(acc);
end

mode = isequal(mode,1);

%---------------------------------------------------

if isempty(x)
   y = {x};
   l = 0;
   return
end

%---------------------------------------------------
% Remove following, leading, ending  NaN-Rows

s2 = size(x,2);

n = find( sum(isnan(x),2) == 2 );

if ~isempty(n)
    [i0,lg] = ind2grp(n);
    i1 = ( ( i0 == 1 ) | ( i0+lg-1 == size(x,1) ) );
    ok = ( i1 | ( lg > 1 ) );
    if any(ok)
       ok = find(ok);
       i0 = i0(ok);
       lg = lg(ok);
       i1 = i1(ok);
       i0 = i0 + ( i1 == 0 );
       i0 = grp2ind(i0,lg);
       x(i0,:) = [];
       if isempty(x) 
          n = [];
       else
          n = find( sum(isnan(x),2) == s2 );
       end
    end
end

if isempty(n)
   y = {x};
   l = size(x,1);
   return
end

%---------------------------------------------------

i0 = cat( 1 , 1 , n+1 );          % StartIndex
i1 = cat( 1 , n-1 , size(x,1) );  %   EndIndex

x0 = x(i0,:);    % StartCoordinates
x1 = x(i1,:);    %   EndCoordinates

n0 = isnan(x0);
n1 = isnan(x1);

isn = any( n0(:) | n1(:) );
if isn
   warning('NaN''s in Start- or EndCoordinates');
end

%-----------------------------------------

n = size(n,1) + 1;

y =  cell(n,1);  % Segents
l = zeros(n,1);  % Length of Segment

z = 0;           % SegmentCounter

%-----------------------------------------
% Fill Closed Segments first

ok = ( abs(x0-x1) <= acc );
if isn
   ok = ( ok | ( n0 & n1 ) );
end

ok = ( sum(ok,2) == s2 );

if any(ok)
   z = sum(ok);
   kk = find(ok);
   for ii = 1 : z
       jj = ( i0(kk(ii)) : i1(kk(ii)) );
       y{ii} = x(jj,:);
       l(ii) = size(jj,2);
   end
end

%-----------------------------------------
% Follow other Segments

while 1

   kk = ( ok == 0 );   % Ok for free Segments

   nk =  sum(kk);      % Number of free Segments

   if ( nk == 0 )
      break
   end

   kk = find(kk);      % Index of free Segments

   ik = 1;             % Start with 1. free Segment

   k1 = kk(ik);        % Index in OK, I0, I1, X0, X1

   z = z+1;            % Erase SegmentCounter

   jj = ( i0(k1) : i1(k1) );
 y{z} = x(jj,:);
 l(z) = size(jj,2);

 ok(k1) = 1;

   y0 = x0(k1,:);
   y1 = x1(k1,:);

   if isn
      m0 = n0(k1,:);
      m1 = n1(k1,:);
   end

   if nk == 1
      break
   end

   while 1

      kk(ik) = [];
      nk     = nk - 1;

      if ( nk == 0 )
         break
      end

      %-------------------------------------------------------------
      % Check with Distance
      % End-Start / End-End / Start-End / Start-Start

      d = cat( 1 , ( y1(ones(1,nk),:) - x0(kk,:) ) , ...
                   ( y1(ones(1,nk),:) - x1(kk,:) ) , ...
                   ( y0(ones(1,nk),:) - x1(kk,:) ) , ...
                   ( y0(ones(1,nk),:) - x0(kk,:) )        );

      d = abs(d);
      c = ( d <= acc );

      k = ( sum(c,2) == s2 );

      %-------------------------------------------------------------
      % Check with NaN's at Start/End

      if ~any(k) & isn

          m = cat( 1 , ( m1(ones(1,nk),:) & n0(kk,:) ) , ...
                       ( m1(ones(1,nk),:) & n1(kk,:) ) , ...
                       ( m0(ones(1,nk),:) & n1(kk,:) ) , ...
                       ( m0(ones(1,nk),:) & n0(kk,:) )       );

          c = ( c | m );

          d(find(m)) = 10*s2*acc;

          k = ( sum(c,2) == s2 );

      end

      %-------------------------------------------------------------
      % Find connected Segment

      if ~any(k)
          break
      end

      k = find(k);
      d = sum( d(k,:) , 2 );

      [d,ii] = min(d);

       ii =  k(ii);                  % Index in Check [ 1 .. 4*NK ]
       ik = ii - nk * floor(ii/nk);  % Index in kk
       ik = ik + nk * ( ik == 0 );
       k1 = kk(ik);                  % Index in OK, I0, I1, X0, X1

      %-------------------------------------------------------------
      % Check Connection-Type
      % End-Start / End-End / Start-End / Start-Start

       ii = ceil(ii/nk);  % [ 1 .. 4 ]
 
       jj = ( i0(k1) : i1(k1) );

       % Check for Flip !!!!
       if mod(ii,2) == 0
          jj = jj(end:-1:1);
          nn = n0(k1,:);
          xx = x0(k1,:);
          x0(k1,:) = x1(k1,:); x1(k1,:) = xx;
          if isn
             n0(k1,:) = n1(k1,:); n1(k1,:) = nn;
          end
       end

       if ii <= 2
          y{z} = cat( 1 , y{z} , x(jj,:) );
          y1   = x1(k1,:);
          if isn
             m1   = n1(k1,:);
          end
       else
          y{z} = cat( 1 , x(jj,:) , y{z} );
          y0   = x0(k1,:);
          if isn
             m0   = n0(k1,:);
          end
       end

       l(z) = l(z) + size(jj,2);

       ok(k1) = 1;

   end

end

y = y(1:z);
l = l(1:z);

if ~mode
    return
end

nn = NaN * zeros(1,s2);

for ii = 1 : z
    if ~all( abs(y{ii}(1,:)-y{ii}(l(ii),:)) <= acc )
        y{ii} = cat( 1 , y{ii} , nn );
    end
end
