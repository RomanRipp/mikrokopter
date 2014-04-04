function [y,l] = joinvec(x,acc,mode)

% JOINVEC  Returns continuous Lines of NaN-Row separated Segments
%
% [ Y , L ] = JOINVEC( X , Accuracy , Mode ) 
%
% Y = { N by 1 }  CellArray with continous Segments
% L = [ N by 1 ]  Length of Segments
%
% Mode = 1  ==> Segments, which are not closed, 
%                 are appended with a NaN-Row.
%
% Slower but more elegant version of JOINVEK
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

i0 = cat( 1 , 1 , n+1 );         % StartIndex
i1 = cat( 1 , n-1 , size(x,1) ); %   EndIndex

x01 = cat( 3 , x(i0,:) , ...    % StartCoordinates
               x(i1,:)       ); %   EndCoordinates

x01 = permute(x01,[1 3 2]);

n01 = isnan(x01);

isn = any(isnan(n01(:)));

if isn
   warning('NaN''s in Start- or EndCoordinates');
end

%-----------------------------------------

n = size(n,1) + 1;      % SegmentsNumber

z = 0;                  % SegmentCounter

%-----------------------------------------
% Fill Closed Segments first

ok = ( abs(diff(x01,1,2)) <= acc );
if isn
   ok = ( ok | ( n01(:,1,:) & n01(:,2,:) ) );
end

ok = ( sum(ok,3) == s2 );

if any(ok)
   z  = sum(ok);
   nr = cumsum(ok,1);   % SegmentNumber
   nr(find(~ok)) = 0;
else
   nr = zeros(n,1);     % SegmentNumber
end

zz = zeros(n,1);        % Counter in Segment

%-----------------------------------------
% Follow other Segments

while 1

   ok = ( nr == 0 );   % Ok for free Segments

   nk =  sum(ok);      % Number of free Segments

   if ( nk == 0 )      % No Segment free
      break
   end

   ok = find(ok);      % Index of free Segments

   ik = 1;             % Start with 1. free Segment

   k1 = ok(ik);        % Index in Segments

   z  = z + 1;         % Erase SegmentCounter

   nr(k1) = z;

   if nk == 1          % Last free Segment
      break
   end

   cc  = zeros(1,2);   % CatCounter

   y01 = x01(k1,:,:);

   if isn
      m01 = n01(k1,:,:);
   end

   %-----------------------------------------------------------
   % Search other Segments

   while 1

      ok(ik) = [];      % Remove Index of previous Segment 
      nk     = nk - 1;

      if ( nk == 0 )    % No Segment free
         break
      end

      %-------------------------------------------------------------
      % Check with Distance
      % End-Start / End-End / Start-End / Start-Start

      d = y01(ones(1,nk),[2 2 1 1],:) - x01(ok,[1 2 2 1],:);

      d = abs(d);
      c = ( d <= acc );

      %-------------------------------------------------------------
      % Check with NaN's at Start/End

      k = ( sum(c,3) == s2 );
      k = k(:);

      if ~any(k) & isn

          m = ( m01(ones(1,nk),[2 2 1 1],:) & n01(ok,[1 2 2 1],:) );

          c = ( c | m );

          d(find(m)) = 10*s2*acc;

          k = ( sum(c,3) == s2 );
          k = k(:);

      end

      %-------------------------------------------------------------
      % Find connected Segment

      if ~any(k)
          break
      end

      d = sum(d,3);                  % Sum of Deviations

      d = d(:);

      k = find(k);

      [d,ii] = min(d(k));

       ii =  k(ii);                  % Index in Check [ 1 .. 4*NK ]
       ik = ii - nk * floor(ii/nk);  % Index in kk
       ik = ik + nk * ( ik == 0 );
       k1 = ok(ik);                  % Index in Segments

      %-------------------------------------------------------------
      % Check Connection-Type
      % End-Start / End-End / Start-End / Start-Start

       ii = ceil(ii/nk);          % [ 1 .. 4 ]
 
       flp = ( mod(ii,2) == 0 );  % Flip   Segment
       app = ( ii <= 2 );         % Append Segment

       app = 1 + app;

       cc(app) = cc(app) - 3 + 2*app;

       zz(k1) = cc(app) + i * flp;
       nr(k1) = z;

       flp = 1.5 - ( 2*flp - 1 ) * ( app - 1.5 );

       y01(:,app,:) = x01(k1,flp,:);

       if isn
          m01(:,app,:) = n01(k1,flp,:);
       end


   end

end

%*****************************************************

y =  cell(z,1);
l = zeros(z,1);

nn = NaN * zeros(1,s2);

for ii = 1 : z

    ok = ( nr == ii );
    nk = sum(ok);

    ok = find(ok);

    [cc,si] = sort(real(zz(ok)));

    ok = ok(si);

    ll = i1(ok) - i0(ok) + 1;
    lg = sum(ll);

    iy = cumsum(cat(1,0,ll));

    yy = zeros(lg,s2);
    
    for ik = 1 : nk
    
        jj = ( i0(ok(ik)) : i1(ok(ik)) );
        
        if imag(zz(ok(ik)))
           jj = jj(ll(ik):-1:1);
        end

        yy(iy(ik)+(1:ll(ik)),:) = x(jj,:);

    end

    if mode
       if ~all( abs(yy(1,:)-yy(lg,:)) <= acc )
           yy = cat( 1 , yy , nn );
       end
    end

    y{ii} = yy;
    l(ii) = lg;

end

