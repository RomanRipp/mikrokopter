function   v = interp0(v,ok,Name,P,schwelle);

% INTERP0   Interpolate linear over gabs 
%
% V = interp0(V,Ok,Name,[Period],[Level]);
%   
% Interpolates over gabs in V signed by ( Ok == 0 )
%
% In case of inputs Period and Level the Data are checked
% for PeriodOverrunnings of Level*Period before and will
% made cyclic (like CompassData)!
%
% see also: CYCLIC
%

if nargin < 3
 Name = '';
end

nl = char(10);

l1 = NaN;
l2 = 50;

NperRow = 500;  % Maximum


i0 = find( diff([1;ok;1]) < 0 ) - 1;
i1 = find( diff([1;ok;1]) > 0 ) ;

   kk  = find( ( i0 < 1 ) | ( i1 > size(v,1) ) );
i0(kk) = [];
i1(kk) = [];

if isempty(i0)
  return
end


if nargin == 5

   if nargin < 6
      schwelle = 0.5;
   end

   is_ok = find(ok);

    w = v(is_ok);

   dw = [ zeros(1,size(w,2)) ; diff(w,1,1) ];

  add = ( abs(dw) > schwelle*P ) .* (-sign(dw));

    w = w + cumsum(P*add,1);
  
   v(is_ok) = w;

end

N = size(i0,1);


l1 = ceil(min(NperRow,N)/l2);

fprintf([ 'Interpolate '   Name  ' %.0f times' ] , N );

if l1 > 1
  fprintf([ ', %.0f per Point' ] , l1 );
end
   
if N > l1*l2
   fprintf([', %.0f per Row' ] , l1*l2 );
end

fprintf(nl)

str = [ '.' nl ];

o2 = ones(1,size(v,2));
  
for ii = 1 : N

    fprintf( str( 1 : 0+(mod(ii,l1)==0)+(mod(ii,l1*l2)==0) ) );

    xi = ( i0(ii)+1 : i1(ii)-1 )';
    xi = ( xi - i0(ii) ) / ( i1(ii) - i0(ii) );  % [ 0 .. 1 ]

    o1 = ones(size(xi));

    dv = ( v(i1(ii),:) - v(i0(ii),:) );

    v(i0(ii)+1:i1(ii)-1,:) = v(o1*i0(ii),:) + dv(o1,:) .* xi(:,o2);

% v(i0(ii)+1:i1(ii)-1,:) =  ...
%    interp1(     [ i0(ii)     i1(ii)   ]        , ...
%              v( [ i0(ii)     i1(ii)   ] , :  ) , ...
%                 ( i0(ii)+1 : i1(ii)-1 )'              );

end

fprintf(nl)


if nargin == 5

  v = v - P*floor(v/P);

end
