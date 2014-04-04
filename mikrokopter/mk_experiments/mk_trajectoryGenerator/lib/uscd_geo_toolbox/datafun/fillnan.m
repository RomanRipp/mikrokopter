function  z = fillnan(z,ww,pw,sm,nz);

% FILLNAN Interpolate over NaN's in Matrice
%
% Z = FILLNAN( Z , Weight , Potenz , Smooth , ZeroInt )
%
% Weight = [ WX  WY ]
% Potenz = [ PX  PY ]
% Smooth = [ NX  NY ]
%
% ZeroInt =  BadInt + i * OkInt for overlap orig + smooth
%
% Interpolates by Column, Row, Diagonals over NaN's
%
% Weight by  1 / ( (Distance*W).^Potenz)
%
 
Nin = nargin;

if Nin < 2, ww = []; end
if Nin < 3, pw = []; end
if Nin < 4, sm = []; end
if Nin < 5, nz = []; end

%**********************************************
% Get Parameter

if prod(size(ww)) < 2
   wx = 1;
   wy = 1;
else
   wx = ww(1);
   wy = ww(2);
end

if isempty(pw)
   px = 1;
   py = 1;
elseif prod(size(pw)) == 1
   px = pw;
   py = pw;
else
   px = pw(1);
   py = pw(2);
end

if isempty(sm)
   sm = [ 1  1 ];
elseif prod(size(sm)) == 1
   sm = sm * [ wy/wx  1 ; 1  wx/wy ];
   sm = max(sm,[],1);
end

sm = 2 * floor(sm/2) + 1;

if isempty(nz)
   nz = ceil( min(sm) / 2 ) * any( sm > 1 );
end

%**********************************************

wd = sqrt( wx.^2 + wy.^2 );
pd = ( px/wx + py/wy ) / ( 1/wx + 1/wy );

n = ( 1 : 100 )';

tw = ( 1 ./ ( (n*wx).^px ) + 1 ./ ( (n*wy).^py ) ) / 2;

td =   1 ./ ( (n*wd).^pd );

fd = polyfit( n(1:20) , tw(1:20)./td(1:20) , 1 );

%%% figure, plot(n,tw./(td.*polyval(fd,n)));

%**********************************************

sc = [ 50 5+i ];
md = 'linear';

%**********************************************

s =  size(z);

n = isnan(z);

if ~any(n(:))
    return
end

z(find(n)) = 0;

w = zeros(s);

x = ( 1 : s(2) );
y = ( 1 : s(1) )';

%**********************************************************
% Run along Dimensions

cl = loopdot(sc,s(2),'Interpolate along Y');

for ii = 1 : s(2)
    nn = n(:,ii);
    if ~all(nn) & any(nn)
        ww = zerosize(~nn,s(1));
        ok = find(~nn);
        nn = find( nn);
        ww = 1 ./ ( (ww(nn)*wy).^py );
        z(nn,ii) = z(nn,ii) + ww .* interp1(y(ok),z(ok,ii),y(nn),md);
        w(nn,ii) = w(nn,ii) + ww;
    end
    loopdot(sc,s(2),ii,0,cl); 
end


cl = loopdot(sc,s(1),'Interpolate along X');

for ii = 1 : s(1)
    nn = n(ii,:);
    if ~all(nn) & any(nn)
        ww = zerosize(~nn,s(2));
        ok = find(~nn);
        nn = find( nn);
        ww = 1 ./ ( (ww(nn)*wx).^px );
        z(ii,nn) = z(ii,nn) + ww .* interp1(x(ok),z(ii,ok),x(nn),md);
        w(ii,nn) = w(ii,nn) + ww;
    end
    loopdot(sc,s(1),ii,0,cl); 
end

%**********************************************************
% Run along Diagonals

id = [ -s(1)+3   s(2)-3 ];

%---------------------------------------------------
if id(1) <= id(2)
%---------------------------------------------------

zz = s(1) * ( x - 1 );
zz = y * ones(1,s(2)) + ones(s(1),1) * zz;

nd = id(2) - id(1) + 1;

id = ( id(1) : id(2) );

 %------------------------------------------
 for kk = [ 1  2 ]
 %------------------------------------------

   if kk == 2
      zz = zz(:,s(2):-1:1);
   end

   cl = loopdot(sc,nd,sprintf('Interpolate along Diag%2.2d',kk));

   for ii = 1 : nd

       jj = diag(zz,id(ii));
       nn = n(jj);
       dn = size(jj,1);

       ok = ~nn;

       if ~all(nn) & ( sum(ok) > 1 )
           ww = zerosize(ok,dn);
           ok = find(ok);
           nn = find(nn);
           kk = jj(ok);
           mm = jj(nn);
           ww = polyval(fd,ww(nn)) ./ ( (ww(nn)*wd).^pd );
           dd = ( 1 : dn )';
           z(mm) = z(mm) + ww .* interp1(dd(ok),z(kk),dd(nn),md);
           w(mm) = w(mm) + ww;
       end

       loopdot(sc,nd,ii,0,cl); 

   end

 %------------------------------------------
 end
 %------------------------------------------

%---------------------------------------------------
end
%---------------------------------------------------

%**********************************************************

nz = [ real(nz) imag(nz) ];

if any( nz > 0 ) & any( sm > 1 )

   if nz(1) > 0

      fprintf(1,'ZEROSIZE( [ %.0f by %.0f ] , +%2.0f )',s,nz(1));

      zz = zerosize( ~n , nz(1) );

      fprintf(1,'\n');

   else
 
      zz = n;

   end

   if nz(2) > 0

      fprintf(1,'ZEROSIZE( [ %.0f by %.0f ] , -%2.0f )',s,nz(2));

      zz = zz + min( 0 , 1 - zerosize( n , nz(2) ) );

      fprintf(1,'\n');

   end

   zz = ( zz + nz(2) ) / ( nz(1) + nz(2) + 1 );

   zz = ( 1 - cos( 2*pi * zz / 2 ) ) / 2;

end

%**********************************************************

n = find(n);

w = w(n);

z(n(find(w==0))) = NaN;

z(n) = z(n) ./ ( w + ( w == 0 ) );

if all( sm <= 1 )
   return
end

%**********************************************************

sm = sm([2 1]);

fprintf(1,'MEANIND2( [ %.0f by %.0f ] , [ %2.0f %2.0f ] )',s,sm);

z = ( 1 - zz ) .* z + zz .* meanind2(z,sm,'cos');

fprintf(1,'\n');
