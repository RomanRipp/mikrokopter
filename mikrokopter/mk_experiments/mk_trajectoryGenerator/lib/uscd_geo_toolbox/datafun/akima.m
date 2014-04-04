function yy = akima(xi,yi,xx);

% AKIMA  Akima spline interpolation
%
% YI = AKIMA(X,Y,XI)
%
% performs Akima spline interpolation (a local interpolant)
% returns vector YI of interpolants at positions XI using data vectors X,Y
% input data will be sorted in ascending X: multiple values are not dealt with.
%

%P. Robbins 2/95

%put all in column vectors
xi = xi(:); yi = yi(:); xx = xx(:);
bad = isnan(yi) | isnan(xi);
if any(bad)
  yi = yi(~bad);
  xi = xi(~bad);
end

%sort xi;
[xi,I] = sort(xi);
yi = yi(I);

n = length(xi);
%extend series by two points as required by this method
xi = [xi; xi(n)+diff(xi(n-2:n-1)); xi(n)+diff(xi([n-2,n]))];
%
% "predict " values at these points assuming a qaudratic fit over last 3 pts
A = [[xi(n-2:n).^2]'; [xi(n-2:n)]' ; 1 1 1];
B = yi(n-2:n)'/A;
yi(n+1:n+2) = polyval(B,xi(n+1:n+2));

%do same at beginning of series

A = [[xi(1:3).^2]'; [xi(1:3)]' ; 1 1 1];
B = yi(1:3)'/A;
xi = [xi(1)-diff(xi([1 3])) ;xi(1)-diff(xi(2:3)); xi];
yi = [polyval(B,xi(1:2)); yi];

%compute  slope for each data pair
h = diff(xi);
m = diff(yi)./h;

%compute t for each data pair
for i = 3:n+2;
  xnum = m(i-1)*abs(m(i+1)-m(i)) + m(i)*abs(m(i-1)-m(i-2));
  denom = abs(m(i+1)-m(i)) + abs(m(i-1)-m(i-2));
  t(i) = xnum / ( denom + ( denom == 0 ) ) .* (~(denom==0));
end

%now do the interpolation
% initialize the prediction vector with nan's
yy = nan*xx;
% loop through all polynomial intervals and interpolate onto any valid
% points

for i= 3:n+1;
  xmin = xi(i); xmax = xi(i+1);  
  fs = ( ( xx >= xmin ) & ( xx <= xmax ) );
  if any(fs)
    %compute coeffecients for interpolating polynomial
    C(4) = yi(i);
    C(3) = t(i);
    C(2) = (3*(yi(i+1)-yi(i))/h(i) - 2*t(i) - t(i+1))/h(i);
    C(1) = (t(i) + t(i+1) - 2*(yi(i+1)-yi(i))/h(i))/(h(i)^2);
    yy(fs) = polyval(C,(xx(fs) - xi(i)));
  end
end

