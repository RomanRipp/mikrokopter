function y = dist(lat,lon)
%DIST Distance between geogralatcal points.
%  Y = DIST(LAT,LON)

transpose = 0;
lat = lat*pi/180;
lon = lon*pi/180;
[m,n] = size(lat);
if m == 1
  lat = lat(:);
  lon = lon(:);
  m = n;
  transpose = 1;
end

y = abs(sin(lat(1:m-1,:)).*sin(lat(2:m,:)) ...
    + cos(lat(1:m-1,:)).*cos(lat(2:m,:)).*cos(lon(2:m,:)-lon(1:m-1,:)));
i = find(y > 1);
if ~isempty(i)
  y(i) = ones(i);
end
y = 180/pi*atan(sqrt((1 - y)./(1 + y)))*222240;
if transpose == 1
  y = y';
end
