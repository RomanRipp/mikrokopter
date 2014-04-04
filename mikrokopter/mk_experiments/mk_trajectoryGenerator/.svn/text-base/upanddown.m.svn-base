function upanddown(centerCoordinates, minAlt, maxAlt, waypointsCount)

delta = 5;
lat = transpose(ones(1,waypointsCount)*centerCoordinates(1));
lon = transpose(ones(1,waypointsCount)*centerCoordinates(2));
alt = transpose(ones(1,waypointsCount)*minAlt);
alt(2:2:waypointsCount) = maxAlt;

for index = 1:2:waypointsCount
    %index
    alt(index) = alt(index) + (index-1) * delta/2;
    alt(index+1) = alt(index+1) + (index-1) * delta/2;
    %alt(index)
end

alt

%plot
figure;
scatter3(lat,lon,alt);
hold on;
plot3(lat,lon,alt);
hold off;

waypoints = [lat lon alt];
hspeed = [2 4 6];
vspeed = [2 4 6]; %this might be too much!!!
writeWPFile(waypoints, -1, 'upanddown',hspeed(1),vspeed(1));
writeWPFile(waypoints, -1, 'upanddown',hspeed(2),vspeed(2));
writeWPFile(waypoints, -1, 'upanddown',hspeed(3),vspeed(3));
