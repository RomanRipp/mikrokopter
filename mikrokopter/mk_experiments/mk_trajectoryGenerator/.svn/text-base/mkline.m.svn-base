function mkline(firstPoint, secondPoint, minAlt, maxAlt, waypointsCount)

lat = transpose(ones(1,waypointsCount)*firstPoint(1));
lat(2:2:waypointsCount) = secondPoint(1);

lon = transpose(ones(1,waypointsCount)*firstPoint(2));
lon(2:2:waypointsCount) = secondPoint(2);

alt = transpose(ones(1,waypointsCount)*minAlt);
delta = 5;
if (minAlt ~= maxAlt)
    alt(2:2:waypointsCount) = maxAlt;
    for index = 1:2:waypointsCount
        %index
        alt(index) = alt(index) + (index-1) * delta/2;
        alt(index+1) = alt(index+1) + (index-1) * delta/2;
        %alt(index)
    end
end

%plot
figure;
scatter3(lat,lon,alt);
hold on;
plot3(lat,lon,alt);
hold off;

waypoints = [lat lon alt];
hspeed = [2 4 6];
vspeed = [2 4 6]; %this might be too much!!!
if (minAlt == maxAlt)
    writeWPFile(waypoints,-1,'hline',hspeed(1),vspeed(1));
    writeWPFile(waypoints,-1,'hline',hspeed(2),vspeed(2));
    writeWPFile(waypoints,-1,'hline',hspeed(3),vspeed(3));
else
    writeWPFile(waypoints,-1,'vline',hspeed(1),vspeed(1));
    writeWPFile(waypoints,-1,'vline',hspeed(2),vspeed(2));
    writeWPFile(waypoints,-1,'vline',hspeed(3),vspeed(3));
end
