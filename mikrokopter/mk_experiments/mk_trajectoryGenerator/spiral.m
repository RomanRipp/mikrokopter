function spiral(center, edge, minAlt, maxAlt, waypoints)
%altitude generator

altRange = (maxAlt - minAlt) / (waypoints-1);
alt = transpose(minAlt:altRange:maxAlt);

%Create circle:
%[latout,lonout] = scircle1(center, 25 / 6371, [], [], [], 10);
[lat,lon] = scircle2(center(1), center(2),edge(1),edge(2), [], [], waypoints);
%Visualize
figure;
scatter3(center(1), center(2), 0,'d');
hold on;
%cylinder(10);
scatter3(center(1), center(2), maxAlt,'d');
plot3(ones(length(alt),1) * center(1), ones(length(alt),1) * center(2), alt, 'color', [0 1 0], 'LineWidth' , 3);

scatter3(lat,lon,alt);
plot3(lat,lon,alt);
hold off;
%waypoints = [lon lat alt];

%write to g maps
%trajectory = 'TreeFlyOverTrajectory.kml';
%path = ge_plot3(lon,lat,alt,'altitudeMode','relativeToGround');
%ge_output(trajectory,path,'name',trajectory);

waypoints = [lat lon alt];
hspeed = [2 4 6];
vspeed = [2 4 6]; %this might be too much!!!
writeWPFile(waypoints, center, 'spiral',hspeed(1),vspeed(1));
writeWPFile(waypoints, center, 'spiral',hspeed(2),vspeed(2));
writeWPFile(waypoints, center, 'spiral',hspeed(3),vspeed(3));
end
