clear all;

centerCoordinates = [44.9704076, -93.2352322];
edgeCoordinates = [44.9703469 , -93.2349586];
minAlt = 10;
maxAlt = 20;
waypointsCount = 10;

% Do the spiral thing:
spiral(centerCoordinates, edgeCoordinates, minAlt, maxAlt, waypointsCount);

% Do the vertical thing:
upanddown(centerCoordinates, minAlt, maxAlt, waypointsCount);

% Do the horisontal line thing:
firstPoint = [44.9704076, -93.2352322];
secondPoint = [44.9703469 , -93.2349586];
maxAlt = 10;
mkline(firstPoint, secondPoint, minAlt, maxAlt, waypointsCount);

% Do the accend decscend line thing:
maxAlt = 20;
mkline(firstPoint, secondPoint, minAlt, maxAlt, waypointsCount);
