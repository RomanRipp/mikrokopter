function [d,b] = horizon(h1,h2,r)

% HORIZON   Calculate Visibility of Object over the Horizon
%
% [ D , B ] = HORIZON( H1 , H2 , [Radius] )
%
% H1 / H2 Height of Point1 / Point2 in Meter
%
% Radius  Radius of Sphere in KiloMeter
%
% D Tangent Distance between Point1 and Point2 
%   
% B Geodetic Distance between Point1 and Point2
%   on Sphere
%
% D and B in KiloMeter
%
 
d = [];
b = [];

Nin = nargin;

if Nin == 0
   error('Input Height is missing.');
end

if Nin < 2
   h2 = 0;
end

if Nin < 3
   r = [];
end

if  ( isempty(h1) | isempty(h2) )
    s1 = size(h1);
    s2 = size(h2);
    if ~( isequal(s1,s2) | ( prod(s1) == 1 ) | ( prod(s2) == 1 ) )
        error('Heights must be single or of same size.');
    end
end

if isempty(r)
   r = 6871;
elseif ~( prod(size(r)) == 1 )
   error('Radius must be a single Element.')
end

if ( isempty(h1) | isempty(h2) )
   return
end

%*******************************************************

h1 = h1/1000 + r;
h2 = h2/1000 + r;

r2 = r^2;

s1 = sqrt( h1.^2 - r2 );
s2 = sqrt( h2.^2 - r2 );

d = s1 + s2;

b = r * ( asin(s1./h1) + asin(s2./h2) );
