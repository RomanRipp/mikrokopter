% gravitational constant
% function gr = g(lat,p)
%
% T.Kanzow

function gr = g(lat,p)

 p     = p / 10.;

 x     = sin(lat/57.29578).^2;

 gr    = 9.780318*(1.0+(5.2788e-3+2.36e-5*x).*x)+1.092e-5*p;
