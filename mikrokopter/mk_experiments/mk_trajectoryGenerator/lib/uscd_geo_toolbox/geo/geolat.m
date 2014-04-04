function x = geolat(x,m,rev);

% GEOLAT  Transforms geodetic Latitudes, using WGS84-Ellipsoid
%
% Y = GEOLAT( LAT , Mode );
%
% The geodetic latitude LAT [deg] is the angle which a line 
%  perpendicular to the surface of the ellipsoid at the given point 
%  makes with the equatorial plane.
%
% Mode:  'a'  Authalic
%        'c'  Conformal
%        'g'  Geocentric
%        'i'  Isometric
%        'p'  Parametric
%        'r'  Rectifying
%
%------------------------------------------------------------------
%
% Use a NonZero 3. Input for inverse Transformation:
%
%   LAT = GEOLAT( Y , Mode , 1 );
%
%------------------------------------------------------------------
%
% See also: GEORAD
%
%------------------------------------------------------------------
% Authalic latitude
%  The authalic latitude is used to map an ellipsoid to a sphere 
%  in such a way that the sphere has equal surface area as the 
%  ellipsoid.  Authalic latitudes are used in place of the
%  geodetic latitudes when projecting the ellipsoid using
%  an equal area projection.
%
%------------------------------------------------------------------
% Conformal latitude
%  The conformal latitude is used to map an ellipsoid conformally 
%  onto a sphere. Conformal latitudes are used in place of the 
%  geodetic latitudes when projecting the ellipsoid using
%  a conformal projection.
%
%------------------------------------------------------------------
% Geocentric latitude
%  The geocentric latitude is the angle made by a line from a point 
%  on the surface of the ellipsoid to the center of the ellipsoid
%  and the equatorial plane.
%
%------------------------------------------------------------------
% Parametric latitude
%  The parametric latitude of a point on the ellipsoid is the latitude
%  on a sphere (of radius r = semimajor axis of ellipsoid) for which
%  the parallel has the same radius as the parallel of geodetic latitude.
%
%------------------------------------------------------------------
% Rectifying latitude
%  The rectifying latitude is used to map an ellipsoid to a sphere 
%  in such a way that correct distances along meridians are preserved.  
%  Rectifying latitudes are used in place of the geodetic latitudes when 
%  projecting the ellispoid using an equal distant projection. 
%
%------------------------------------------------------------------
% Isometric latitude
%  The isometric latitude is a nonlinear function of the geodetic
%  latitude.  It is directly proportional to the spacing of parallels
%  of geodetic latitude from the Equator on the ellipsoidal Mercator
%  projection.
%
%


Nin = nargin;

if Nin < 2
   error('Not enough Input Arguments.');
end

if isempty(x)
   return
end

if Nin < 3
   rev = [];
end

if isempty(rev)
   rev = 0;
end

rev = ~isequal(rev,0);

%******************************************************************
% WGS84

a = 6378137.;
b = 6356752.314;        %       A * sqrt(1-E*E);
e = sqrt(a*a - b*b)/a;  %   E = 0.081819191;

%******************************************************************

m = upper(m(1));

acc  = 1e3 * eps;  % Accuracy
p180 = pi/180;

if ~strcmp(m,'I')
    x = x * p180;
end

%******************************************************************

switch m

  %****************************************************************
  case 'A'  % Authalic
  %****************************************************************

    r = a/2 * sqrt( 2 + ( 1/e - e ) * log((1+e)/(1-e)) );

    %--------------------------------------------------------------
    if rev
    %--------------------------------------------------------------

       f1 = 1/3*e^2 + 31*e^4 / 180 + 517*e^6 /  5040;
       f2 =           23*e^4 / 360 + 251*e^6 /  3780;
       f3 =                          761*e^6 / 45360;

    %--------------------------------------------------------------
    else
    %--------------------------------------------------------------

       f1 = 1/3*e^2 + 31*e^4 / 180 + 59*e^6 /   560;
       f2 =           17*e^4 / 360 + 61*e^6 /  1260;
       f3 =                         383*e^6 / 45360;

       f1 = -f1;
       f3 = -f3;

    %--------------------------------------------------------------
    end
    %--------------------------------------------------------------

    x = x + f1*sin(2*x) + f2*sin(4*x) + f3*sin(6*x);

  %****************************************************************
  case 'C'   % Conformal
  %****************************************************************

    %--------------------------------------------------------------
    if rev
    %--------------------------------------------------------------

       f1 = e^2 / 2 + 5*e^4 / 24 +    e^6 /  12 +   13*e^8 /    360;
       f2 =           7*e^4 / 48 + 29*e^6 / 240 +  811*e^8 /  11520;
       f3 =                         7*e^6 / 120 +   81*e^8 /   1120;
       f4 =                                       4279*e^8 / 161280;

        x = x +  f1*sin(2*x) + f2*sin(4*x) + ...
                 f3*sin(6*x) + f4*sin(8*x);


    %--------------------------------------------------------------
    else
    %--------------------------------------------------------------

       f1 = 1 - e*sin(x);
       f2 = 1 + e*sin(x);
       f3 = 1 -   sin(x);
       f4 = 1 +   sin(x);

          ii  = find( f3 < acc );
       f3(ii) = 1;

       x = 2 * atan(sqrt((f4./f3) .* ((f1./f2).^e)) ) - pi/2;

       x(ii) = pi/2;

    %--------------------------------------------------------------
    end
    %--------------------------------------------------------------


  %****************************************************************
  case 'R'   % Rectifying
  %****************************************************************

        n = (a-b) / (a+b);  % (1-sqrt(1-e.^2)) / (1+sqrt(1-e.^2));

    %--------------------------------------------------------------
    if rev
    %--------------------------------------------------------------

       f1 = 3*n / 2 - 27*n^3 / 32;
       f2 =           21*n^2 / 16 -  55*n^4 /  32;
       f3 =          151*n^3 / 96;
       f4 =                        1097*n^4 / 512;

    %--------------------------------------------------------------
    else
    %--------------------------------------------------------------

        f1 = 3*n / 2 - 9*n^3 / 16;
        f2 =          15*n^2 / 16 -  15*n^4 /  32;
        f3 =          35*n^3 / 48;
        f4 =                        315*n^4 / 512;

        f1 = -f1;
        f3 = -f3;

    %--------------------------------------------------------------
    end
    %--------------------------------------------------------------


        x = x +  f1*sin(2*x) + f2*sin(4*x) + ...
                 f3*sin(6*x) + f4*sin(8*x);

  %****************************************************************
  case 'I'   % Isometric
  %****************************************************************

    %--------------------------------------------------------------
    if rev
    %--------------------------------------------------------------
 
        x = 2 * atan(exp(x)) - pi/2;

        x = geolat(x*p180,'c',1);

    %--------------------------------------------------------------
    else
    %--------------------------------------------------------------

       x = geolat(x,'c',0) * p180;

       f = ( abs(x+pi/2) < acc );
       f = ( 1 - 2*f );

       x = f.* log( abs( tan( pi/4 + f.*x/2 ) ) );

       x = x / p180;

    %--------------------------------------------------------------
    end
    %--------------------------------------------------------------

    return

  %****************************************************************
  case { 'G'  'P' }   % Geocentric | Parametric
  %****************************************************************

    ok = ( abs( ( abs(x) - pi/2 ) ) >= acc );

    f = ( 1 - e^2 ) ^ ( 1 - 0.5*strcmp(m,'P') );

    f = f .^ (1-2*rev);
 
    x = atan( f * tan(x.*ok) ) + pi/2 * sign(x).* (~ok);

  %****************************************************************
  otherwise
  %****************************************************************

    error('Invalid Mode.');

end

%******************************************************************

x = x / p180;
