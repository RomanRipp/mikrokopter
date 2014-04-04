function [x,y,z] = obs_proj(c,x,y,z,R,scl)

% OBS_PROJ  Projection of spherical Coordinates to cartesian Coordinates
%
%  [X,Y,Z] = OBS_PROJ( CENTER , LON , LAT , <ALT> )
%
%  CENTER = [ CenterLongitude  CenterLatitude  <CenterAltitude> <Rotation> ]
%
%  The Rotation gives the Angle of the Y-Axis clockwise to North
%
%  default:  CenterAltitude = 0
%                  Rotation = 0
%                       ALT = 0
%
% The Direction of the Cartesian System are:
%
%    X-Axis: East  + Rotation
%    Y-Axis: North + Rotation
%    Z-Axis: Zenit ==  [ 0  0  0 ] --> CENTER
%
% Units of X, Y, Z and ALT: meter
% Units of LON , LAT: degree
%
%-----------------------------------------------------------------------
% Optional:
%
%    ... = OBS_PROJ( ... , Radius , Scale )
%
%  use Radius in Kilometer and  Projects the Data
%   by Distance of Observer to New Spheroid Radius = Radius * Scale
%
%  defaults: Radius = 6371;  % Units: [km]
%            Scale  = 1;     % Flat Surface
%
%   Radius == NaN  ==>  use DefaultValue
%   Scale  == NaN  ==>  total Reduction of Curvation of the Spheroid
%                       (flat surface)
%
%   ... = OBS_PROJ( LON , LAT , ALT , ...)
%
%   use: CENTER = [ 0  90  -Radius  0 ]
%
%
%-----------------------------------------------------------------------
% Give an Imaginary EarthRadius for reverse Projection.
%
%  [ LON , LAT , ALT ] = OBS_PROJ( CENTER , X , Y , Z , i*Radius , Scale )
%
% Note:  Values for LON will returned in [ -180 .. 180 ],
%                   NaN-Values for LAT at [ -90  90 ]
%
%-----------------------------------------------------------------------
%
% see also: SPH_PROJ, DST_PROJ
%
%-----------------------------------------------------------------------
% Equivalent to SPH_PROJ
%
%  [ X , Y , Z ] = OBS_PROJ( [ CLON CLAT -1 ] , LON , LAT , 0 , 1e-3 , 1 )
%  [ X , Y , Z ] = SPH_PROJ( [ CLON CLAT    ] , LON , LAT )
%
%-----------------------------------------------------------------------
% Equivalent to DST_PROJ
%
%  [ X , Y ] = OBS_PROJ( [ CLON CLAT ] , LON , LAT , 0 , 180/pi*1e-3 , NaN )
%  [ X , Y ] = DST_PROJ( [ CLON CLAT ] , LON , LAT )
%


Nin = nargin;

if Nin < 2
   error('Inputs [ CENTER LON LAT ] or [ LON LAT ] are missing.');
end
 
if Nin < 4
   z = [];
end

if Nin < 5
   R = NaN;
end

if Nin < 6
   scl = 1;
end

%*******************************************************
% Check for Center

sc = size(c);  pc = prod(sc);

if pc == 2
   c  = cat( 2 , c(:)' , 0 );
   pc = 3;
end

if pc == 3
   c  = cat( 2 , c(:)' , 0 );
   pc = 4;
end

if ~( pc == 4 );

    if Nin == 6
       error('Input for CENTER must have min. 2 Elements: [ Lon  Lat  <Alt>  <ROT> ].');
    end

    if Nin > 4
       scl = R;
    end
    if Nin > 3
       R = z;
    end
    if Nin > 2
       z = y;
    end

    y = x;
    x = c;
    c = [];

elseif Nin == 2 

    if ~isequal(sc,size(x))
       error('Input LAT is missing.')
    end

    y = x;
    x = c;
    c = [];

end

if isempty(z)
   z = 0;
end

is_inv = ~( imag(R) == 0 );

R = abs(R);

if isnan(R)
   R = 6371;
end

R  = R  * 1e3;  % [km]  -->  [m]

if isempty(c)
   c = [ 0  90  -R  0 ];
end

%*******************************************************
% Check Size of X, Y and Z

sx = size(x);  px = prod(sx);  vx = ( px == max(sx) );
sy = size(y);  py = prod(sy);  vy = ( py == max(sy) );
sz = size(z);  pz = prod(sz);

if ~isequal(sx,sy)

    if     vx & vy
       x = ones(py,1) * x(:)';
       y = y(:) * ones(1,px);
    elseif vx & ( px == sy(2) ) & ( py == sy(1)*sy(2) )
       x = ones(sy(1),1) * x(:)';
    elseif vy & ( py == sx(1) ) & ( px == sx(1)*sx(2) )
       y = y(:) * ones(1,sx(2));
    else
       error('Size of X and Y must be agree.');
    end

    sx = size(x);
    sy = size(y);

end

if pz == 1
   z = z*ones(sx);
elseif ~isequal(sx,sy,sz)
   error('Size of X, Y and Z must be agree.');
end

%*******************************************************

p180 = pi/180;

x = x(:); 
y = y(:);
z = z(:);


%-------------------------------------
% Basis CoordinateSystem

% ex0 = [ 1 ; 0 ; 0 ];
% ey0 = [ 0 ; 1 ; 0 ];
% ez0 = [ 0 ; 0 ; 1 ];


%-------------------------------------
% New CoordinateAxes in Basis System

z0 = c(3);

c = c([1 2 4]) * p180;

c(3) = -c(3);  % ClockWise

sc = sin(c);
cc = cos(c);

ey = [  -sc(1)*sc(2)
         cc(1)*sc(2)
         cc(2)       ];


ez = [   sc(1)*cc(2)
        -cc(1)*cc(2)
         sc(2)       ];

ex = cross(ey,ez);

%-------------------------------------
% Rotation

ry = [  -sc(3)
         cc(3) 
          0    ];

rx = [  cc(3)
        sc(3) 
         0    ];

rz = [ 0 ; 0 ; 1 ];

%--------------------------------------------
% TransformationMatrice

% e0 = [  ex0(:)  ey0(:)  ez0(:)  ];
% e1 = [  ex1(:)  ey1(:)  ez1(:)  ];
% T = e0  * inv( e1 );  % inv(e1) = e1';  inv(T)  = e1

  T = [ ex ey ez ] * [ rx ry rz ];

%********************************************************
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%  Length > 100000  ==>  Use parts off Length == 100000  

nn = prod(sx);
mm = 1e5;

m = nn + ( mm - nn ) * ( nn > mm );

n = ceil( nn / m );

for ii = 1 : n

 ende = m*ii + ( nn - m*ii ) * ( nn < m*ii );

 ind = ( m*(ii-1)+1 :  ende );

 if is_inv

   [x(ind),y(ind),z(ind)] = invproj(x(ind),y(ind),z(ind),T,R,scl,z0);

    x(ind) = x(ind) / p180;  % [rad] --> [deg]
    y(ind) = y(ind) / p180;  % [rad] --> [deg]

 else

   [x(ind),y(ind),z(ind)] = doproj(x(ind),y(ind),z(ind),T,R,scl,z0);

 end

end


%--------------------------------------------
% OutPut

x = reshape(x,sx);
y = reshape(y,sx);
z = reshape(z,sx);

if ~is_inv
    return
end

%--------------------------------------------
% Check for NaN's in equal X-Grid

x = chk_nan(x,1);
y = chk_nan(y,2);

%*******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function x = chk_nan(x,d);

s = size(x);
n = prod(s);

% Check for Grid
if ~( ( s(1)*s(2) == n ) & ( s(d) > 1 ) );
    return
end

nn = isnan(x);
sn = sum(nn,3-d);
   
if ~( any(nn(:)) & any( sn == 0 ) )
    return
end

nn = diff( round( x * 1e10 ) / 1e10 , 1 , d );
      
if ~all( isnan(nn(:)) | ( nn(:) == 0 ) )
    return
end

ii = find( sn == 0 );

if d == 1    
   x = round( x(ii(1),:) * 1e10 ) / 1e10;
   x = x(ones(s(d),1),:); 
else
   x = round( x(:,ii(1)) * 1e10 ) / 1e10;
   x = x(:,ones(1,s(d))); 
end

%*******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [x,y,z] = doproj(x,y,z,T,R,sc,z0)

% function [x,y,z] = doproj(x,y,z,T,R,sc,z0)
%
% Input: x,y,z  ...  [ N by 1 ];  % Column
% Output x,y,z  ...  [ N by 1 ];  % Column
%

p180 = pi/180;

% Cartesic Coordinates of Grid in Basis System

xyz = ~cat( 2 , ( round((x+00)/180) == (x+00)/180 ) , ...
                ( round((x+90)/180) == (x+90)/180 ) , ...
                ( round((y+90)/180) == (y+90)/180 )        );

xyz(:,1) = ( xyz(:,1) & xyz(:,3) );
xyz(:,2) = ( xyz(:,2) & xyz(:,3) );
xyz(:,3) = ~( round((y+00)/180) == (y+00)/180 );

xyz = double(xyz);

xyz(:,1) =  sin(x*p180) .* cos(y*p180) .* xyz(:,1);
xyz(:,2) = -cos(x*p180) .* cos(y*p180) .* xyz(:,2);
xyz(:,3) =  sin(y*p180) .* xyz(:,3);


%--------------------------------------------
% Projection 

xyz  = ( xyz * T );   %  == (e0*inv(e1)) * xyz'

%--------------------------------------------

is_nan = find(isnan(xyz));

xyz = min(max(xyz,-1),1);

xyz(is_nan) = NaN;

%--------------------------------------------
% Transform: 1 --> sc

if ~isequal(sc,1)

  %--------------------------------------------
  % TRUE for ( x == 0 )  &  ( y == 0 )

  zero = ( abs(xyz(:,[1 2])) < 1e3*eps );

  zero = ( sum( double(zero) , 2 ) == 2 );
 
  zero = double(zero);

  zero( find( zero & ( xyz(:,3) < 0 ) ) ) = NaN;

  %--------------------------------------------
  if isnan(sc)

    % acos(xyz(:,3)) ==  ArcLength

    k = acos(xyz(:,3)).^2;

    xyz(:,3) = 1;

    sc = 1;

  %--------------------------------------------
  else

    % (0,0,1) * (xyz1,xyz2,xyz3 ) = xyz3 = ( R + z ) * cos(phi)
    %
    % |xyz| = R + z
    %
 
    xyz(:,3) = cos( 1/sc * acos( xyz(:,3) ) );
 
    k = ( 1 - xyz(:,3).^2 );

  end

  %--------------------------------------------
  % Project X & Y

    k = k ./ ( xyz(:,1).^2 + xyz(:,2).^2 + zero );

    k = sqrt( k .* ( 1 - zero ) );

    xyz(:,[1 2]) = k(:,[1 1]) .* xyz(:,[1 2]);

end

%--------------------------------------------
% Scale with Distance from Center: ( sc*R + z )

xyz = ( sc*R + z(:,ones(1,3)) ) .* xyz;
  
%--------------------------------------------

x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3) - ( sc*R + z0 );  % ( Z == ZERO ) at  ( Z == z0 )


%*******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [x,y,z] = invproj(x,y,z,T,R,sc,z0)

% function [x,y,z] = invproj(x,y,z,T,R,sc,z0)
%
% Input: x,y,z  ...  [ N by 1 ];  % Column
% Output x,y,z  ...  [ N by 1 ];  % Column
%

xyz = cat( 2 , x , y , z );

xyz = xyz .* ( abs(xyz) > 1e3*eps );

sc_nan = isnan(sc);

if sc_nan
   sc = 1;
end

%--------------------------------------------
% Z == ZERO at Z == z0

xyz(:,3) = xyz(:,3) + ( sc*R + z0 );

if sc_nan
   z = xyz(:,3);
else
   z = sqrt( sum( xyz.^2 , 2 ) );
end

%--------------------------------------------
% Scale from Distance from Center to 1: ( sc*R + z ) --> 1

xyz = xyz ./ z(:,ones(1,3));  % normalized

  z = z - sc * R;

%--------------------------------------------
% reverse Transformation: sc --> 1

if ~isequal(sc,1) | sc_nan

  %--------------------------------------------
  % TRUE for ( x == 0 )  &  ( y == 0 )

  zero = ( abs(xyz(:,[1 2])) < 1e3*eps );

  zero = ( sum( double(zero) , 2 ) == 2 );
 
  zero = double(zero);


  %--------------------------------------------
  if sc_nan

   xyz(:,3) = cos( sqrt( xyz(:,1).^2 + xyz(:,2).^2 ) );

  %--------------------------------------------
  else

    xyz(:,3) = cos( sc * acos( xyz(:,3) ) );

  end

  %--------------------------------------------
  % Project X & Y

    k = ( 1 - xyz(:,3).^2 ) ./ ( xyz(:,1).^2 + xyz(:,2).^2 + zero );

    k = sqrt( k .* ( 1 - zero ) );

    xyz(:,[1 2]) = k(:,[1 1]) .* xyz(:,[1 2]);

end

%--------------------------------------------
% inverse Projection 

xyz  = xyz * T';  % == ( inv(e0*inv(e1)) * xyz ) == ( xyz *  e1' );

%--------------------------------------------

is_nan = find(isnan(xyz));

xyz = min(max(xyz,-1),1);

xyz(is_nan) = NaN;

%--------------------------------------------

y = asin( xyz(:,3) );

xyz(:,3) = cos(y);                 % cos(lat)

xyz( find( abs(xyz(:,3)) < 1e8*eps ) , 3 ) = NaN;

xyz(:,1) =  xyz(:,1) ./ xyz(:,3);  % sin(lon)
xyz(:,2) = -xyz(:,2) ./ xyz(:,3);  % cos(lon)

x = atan2( xyz(:,1) , xyz(:,2) );

