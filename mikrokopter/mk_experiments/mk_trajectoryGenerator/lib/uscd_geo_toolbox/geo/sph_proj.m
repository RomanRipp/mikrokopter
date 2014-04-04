function [x,y,z] = sph_proj(c,x,y,z)

% SPH_PROJ  Projection of spherical Coordinates to cartesian Coordinates of Sphere
%
%  [ X , Y , Z ] = SPH_PROJ( CENTER , LON , LAT )
%  [ LON , LAT ] = SPH_PROJ( CENTER , X , Y , Z );
%
%  CENTER = [ CenterLongitude  CenterLatitude <Rotation> ]
%
%  The Rotation gives the Angle of the Y-Axis clockwise to North
%      default: Rotation = 0
%
%  [ X , Y , Z ] = SPH_PROJ( LON , LAT )
%  [ LON , LAT ] = SPH_PROJ( X , Y , Z );
%
%   use: CENTER = [ 0  0  0 ]
%
% The Direction of the Cartesian System are:
%
%    X-Axis: East  + Rotation
%    Y-Axis: North + Rotation
%    Z-Axis: Zenit ==  [ 0  0  0 ] --> CENTER
%
% Units of LON, LAT are in degree.
% Units of X, Y and Z are normalized: [ -1 .. 1 ]
%
% Note:  Values for LON will returned in [ -180 .. 180 ],
%                   NaN-Values for LAT at [ -90  90 ]
%
%-----------------------------------------------------------------------
%
% see also: OBS_PROJ, DST_PROJ
%
%-----------------------------------------------------------------------
% Equivalent to OBS_PROJ
%
%  [ X , Y , Z ] = OBS_PROJ( [ CLON CLAT -1 ] , LON , LAT , 0 , 1e-3 , 1 )
%  [ X , Y , Z ] = SPH_PROJ( [ CLON CLAT    ] , LON , LAT )
%


Nin = nargin;

if Nin < 2
   error('Inputs [ CENTER LON LAT ] or [ LON LAT ] are missing.');
end
 
%*******************************************************
% Check for Center

sc = size(c);  pc = prod(sc);

if pc == 2
   c  = cat( 2 , c(:)' , 0 );
   pc = 3;
end

if ~( pc == 3 );

    if Nin == 4
       error('Input for CENTER must have min. 2 Elements: [ Lon Lat <ROT> ].');
    end

    if Nin > 2
         z = y;
       Nin = 4;
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

is_inv = ( Nin == 4 );

if isempty(c)
   c = [ 0  0  0 ];
end

%*******************************************************
% Check Size of Inputs

sx = size(x);  px = prod(sx);  vx = ( px == max(sx) );
sy = size(y);  py = prod(sy);  vy = ( py == max(sy) );

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

if is_inv
   sz = size(z);  pz = prod(sz);
   if pz == 1
       z = z * ones(sx);
   elseif ~isequal(sx,sy,sz)
       error('Size of X, Y and Z must be agree.');
   end
end
     
%*******************************************************

p180 = pi/180;

x = x(:); 
y = y(:);

if is_inv
   z = z(:);
end

%-------------------------------------
% Basis CoordinateSystem

% ex0 = [ 1 ; 0 ; 0 ];
% ey0 = [ 0 ; 1 ; 0 ];
% ez0 = [ 0 ; 0 ; 1 ];


%-------------------------------------
% New CoordinateAxes in Basis System

c = c * p180;

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

if ~is_inv
    z = zeros(nn,1);
end

for ii = 1 : n

 ende = m*ii + ( nn - m*ii ) * ( nn < m*ii );

 ind = ( m*(ii-1)+1 :  ende );

 if is_inv

   [x(ind),y(ind)] = invproj(x(ind),y(ind),z(ind),T);

    x(ind) = x(ind) / p180;  % [rad] --> [deg]
    y(ind) = y(ind) / p180;  % [rad] --> [deg]

 else

   [x(ind),y(ind),z(ind)] = doproj(x(ind),y(ind),T);

 end

end


%--------------------------------------------
% OutPut

x = reshape(x,sx);
y = reshape(y,sx);

if ~is_inv
    z = reshape(z,sx);
    return
else
    z = [];
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

function [x,y,z] = doproj(x,y,T)

% function [x,y,z] = doproj(x,y,T)
%
% Input: x,y    ...  [ N by 1 ];  % Column
% Output x,y,z  ...  [ N by 1 ];  % Column
%
% T = [ ex ey ez ]
%

p180 = pi/180;

%--------------------------------------------
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

is_nan = find(isnan(xyz));

xyz = min(max(xyz,-1),1);

xyz(is_nan) = NaN;

%--------------------------------------------

x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);

%*******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [x,y] = invproj(x,y,z,T)

% function [x,y] = invproj(x,y,T)
%
% Input: x,y,z  ...  [ N by 1 ];  % Column
% Output x,y    ...  [ N by 1 ];  % Column
%

xyz = cat( 2 , x , y , z );

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
