function [h,v,f,c,x,y] = bar3d(varargin)

% BAR3D    3-Dimensional Bar-Plot of a Matrice, based to ZERO
%
%------------------------------------------------------------
% Inputs:
%
%  BAR3D(Z)
%  BAR3D(Z,C)
%  BAR3D(X,Y,Z)
%  BAR3D(X,Y,Z,C)
%
%  BAR3D( ... , Property , Value , ... )
%
%  Z = [ Nx by Ny ]
%  C = [ Nx by Ny ]
%
%  X = [ 1  by Nx ] | [    1 by Nx+1 ];  default: ( 0 : Nx )
%  Y = [ Ny by 1  ] | [ Ny+1 by    1 ];  default: ( 0 : Ny )
%
% The Data will shown as Bars based to ZERO in Z.
% If a logarithmical Z-Scale is used, the Base is ONE. 
% Use not negative Values for Z in that case. 
%
%------------------------------------------------------------
% OutPuts:
%
% [PatchHandle,Vertices,Faces,CData,X,Y] = BAR3D( ... );
%
% The CData depends on the FaceColor of the Patch:
%
% [ ... , VertexCData , ... ] = BAR3D( ... , 'FaceColor' , 'interp' );
% [ ... ,   FaceCData , ... ] = BAR3D( ... , 'FaceColor' , 'flat'   );
%
%------------------------------------------------------------
% if NO Plot is required, set the first Property to 'none'
%
% [PatchHandle,Vertices,Faces,VertexCData] = BAR3D( ... , 'none' );
%  
% The PatchHandle is empty in that case.
%
%------------------------------------------------------------
% example:
%
%  a = rand(4,5);
%  
%  a(3:4,2) = NaN;
%  a(  2,5) = NaN;
%  a(  4,4) = -a(4,4);
%  
% fig = figure;
% axe = axes( 'parent' , fig , ...
%             'view'   , [ 160 30 ]   );
%
% [h,v,f,c] = bar3d( a ,'parent'       , axe      , ...
%                       'cdatamapping' , 'scaled' , ...
%                       'facecolor'    , 'flat' , ...
%                       'edgecolor'    , [ 0  0  0 ]     );
%
% xlabel('x');
% ylabel('y');
%

Nin  = nargin;
Nout = nargout;

nl = char(10);

if Nin == 0
   error('Not enough InputArguments.')
end

%---------------------------------------
% Look for the First String-Property

Vin = varargin;

i0 = Nin+1;

for ii = 1 : Nin
    if ischar(Vin{ii})
       i0 = ii;
       break
    end
end

pind = ( i0 : Nin );  % Index for PatchInputs

if i0 <= Nin
 
  Nin = i0 - 1;

  if ( Nin == 0 )
    error('Input Z is missing.')
  end

end

FaceIn = Vin(1:Nin);  % Inputs for FACE3D
Vin    = Vin(pind);   % Inputs for PATCH

%-----------------------------------------------------------
% Check additional Inputs, using for Patch

show_patch = 1;

if ~isempty(Vin)

   if ischar(Vin{1});
      show_patch = ~isequal(lower(Vin{1}),'none');
   end

   if show_patch
   % Check Property-Value-Inputs for PATCH

     Nin = prod(size(Vin));
     ok = ( mod(Nin,2) == 0 );
     if ok
        ok = iscellstr(Vin(1:2:(Nin-1)));
        if ok
           try
             cat(2,Vin{1:2:(Nin-1)});
           catch
             ok = 0;
           end
        end
     end

     if ~ok
        error(['Following Inputs for Patch must be Property-Value-Pairs.' ...
               nl  'Properties must be Strings.' ]);
     end

   end
   % show_patch

end

%-----------------------------------------------------------
% Create Patch first

if show_patch

  try
    h = patch( Vin{:} );
  catch
    error([ 'Invalid Inputs for PATCH.'  nl  lasterr ]);
  end

else

  h = [];

end

%-----------------------------------------------------------
% Create Vertices, Face, FaceVertexCData

[v,f,c,x,y] = face3d(FaceIn{:});

if ~show_patch
   return
end

%-----------------------------------------------------------
% Set PATCH

cm = get( h , 'cdatamapping' );
fc = get( h , 'facecolor'    );

%--------------------------------------
% Set CData  with  Z == ZERO 

if strcmp( cm , 'direct' )

   c(find(isnan(c))) = 1;

else

   if strcmp( get(gca,'climmode') , 'auto' )

     cmin = min(c(:));

   else

     cmin = min( get(get(h,'parent'),'clim') );
  
   end

   c(find(isnan(c))) = cmin;

end

%--------------------------------------
% FaceVertexCData

if isequal( fc , 'flat' )
% One Color per Face

  % search for maximum Elevation of Face

  z = abs(v(:,3));
  z = z(f);

  [z,si] = sort( z , 2 );

   nf = size(f,1);

   si = (1:nf)' + nf * ( si(:,4) - 1 );
 
   c = c(f);
   c = c(si);

end


axe = get(h,'parent');

if strcmp( get(axe,'zscale') , 'log' );
   v(:,3) = max( v(:,3) ,1 );
end


set( h , 'vertices'        , v , ...
         'faces'           , f , ...
         'facevertexcdata' , c        );


axe = get(h,'parent');


if strcmp( get(axe,'xlimmode') , 'auto' )
   set( axe , 'xlim' , [ min(x)  max(x) ] );
end

if strcmp( get(axe,'ylimmode') , 'auto' )
   set( axe , 'ylim' , [ min(y)  max(y)] );
end
   


%*******************************************************************
%-------------------------------------------------------------------

function [v,f,d,x,y] = face3d(x,y,z,c)

% FACE3D    Creates Vertices and Faces for BAR3D
%
% [Vertices,Faces,CData] = FACE3D( Z )
% [Vertices,Faces,CData] = FACE3D( Z , C )
% [Vertices,Faces,CData] = FACE3D( X , Y , Z )
% [Vertices,Faces,CData] = FACE3D( X , Y , Z , C )
%
%  Z = [ Nx by Ny ]
%  C = [ Nx by Ny ]
%
%  X = [ 1  by Nx ] | [    1 by Nx+1 ];  default: ( 0 : Nx )
%  Y = [ Ny by 1  ] | [ Ny+1 by    1 ];  default: ( 0 : Ny )
%
%


Nin  = nargin;
Nout = nargout;

if Nin == 0
   error('Not enough InputArguments.')
end

v = zeros(0,3);
f = zeros(0,4);
d = zeros(0,1);

is_xy = ( Nin > 2 );

is_c = ( mod(Nin,2) == 0 );
is_d = ( Nout >= 3 );

is_c = ( is_c & is_d );

%-----------------------------------------------------------
% Make X-Y-Vectors, if not in Inputs

if ~is_xy

  z = x;
  
  if is_c
     c = y;
  end
     
  x = ( 1 : size(z,2) );
  y = ( 1 : size(z,1) );

end

%-----------------------------------------------------------
% Check X-Y-Z-C-Size

if isempty(z)
   return
end

si = size(z);

if ~( prod(si) == si(1)*si(2) )
   error('Number of Dimensions of Z must be 2.')
end

  x = x(:)';
  y = y(:);

  if size(x,2) == si(2)

     x = cat( 2 , ( x(1) - 0.5 * ( x(2) - x(1) ) )     , ...
                  ( 0.5 * ( x( 1 : (si(2)-1) ) + x( 2 : si(2) ) ) ) , ...
                  ( x(si(2)) + 0.5 * ( x(si(2)) - x(si(2)-1) ) )           );

  end

  if size(y,1) == si(1)

     y = cat( 1 , ( y(1) - 0.5 * ( y(2) - x(1) ) )     , ...
                  ( 0.5 * ( y( 1 : (si(1)-1) ) + y( 2 : si(1) ) ) ) , ...
                  ( y(si(1)) + 0.5 * ( y(si(1)) - y(si(1)-1) ) )           );

  end

  if ~isequal( [ size(y,1) size(x,2) ] , si+1 )
     error( 'Length of X and Y must match or ONE more then Size of Z.' );
  end


if is_c

   if ~isequal( si , size(c) )
      error('Size of Z and C must be agree.')
   end
   
   z( find(isnan(c)) ) = NaN;

end

%-----------------------------------------------------------

s1 = size(z,1);
s2 = size(z,2);

g1 = s1+2;
g2 = s2+2;
g3 =  4;    % 4 Vertices per Face
g4 =  3;    % 3 Coordinates: [ X Y Z ]


%---------------------------------------------------

i1 = ( 1 : s1 )' + 1;
i2 = ( 1 : s2 )  + 1;

i3 = ( 1 : g3 )  - 1;

j1 = ( 1 : ( s1 + 1 ) )';
j2 = ( 1 : ( s2 + 1 ) ) ;


o1 = ones( s1 ,  1 );
o2 = ones(  1 , s2 );

o3 = ones(  1 , g3 );

q1 = ones( s1+1 ,    1 );
q2 = ones(    1 , s2+1 );


s12 = s1 * s2;

g12 = g1 * g2;

%***************************************************
% Vertices, append a ZERO-Point arround  (Eckpunkte)
%
%  4 <------ 3
%            |
%     Z(ii)  |
%            |
%  1 ------> 2
%
%---------------------------------------------------

v = NaN * zeros( g1 , g2 , g3 , g4 );

if is_d
   d = NaN * zeros( g1 , g2 , g3 );
end


%---------------------------------------------------
% Fill Points around Z

v( i1 , i2 , 1 , 1 ) = x( o1 , i2-1 );
v( i1 , i2 , 2 , 1 ) = x( o1 , i2+0 );
v( i1 , i2 , 3 , 1 ) = x( o1 , i2+0 );
v( i1 , i2 , 4 , 1 ) = x( o1 , i2-1 );

v( i1 , i2 , 1 , 2 ) = y( i1+0 , o2 );
v( i1 , i2 , 2 , 2 ) = y( i1+0 , o2 );
v( i1 , i2 , 3 , 2 ) = y( i1-1 , o2 );
v( i1 , i2 , 4 , 2 ) = y( i1-1 , o2 );

v( i1 , i2 , : , 3 ) = z( : , : , o3 );

if is_d
   if is_c
      d(i1,i2,:) = c(:,:,o3);
   else
      d(i1,i2,:) = z(:,:,o3);
   end
end

%---------------------------------------------------
% Fill Points around ZERO

%-------------------------------------- 
% Upper Row, close at  1 -- 2 
%-------------------------------------- 
v(  1 , i2 , 1 , 1 ) = x( i2-1 );
v(  1 , i2 , 2 , 1 ) = x( i2+0 );

v(  1 , i2 , [1 2] , 2 ) = y(1);
v(  1 , i2 , [1 2] , 3 ) = 0;

%-------------------------------------- 
% Lower Row, close at  4 -- 3
%-------------------------------------- 
v( g1 , i2 , 3 , 1 ) = x( i2+0 );
v( g1 , i2 , 4 , 1 ) = x( i2-1 );

v( g1 , i2 , [3 4] , 2 ) = y( s1+1 );
v( g1 , i2 , [3 4] , 3 ) = 0;

%-------------------------------------- 
% Left Column, close at   2 -- 3
%-------------------------------------- 
v( i1 ,  1 , 2 , 2 ) = y( i1+0 );
v( i1 ,  1 , 3 , 2 ) = y( i1-1 );

v( i1 ,  1 , [2 3] , 1 ) = x(1);
v( i1 ,  1 , [2 3] , 3 ) = 0;

%-------------------------------------- 
% Right Column, close at   1 -- 4
%-------------------------------------- 
v( i1 , g2 , 4 , 2 ) = y( i1-1 );
v( i1 , g2 , 1 , 2 ) = y( i1+0 );

v( i1 , g2 , [1 4] , 1 ) = x( s2+1 );
v( i1 , g2 , [1 4] , 3 ) = 0;


%---------------------------------------------------
% 4D --> 3-Column-Matrice:  v == [ X Y Z ]

v = reshape( v , g1*g2*g3 , g4 );

%-------------------------------------- 
% Remove not used Points arround

bad = ( sum(isnan(v),2) == g4 );

v(find(bad),:) = [];

%-------------------------------------- 
% Z == NaN --> ZERO

n3n = isnan( v(:,3) );

v( find(n3n) , 3 ) = 0;

%-------------------------------------- 

if ( Nout == 1 ) 

   return

end

%----------------------------------------------------

if is_d

   d            = reshape( d , g1*g2*g3 , 1 );
   d(find(bad)) = [];

end

 
%***************************************************
% Faces (Flaechen)

% f = zeros( ( s1*s2 + (s1+1)*s2 + (s2+1)*s1 ) , g3 );

f = zeros( ( 3*s12 + s1 + s2 ) , g3 );


%---------------------------------------------------
% Horizontal Faces around Z-Points:  
%
%   1 ------- 2
%             |
%      Z(ii)  |
%             |
%   4 ------- 3
%
%-------------------------------------- 

ni = s12;
ii = ( 1 : ni );

f(ii,:) = reshape( ( i1*o2 + o1*(i2-1)*g1 ) , ni , 1 ) * o3;

f(ii,:) = f(ii,:) + ones(ni,1) * i3 * g12;


%---------------------------------------------------
% Vertical Faces along Column:
%
%       Z(ii)
%   1 ---------> 2
%                |
%                |
%                v
%   4 <--------- 3
%       Z(ii+1)
%
%-------------------------------------- 

ni = s12 + s2;

ii = ( 1 : ni ) + s12;

f(ii,:) = reshape( ( j1*o2 + q1*(i2-1)*g1 ) , ni , 1 ) * o3;

f(ii,[3 4]) = f(ii,[3 4]) + 1;

f(ii,:) = f(ii,:) + ones(ni,1) * i3 * g12;


%---------------------------------------------------
% Vertical Faces along Row:
%
%         3 -----> 4
%         |        |
%   Z(ii) |        |  Z(ii+1)
%         |        | 
%         2        1
%
%-------------------------------------- 

ni = s12 + s1;

ii = ( 1 : ni ) + s12 + s12 + s2;

f(ii,:) = reshape( ( i1*q2 + ( o1*j2 - 1 )*g1 )' , ni , 1 ) * o3;

f(ii,[1 4]) = f(ii,[1 4]) + g1;

f(ii,:) = f(ii,:) + ones(ni,1) * i3 * g12;


%---------------------------------------------------
% Remove Faces with not used Points arround 

bad = cumsum(bad);

f = f - bad(f);

%---------------------------------------------------
% Remove Faces with only NaN

if is_d

  n3n = isnan(d);

end

f( find( sum(n3n(f),2) == g3 ) , : ) = []; 


if is_d & ~is_c

   d(find(isnan(d))) = 0;

end