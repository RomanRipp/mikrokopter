function  c = xy_rot(c,varargin)

% XY_ROT   Rotates Matrix in X-Y-LAYER
%
%   C = XY_ROT( C , Angle , CenterXY , OutsideValues )
%
%   C             = [ nY by nX by n3 ]
%
%   CenterXY      = [  1 by  2 ] = [ CenterX CenterY ]
% 
%   OutsideValues = [  1 by n3 ] 
%

Nin = nargin;

if Nin == 0
   error('Not enough InputArguments.')
end

if ( Nin < 2 ) | isempty(c)
   return
end

s = size(c);

s12 = s(1) * s(2);

if any( s([1 2]) == 1 )
   error('C must be at least [ 2 by 2 ].')
end
 
ns = size(s,2);

if ns > 3

   c = reshape( c , s(1) , s(2) , prod(s(3:ns)) );

elseif ns < 3

   s = cat( 2 , s , -1 );
 
end

%-------------------------------------

s3  = size(c,3);

cc  = s([2 1]) / 2;
phi = 0;
bg  = NaN*ones(1,s3);

for ii = 1 : Nin-1

    v = varargin{ii};

    sv = size(v);

    ok = ( isnumeric(v)  &  ~isempty(v)  &  ( prod(sv) == sv(2) ) );
    if ok & any( sv(2) == [ 1  2 ] )
       ok = all(isfinite(v)) ;
    end

    if ok
       switch sv(2)
         case 1
              phi = v;
         case 2
              cc = v;
         case s3
              bg = v;
         case s(3)
              bg = v;
              ind = ( 1 : prod(s(3:ns)) );
              ind = ind - s(3) * floor( ind/s(3) );
              ind = ind + s(3) * ( ind == 0 );
              bg  = bg(ind);
       end
    end
                
end

%-------------------------------------

phi = phi - 360 * floor(phi/360);
 
if phi == 0
   return
end

phi = phi * pi / 180;

%-------------------------------------
% Basis CoordinateSystem

ex0 = [ 1 ; 0 ];
ey0 = [ 0 ; 1 ];

%-------------------------------------
% New CoordinateAxes in Basis System

ex1 = [  cos(phi) ; sin(phi) ];
ey1 = [ -sin(phi) ; cos(phi) ];

%--------------------------------------------
% TransformationMatrice

E0 = [  ex0(:)  ey0(:)  ];
E1 = [  ex1(:)  ey1(:)  ];

T = E0  * inv( E1 );

%--------------------------------------------

x = ( 1 : s(2) ) - cc(1);
y = ( 1 : s(1) ) - cc(2);

xy = T * cat( 1 , reshape( ones(s(1),1) * x    , 1 , s12 ) , ...
                  reshape( y(:) * ones(1,s(2)) , 1 , s12 ) );


xy = reshape( permute(xy,[2 1]) , s(1) , s(2) , 2 );
         
%--------------------------------------------
% Interpolate x,y,c --> xy
%
%  for ii = 1 : s3
%      c(:,:,ii) = interp2( x , y , c(:,:,ii) , xy(:,:,1) , xy(:,:,2) , 'linear' );
%  end
%

  xy(:,:,1) = 1 + ( xy(:,:,1) - x(1) ) / ( x(s(2)) - x(1) ) * (s(2)-1); % Column
  xy(:,:,2) = 1 + ( xy(:,:,2) - y(1) ) / ( y(s(1)) - y(1) ) * (s(1)-1); % Row

  % Check for out of range values and set to 1
  bad = find( ( xy(:,:,1) < 1 ) | ( xy(:,:,1) > s(2) ) | ...
              ( xy(:,:,2) < 1 ) | ( xy(:,:,2) > s(1) )       );

  if ~isempty(bad)
      xy(bad    ) = 1;
      xy(bad+s12) = 1;
  end

  % Matrix element indexing
  ind = floor(xy(:,:,2)) + floor(xy(:,:,1)-1) * s(1);

  % Check for BorderValues
  b        = ones(s(1),s(2),2);
  b(:,:,1) = s(2);
  b(:,:,2) = s(1);

  bb = find( xy == b );

  xy = xy - floor(xy);

  if ~isempty(bb)

    b(:,:,2) = 1;

      xy(bb) =  xy(bb) + 1;
     ind(bb) = ind(bb) - b(bb);

  end

  xy = xy(:,:,[1 2 1 2]);

  xy(:,:,[3 4]) = 1 - xy(:,:,[3 4]);
  
  xy = xy(:,:,[4 2 4 2]) .* xy(:,:,[3 3 1 1]);
   
  for ii = 1 : s3

     ind = ind + ( ii > 1 ) * s12;

       d =  c(ind)      .* xy(:,:,1) + c(ind+1)      .* xy(:,:,2) + ...
            c(ind+s(1)) .* xy(:,:,3) + c(ind+s(1)+1) .* xy(:,:,4);

    d(bad)     = bg(ii);

    c(:,:,ii ) = d;

  end

%--------------------------------------------
% Reshape back

if ns > 3

   s = num2cell(s);

   c = reshape( c , s{:} );

end
