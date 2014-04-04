function hm = getcdata(h,action,mode);

% GETCDATA creates a Menu and Show X-Y_Position
%
% MenuHandle = GETCDATA( SurfaceHandle , Form , Mode )
%
%   creates Menu to Display [ X Y C ] of the PointerPosition
%
%
%  GETCDATA( MenuHandle , 'on'    )    
%  GETCDATA( MenuHandle , 'off'   ) 
%  GETCDATA( MenuHandle , 'delete') 
%
% internal for ButtonDownFcn:
%
%    GETCDATA('down',SurfaceHandle,MenuHandle)
%
%

 hm = [];

 Nin = nargin;
 if Nin < 1
   return
 end

 ok = ( isnumeric(h) & ( prod(size(h)) == 1 ) );
 if ok
    ok = ishandle(h);
    if ok
       typ = get(h,'type');
       ok = any( strcmp( typ , { 'uimenu'  'surface'  'image' } ) );
    end
 end

 if ~ok
    return
 end

 
switch typ

 %**********************************************
 case { 'surface'  'image' }

 
   if Nin < 2
      form = ' %g ';
   else
     if ( ischar(action) & ~isempty(action) & ...
         ( prod(size(action)) == size(action,2) )   )
       form = action;
     else
       form = ' %g ';
     end
   end

   if Nin < 3
      mode = 'none';
   else
     if ~( ischar(mode) & ~isempty(mode) & ...
         ( prod(size(mode)) == size(mode,2) )   )
       mode = 'none';
     end
   end

   %------------------------------------------------
   % Create Menu

   axe = get( h , 'parent' );
   fig = get( axe , 'parent' );

    hm = uimenu('parent'   , fig , ...
                'label'    , ''  , ...
                'tag'      , 'CDATA_MENU'     );


   %------------------------------------------------

   Hform = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 );

   hh    = sprintf( Hform , hm );

   DownFcn = [ 'getcdata(' hh ',''down'');' ];

   ud = struct( 'Handle'  , { h    }    , ...
                'Type'    , { typ  }    , ...
                'Form'    , { form }    , ...
                'Mode'    , { mode }    , ...
                'DownFcn' , { DownFcn }       );


    set( hm , 'userdata' , ud );

    getcdata(hm,'on');

    
 %**********************************************
 case 'uimenu'

   if ~strcmp( get(h,'tag') , 'CDATA_MENU' )
      return
   end

   if ~( ischar(action) & ~isempty(action) & ...
         ( prod(size(action)) == size(action,2) )   )
      return
   end

   ud = get(h,'userdata');

   switch action

     %-------------------------------------------
     case 'on'

        set( ud.Handle , 'ButtonDownFcn' , ud.DownFcn );

     %-------------------------------------------
     case 'off'

        set( ud.Handle , 'ButtonDownFcn' , '' );
   
        set( h , 'label' , '' );

     %-------------------------------------------
     case 'delete'

        getcdata( h , 'off' );

        delete( h );

     %-------------------------------------------
     case 'down'

        axe = get( ud.Handle , 'parent' );
        
         cp = get( axe , 'CurrentPoint' );

         cp = cp(1,[1 2]);

          x = get( ud.Handle , 'xdata' );
          y = get( ud.Handle , 'ydata' );
          c = get( ud.Handle , 'cdata' );

          hm = [ cp NaN ];

          if isempty(c)
             return
          end

          if strcmp( ud.Type , 'image' )
             x = linspace( x(1) , x(end) , size(c,2) );
             y = linspace( y(1) , y(end) , size(c,1) );
          end

          %-----------------------------------
          % Check if CurrentPoint in Range

          if  ( cp(1) < min(x(:)) )  |  ...
              ( cp(1) > max(x(:)) )  |  ...
              ( cp(2) < min(y(:)) )  |  ...
              ( cp(2) > max(y(:)) )
              return
          end

          %-----------------------------------
          % Check for equal Grid

          if isequal( size(x) , size(c) )
             d = x - x(ones(1,size(x,1)),:);                      
             if all( d(:) == 0 );
                x = x(1,:); 
             end
          end

          if isequal( size(y) , size(c) )
             d = y - y(:,ones(1,size(y,2)));                      
             if all( d(:) == 0 );
                y = y(:,1); 
             end
          end

          if isequal( [ prod(size(y)) prod(size(x)) ] , size(c) )

             hm(3) = interp21( x , y , c , cp(1) , cp(2) , ud.Mode );
 
          else

            if strcmp(mode,'geo')

                y     = deg2merc( y     , 1 );
                cp(2) = deg2merc( cp(2) , 1 );

            end

            if prod(size(x)) == size(c,2)
               x = x(:)';
            end
            if prod(size(y)) == size(c,1)
               y = y(:);
            end

            try
              hm(3) = interp2( x , y , c , cp(1) , cp(2) );
            end

          end          

          set( h , 'label' , sprintf(ud.Form,hm) );
   end

end


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function zi = interp21(x,y,z,xi,yi,mode);

% INTERP21  interpolates a Single Point into a 2-D function
%
%  ZI = INTERP21( X , Y , Z , XI , YI );
%
%   X = [  1 by Nx ]   GridVector for X 
%   Y = [ Ny by  1 ]   GridVector for Y
%   Z = [ Ny by Nx ]   MatriceFunction
%
%  XI = [ 1 by 1 ]  Target X
%  YI = [ 1 by 1 ]  Target Y
% 
%  ZI = [ 1 by 1 ]  
%
% If XI or YI are out of range, ZI will set to NaN.
%
% The Interpolation Method is an invers weight.
%
%  ZI = INTERP21( X , Y , Z , XI , YI , 'geo' );
%  
%  Interpretes X as Longitude, Y as Latitude to determine
%   invers weighted Distance
%

zi = [];

Msg = '';

nl = char(10);


if nargin < 5
 error('Number of Inputs must be 5.' )
end


if nargin < 6
  mode = 'none';
end


x = x(:);
y = y(:);



if ~isequal( size(z) , [ size(y,1)  size(x,1) ] )
  Msg = 'Size of Z must be  [ length(Y)  by  length(X) ].';
end

if ~all( ( diff(x) < 0 )  |  ( diff(x) > 0 ) )  |  ...
   ~all( ( diff(y) < 0 )  |  ( diff(y) > 0 ) ) 
  Msg = [ Msg nl(1:(end*(~isempty(Msg))))   ...
           ' X and Y must be monotonic. ' ];
end

if ( prod(size(xi)) ~= 1  )  |  ...
   ( prod(size(yi)) ~= 1  )
  Msg = [ Msg nl(1:(end*(~isempty(Msg))))   ...
          'XI, YI must define a single Target.' ];
end

mode = mode(:)';

if ~ischar(mode) 
  Msg = [ Msg nl(1:(end*(~isempty(Msg))))   ...
            ' Additional Inputs  could be String ''geo''.' ];  
end


if ~isempty(Msg)
  error([ ' INTERP21: ' Msg ])
end

zi = NaN;


if strcmp(mode,'geo')

  y  = deg2merc( y  , 1 );
  yi = deg2merc( yi , 1 );

end


 fx = 1 - 2 * ( x(1) > x(end) );  % Monotonic
 fy = 1 - 2 * ( y(1) > y(end) );

 ix1 = sum( fx*x <= fx*xi );
 iy1 = sum( fy*y <= fy*yi );

 ix2 = size(x,1) - sum( fx*x >= fx*xi ) + 1;
 iy2 = size(y,1) - sum( fy*y >= fy*yi ) + 1;

if ( ix1 >= 1 ) & ( iy1 >= 1 ) & ...
   ( ix2 <= size(x,1)  ) &  ( iy2 <= size(y,1)  )      
          
  if ( ix1 == ix2 ) & ( iy1 == iy2 )

    zi = z( iy1 , ix1 );

  else
       
    ix = [ ix1 ; ix2 ; ix2 ; ix1 ];
    iy = [ iy1 ; iy1 ; iy2 ; iy2 ];

    x = x(ix);
    y = y(iy);
    z = z( (ix-1)*size(z,1) + iy );


    % Normalized to Grid !!!

    gx = ( x(1) - x(2) );
    gy = ( y(2) - y(3) );

    dx = abs( ( x - xi ) / ( gx + ( gx == 0 ) ) );  
    dy = abs( ( y - yi ) / ( gy + ( gy == 0 ) ) );

    d  = sqrt( dx.^2 + dy.^2 );
    d  = ( sqrt(2) - d ) .* ( 1 - dx ) .* ( 1 - dy );

    zi = sum( z.*d ) / sum(d);

  end
         
end


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function  y = deg2merc(y,ext)

 %  function  y = deg2merc(y,ext)

 if ext 
  y = 180/pi * log(abs(tan(pi/4 + y*ext*pi/180/2))) / ext;
 end
