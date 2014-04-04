function  [p,d] = geosect(p,varargin)

% GEOSECT Returns Coordinates of a Geodesic Section
%
% [Points,Distance] = GEOSECT( POS , Mode , N );
%
% POS = Start- and EndCoord of the Section
% POS = [ Lat1 Lon1 ; Lat2 Lon2 ]; 
%
% Mode = Mode for N: '#' | '<' | '>'
%
% '#'  N == Number of Points
% '<'  N == max. Distance [nm] between points
% '>'  N == min. Distance [nm] between points
%
% The Coordinates of the SectionPints are estimated
% on a geodesic line on a ellipsoid (see GEODIST).
%
% Use a imaginary Value of N for a linear Section in
% a cylindrical Projection. 
% In this case, the real Part is the MercatorParameter.
%
%   0.0 + N*i => linear in equidistant cylindrical Projection
%   0.8 + N*i => linear in Millers Projection
%   1.0 + N*i => linear in Mercator Projection
%
% Points   = Coordinates of Section in [ Lat Lon ]
% Distance = Distance [nm] along Section
%
% See also: GEODIST, MERCATOR
%

mode = { '#' '<' '>' };

N   = 2;
ext = [];

Nin  = nargin;
Nout = nargout;

%****************************************************************
% Check Inputs

if Nin < 1
   error('Not enough InputArguments.')
end

[msg,N,ext,mode] = checkin(p,N,ext,mode,varargin{:});

if ~isempty(msg)
   error(msg);
end

%****************************************************************
% convert to [ -180 .. 180 )

% p1(2) = p1(2) - 360 * floor( (p1(2)+180) / 360 );
% p2(2) = p2(2) - 360 * floor( (p2(2)+180) / 360 );


%****************************************************************

is_lin = ~isempty(ext);
is_num =  strcmp(mode,'#');

is_dist = ( ~is_num | ( Nout == 2 ) );

p0 = p;

%****************************************************************
% Estimate linear Distance

if is_lin

   p(:,1) = deg2merc(p(:,1),ext);

   if  is_dist
 
      nn = 1000;

      nq = ( 0 : (nn-1) ) / (nn-1);
   
      nq = nq(:);

       q = zeros(nn,2);

       q(2:(nn-1),:) = p(ones(1,nn-2),:) + nq(2:(nn-1)) * diff(p,1,1);

       q(:,1)      = merc2deg(q(:,1),ext);

       q([1 nn],:) = p0;

       dq         = zeros(nn,1);

       dq(2:nn)   = geodist(q(1:(nn-1),:),q(2:nn,:));
 
       dq = cumsum(dq) / 1852;

       d = dq(nn);

   end

%****************************************************************
elseif is_dist

   d = geodist(p(1,:),p(2,:)) / 1852;

end

%****************************************************************
% Get Number of Points

if ~is_num 

   N = d / N;           % Number of Intervalls

   switch mode
     case '<'
         N = ceil(N);
     case '>'
         N = floor(N);
     otherwise
         N = round(N);
   end

   N = N + 1;           % Number of Points

end

%****************************************************************

if N <= 2

   p = p0;

   if is_dist
      d = cat( 1 , 0 , d );
   end

   return

end

%****************************************************************
  
if ~is_lin

   [d,lat,lon] = geodist(p(1,:)',p(2,:)',N);

    p = cat( 2 , lat , lon );

    p([1 N],:) = p0;

    d = d / 1852;

    return

end

%****************************************************************

  np = ( 0 : (N-1) ) / (N-1);
   
  np = np(:);

   p1 = p;

   p = zeros(N,2);

   p(2:(N-1),:) = p1(ones(1,N-2),:)+np(2:(N-1)) * diff(p1,1,1);

   p(:,1)     = merc2deg(p(:,1),ext);

   p([1 N],:) = p0;


if ( Nout == 2 )

   d = interp1( nq , dq , np );

end

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  y = deg2merc(y,ext)

if ext ~= 0

  y = 180/pi * log(abs(tan(pi/4 + y*ext*pi/180/2))) / ext;

end

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  y = merc2deg(y,ext)
 
if ext ~= 0

 y = 2*180/pi * ( atan(exp(y*ext*pi/180)) - pi/4 ) / ext; 

end

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,N,ext,mode] = checkin(p,N,ext,modes,varargin);

% CHECKIN  Check Inputs
%

Nin = nargin - 4;

msg = '';
nl  = char(10);

mode = '';
par  = [];

%****************************************************************
% Check Coordinates

if ~isequal(size(p),[2 2])

   msg = 'POS must define 2 Points in [ Lat Lon ].';

end

%****************************************************************
% Check Npoints and Method

if Nin > 0

  Nin = min( Nin , 2 );

  ok = zeros(1,Nin);

  for ii = 1 : Nin
      v = varargin{ii};
      ok(ii) = chkstr(v,0);
      if ok(ii)
         mode = v;
      else
         ok(ii) = isnumeric(v) & ( prod(size(v)) == 1 );
         if ok(ii)
            ok(ii) = 2 * isfinite(v);
            if ok(ii)
               par = v;
            end
         end
      end
  end

  if ( any( ok == 0 ) | ( sum(ok==1) > 1 ) | ( sum(ok==2) > 1 ) )

     msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                'Following Inputs must be any of DivisionParameter and Mode.' , nl , ...
                'Mode must be a String.'  )

  end


  if isequal(ok,1) & ~isempty(mode)
     [msg,v] = str2vec( mode , '012345679.i*+-/()' ); 
     if isempty(msg) & ~isempty(v)
        par = v(1);
     end
  end
      
  if ~isempty(par)
     if imag(par) == 0
        N = par;
     else
        N = imag(par);
      ext = real(par);
     end
  end

  if ~isempty(ext)
     if ~( ( 0 <= ext ) & ( ext <= 1 ) )
        msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                 'Extension must between 0 and 1.' );
     end
  end

  if ~( N > 0 )
        msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                 'PointNumber or Distance must be larger 0.' );
  end

end

if isempty(mode)  
   mode = modes{1};
else
   mode = mode(1);
end


%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = chkstr(str,opt)


% CHKSTR  Checks Input for String
%
%  ok = chkstr(str,Option)
%
%  Option ~= 0 ==> only true for nonempty Strings
%
%   default:   Option == 0  (true for empty Strings)
%

 
if nargin < 2
   opt = 0;
end

ok = ( strcmp( class(str) , 'char' )      & ...
       ( prod(size(str)) == size(str,2) ) & ...
       ( isequal(opt,0) | ~isempty(str) )         );

