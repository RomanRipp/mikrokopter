function [az,el,hdl] = lightang(varargin);

% LIGHTANG  Spherical position of a light
%
% Improved version of Matlab's LIGHTANGLE.
%
% LIGHTANG creates or positions a light using Azimuth and 
%  Elevation. The Azimuth (AZ) is the horizontal rotation and
%  Elevation (EL) is the vertical rotation (both in degrees).
%  The interpretation of Azimuth and Elevation are exactly 
%  the same as with the VIEW command. 
% 
%  If the light passed into lightangle is a local light, the 
%   distance between the light and the camera target is preserved.
%
%---------------------------------------------------------------------
% InputOrder
%
% LIGTHANG( [Handle] , [ AZ     EL ] , ... )
% LIGTHANG( [Handle] ,   AZ  ,  EL   , ... )
%
% Handle can be a single GraphicHandle of Type:
%
%   'light' | 'axes' | 'figure' 
%
% In case of missing Handle, a new Object(s) are created
%  in the specified Axes/CurrentAxes or Figure/CurrentFigure.
%
% Following Inputs can be Property-Value-Pairs for a Light-Object.
%
%---------------------------------------------------------------------
% OutputOrder in case of LightHandleInput
%
% A Single Output returns the Vector: [ AZ  EL ]
%
%     AZEL = LIGHTANG( LightHandle , ... )
%
% The multiple OutputOrder is:  AZ,  EL,  LightHandle
%
%     [ AZ , EL , LightHandle ] = LIGHTANG( LightHandle , ... )
%
%---------------------------------------------------------------------
% OutputOrder in case missing LightHandleInput
%
% A Single Output returns the new LightHandle
%
%     LightHandle = LIGHTANG( ... )
%
% The multiple Outputs return the Azimuth and Elevation:
%
%     [ LightHandle , [ AZ   ,   EL ] ] = LIGHTANG( ... )
%     [ LightHandle ,   AZ   ,   EL   ] = LIGHTANG( ... )
%
%---------------------------------------------------------------------
%
%   See also LIGHT, CAMLIGHT, LIGHTING, MATERIAL, VIEW
%
 
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 1.5 $  $Date: 1998/07/13 16:41:11 $
%   
%   Modified by Christian Begler 2001/09/15
%                                2005/02/04
%

acc  = 1e-3;      % Accuracy for Angles

Nin  = nargin;
Nout = nargout;

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end

app = upper(fcn);
tag = sprintf('CREATED_BY_%s',upper(fcn));

%********************************************************
% Check Inputs

% 1: [ AZ  EL ]
% 1:   AZ        2:   EL
% 1:   H         2: [ AZ  EL ]
% 1:   AZ        2:   EL
% 1:   H         2:   AZ        3:  EL

% Valid Types of 1. HandleInput
tps = { 'figure'  'axes'  'light'  };

typ = '';
hdl = [];

az  = [];
el  = [];

new = zeros(0,1);

iv  = 0;

for ii = 1 : min(Nin,3)

    iv = ii + 1;

     v = varargin{ii};
     n = prod(size(v));

     if ~( isnumeric(v) & ( 1 <= n ) & ( n <= 2-(ii==3) ) );
         iv = iv - 1;
         break
     end

     %-----------------------------------------------
     % Check for [ AZ  EL ] (1. or 2. Input)

     if ( n == 2 )
        if isempty(az)
           az = v(1);
           el = v(2);
        end
        break
     end
 
     %-----------------------------------------------
     % Check 1. Input for Handle: Light / Axes / Figure

     ok = ( ( ii == 1 ) & ishandle(v) );
     if ok
        typ = get(v,'type');
        ok  = any(strcmp(typ,tps));
     end

     if ok
        hdl = v;
     elseif isempty(az)
        az = v;
     elseif isempty(el)
        el = v;
     else
        iv = iv - 1;
        break
     end

end

varargin = varargin( iv : Nin );

Nin = Nin - iv + 1;

if xor(isempty(az),isempty(el))
   error('Azimuth AND Elevation must be defined.');
end

if ~( isfinite(az) & isfinite(el) )
   error('Azimuth and Elevation must be finite.');
end

%********************************************************
% Check to Create new Handle

isg = isempty(az);   % Get AzEl

if isempty(typ)
   typ = 'root';
   hdl = 0;
   tps = cat( 2 , {typ} , tps );
end

isn = ~strcmp(typ,'light');   % New Light

nt  = size(tps,2);

it  = find(strcmp(typ,tps)) + 1;

for ii = it : nt

    typ = tps{ii};

    isl = strcmp(typ,'light');
    par = hdl;

    if isl
       hdl = findobj( par , 'type' , typ );
    else
       hdl = get( par , [ 'current'  typ ] );
    end

    if isempty(hdl) 

       if isg
          return
       end

       hdl = feval( typ , 'parent' , par , 'tag' , tag );

       if ~isl
           set( hdl , 'nextplot' , 'add' );
           set( par , [ 'current' typ ] , hdl );
       end

       new = cat( 1 , new , hdl );

    elseif isl

       hdl = hdl(1);   % CurrentLight

    end

end

%********************************************************
% Set VARARGIN

if Nin > 0
   try
      set( hdl , varargin{:} );
   catch
      if ~isempty(new)
          try, delete(new); end
      end
      error(sprintf('Error set Properties to "%s".\n%s\n',epsstr(hdl),lasterr));
   end
end

%********************************************************
% Set / Get LightPosition

ax = get( hdl , 'parent' );

ct = get( ax  , 'cameratarget'    );
rt = get( ax  , 'dataaspectratio' );

ct = ct * strcmp( get(hdl,'style') , 'local' );

dp = ( get( hdl , 'position' ) - ct ) ./ rt;

dn = norm(dp);    % [ 1 by 1 ]

%--------------------------------------
if ~isg     % Set

   if mod( el-90 , 180 ) == 0
      pos = [ 0  0  sin(el) ];
   else
      d2r = pi/180;
      azr = az * d2r;
      elr = el * d2r;
      pos = [ sin(azr)*cos(elr)  -cos(azr)*cos(elr)  sin(elr) ];
   end

   pos = (pos*dn) .* rt + ct;

   set( hdl , 'position' , pos );

   dp = ( pos - ct ) ./ rt;

   dn = norm(dp);    % [ 1 by 1 ]

end

%--------------------------------------
% Get

  r2d = 180/pi;

  az =  r2d * ( atan2( dp(2) , dp(1) ) + pi/2 );
  el =  r2d *    asin( dp(3) / dn );

  az = acc * round( az / acc );
  el = acc * round( el / acc );

  az = az - 360 * floor(az/360); % [    0 .. 360 )
  az = az - 360 * ( az > 180 );  % ( -180 .. 180 ] 

%********************************************************
% Check OutputOrder for New Light

if isn

   h = hdl;

   if Nout <= 2
      el = [ az  el ];
   else
      hdl = el;
       el = az;
   end

   az = h;

elseif Nout <= 1

   az = [ az  el ];

end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [ h , form ] = epsstr( h );

% EPSSTR  Transform Number exact into String,
%          
% using Matlab's floating point relative accuracy.
%
%  Form = EPSSTR;   
%    returns Format for using with SPRINTF
%
%  [ String , Form ] = EPSSTR( Number ); 
%    returns exact String for Number
%
%  Form = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 )
%


form = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 );

if nargin < 1

  h = form;

  return

end


if ~isnumeric(h)  |  ( prod(size(h)) > 1 )
 error('Handle must be a single Numeric.');
end


  h = sprintf(form,h);
