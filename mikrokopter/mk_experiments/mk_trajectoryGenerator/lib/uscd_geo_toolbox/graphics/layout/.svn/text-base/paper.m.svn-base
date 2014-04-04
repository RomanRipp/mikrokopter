function [pappos,figpos] = paper(varargin);

% PAPER   Set nice PaperFilling Position of Figure
%
% PAPER( Fig , Orient , Offset , [PaperType] , [PaperUnits] )
%
%  Orient = {'portrait'}  |  'landscape'
%  Offset = [ Left  Bottom  Right  Top ]
%
%  PaperTypes: disp(set(0,'DefaultFigurePaperType'))
%  PaperUnits: disp(set(0,'DefaultFigurePaperUnit'))
%
% PAPER( Fig , Axe , Offset , ... )
% 
%  Determines Orientation from Extension of Axes
%
% [ PaperPosition , FigurePosition ] = PAPER( ... )
% 
% Returns PaperPositon in [PaperUnits] and FigurePosition in [pixels]
%

Vin = varargin;
Vin = Vin(:);

pappos = NaN*ones(1,4);
figpos = NaN*ones(1,4);

offset  = [ 1.0  1.5  0.8  1.5 ];       % PaperOffset: [ Left Bottom Top Right ]

modus = lower(set(0,'defaultfigurepaperorientation')); % PaperOrientation
types = lower(set(0,'defaultfigurepapertype'));        % PaperTypes
units = lower(set(0,'defaultfigurepaperunits'));       % PaperUnits

%*****************************************************************
% Check Inputs

[fig,Vin]  = hndlvarg(Vin,size(Vin,1),'figure');
if isempty(fig)
   return
end

[axe,Vin]  = hndlvarg(Vin,prod(size(Vin)),'axes');

if ~isequal( get(axe,'parent') , fig )
    axe = [];
end

off  = [];
mode = '';
typ  = '';
uni  = '';

for ii = 1 : prod(size(Vin))
    if chkstr(Vin{ii})
       str = Vin{ii};
       if ~isempty(str)
           str = lower(str);
           jj  = strcmp(str,modus);
           kk  = strcmp(str,types);
           ll  = strcmp(str,units);
           if     any(jj)
               mode = modus{min(find(jj))};
           elseif any(kk)
               typ  = types{min(find(kk))};
           elseif any(ll)
               uni  = units{min(find(ll))};
           else
               jj  = strmatch(str,modus);
               kk  = strmatch(str,types);
               ll  = strmatch(str,units);
               if     ~isempty(jj)
                   mode = modus{min(min(jj))};
               elseif ~isempty(kk)
                   typ  = types{min(min(kk))};
               elseif ~isempty(ll)
                   uni  = units{min(min(ll))};
               else
                   error(sprintf('Invalid Input: %s',Vin{ii}))
               end
           end
       end
    elseif isnumeric(Vin{ii})
       off = Vin{ii};
       if ~isempty(off)
           if ~all( ~isnan(off) & isfinite(off) )
               error('Value for Offset must be finite and not NaN.');
           end
       end
    elseif ~isempty(Vin{ii})
       error(sprintf('Invalid Type "%s" of Input.',class(Vin{ii})));
    end
end                    

%*****************************************************************
% Check Offset

if isempty(off)

   off = offset;

else

   off = off(:)';

   noff = size(off,2);

   off = off(1:min(4,noff));

   switch noff
     case 1
       off = off([1 1 1 1]);
     case 2
       off = off([1 2 1 2]);
     case 3
       off = off([1 2 3 2]);
   end

end

%*****************************************************************
% Check Mode

if isempty(mode)
   is_land = ~isempty(axe);
   if is_land
      [upp,psi] = ppunit(axe);
      is_land = ( psi(3) > psi(4) );
   end
   mode = modus{ 1 + is_land };
end

if isempty(typ)
   typ = get(fig,'papertype');
end

if isempty(uni)
   uni = get(fig,'paperunits');
end
%*****************************************************************
% Determine PaperPosition

set( fig , 'papertype'        ,  typ     , ...
           'paperorientation' ,  mode    , ...
           'paperunits'       ,  uni           );
 
papsi = get(fig,'papersize');

ind = [ 1  2
        3  4 ];

pappos = cat( 2 , off([1 2]) , papsi - sum(off(ind),1) );

set(fig,'paperposition',pappos)

%*****************************************************************
% Set FigurePosition

figuni = get(fig,'units');

set( fig , 'units' , 'pixels' );
 
fp = get(fig,'position');

wysiwyg(fig);

figpos = get(fig,'position');

figpos(1) = fp(1);
figpos(2) = fp(2) + fp(4) - figpos(4);

if figpos(2) < 0

   ssi = get(0,'screensize');

   mb = get(fig,'menubar');
   tb = get(fig,'toolbar');

   mb =   strcmp(mb,'figure');
   tb = ( strcmp(tb,'figure') | ( strcmp(tb,'auto') & mb ) );

   if ~mb
      shh = get(0,'showhiddenhandles');
            set(0,'showhiddenhandles','on');
      hm = findobj(fig,'type','uimenu','visible','on','parent',fig);
            set(0,'showhiddenhandles',shh);
      mb = ~isempty(hm);
   end

   top = ssi(4) - 25 * ( 1 + mb + tb );

   if figpos(2)+figpos(4) < top

      figpos(2) = top - figpos(4);

   end

end

set( fig , 'position' , figpos );

set( fig , 'units' , figuni );


if nargout == 0
   clear('pappos','figpos');
end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function wysiwyg(fig)

% WYSIWYG  WhatYouSeeIsWhatYouGet
%
% WYGIWYS( FigureHandle )
%
%WYSIWYG -- changes the size of the figure on the screen to equal
%       the size of the figure that would be printed, 
%       according to the papersize attribute.  Use this function
%       to give a more accurate picture of what will be 
%       printed.
%       Dan(K) Braithwaite, Dept. of Hydrology U.of.A  11/93
 
if nargin < 1
  fig = get(0,'currentfigure');
end

if isempty(fig)
 return
end

ok = ( isnumeric(fig)  &  ( prod(size(fig)) == 1 ) );
if ok
   ok = ishandle(fig);
   if ok
      ok = strcmp( get(fig,'type') , 'figure' );
   end
end

if ~ok
   error('Input must be a FigureHandle.');
end

unis = get(fig,'units');
ppos = get(fig,'paperposition');

org = get(fig,'position');

set(fig,'units',get(fig,'paperunits'));

pos = get(fig,'position');
pos(3:4) = ppos(3:4);

set(fig,'position',pos);
set(fig,'units',unis);

new = get(fig,'position');

if ~all( new([1 2]) == org([1 2]) )
    new([1 2]) = org([1 2]);
    set(fig,'position',new);
end
  
%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [upp,axesize,x,y] = ppunit(axe);

% PPUNIT  Returns AxesUnitPerPixel and AxesPixelSize
%
%  [ UnitPerPixel, PixelSize, X , Y ] = PPUNIT( AxesHandle )
%
%  UnitPerPixel = [ XUnitsPerPixel  YUnitsPerPixel  ZUnitsPerPixel ] ;
%  PixelSize    = [ PixelLeft  PixelBottom  PixelWidth  PixelHight ];
%
%
%  Units --> Points:   
% 
%    [points] = [unit] / UnitsPerPixel / ScreenPixelPerInch * PointsPerInch
%
%   PointsPerInch = 72
%
%

if nargin == 0
   fig = get(0,'currentfigure');
   axe = get(fig,'currentaxes');
end

if isempty(axe)

   upp     = zeros(0,2);
   axesize = zeros(0,4);

   x = [];
   y = [];

   return

end


ok = ( isnumeric(axe) & ( prod(size(axe)) == 1 ) );
if ~ok
    ok = ishandle(axe);
    if ok
       ok = strcmp( get(axe,'type') , 'axes' );
    end
end

if ~ok
   error('Input must be an AxesHandle.');
end


  
axeuni  = get(axe,'units');
          set(axe,'units','pixels')
axesize = get(axe,'position');
          set(axe,'units',axeuni);


%********************************************
if ~isequal( get(axe,'view') , [ 0  90 ] );

  [upp,axesize,x,y] = ppunit3(axe,axesize);

  return

end
%********************************************


dx = diff(get(axe,'xlim'));
dy = diff(get(axe,'ylim'));

mode   = get(axe,'dataaspectratiomode');

if strcmp(mode,'manual')

   aspect = get(axe,'dataaspectratio');

   %  w/dx*ax == h/dy*ay 
   %
   % w = h * dx/dy * ay/ax;
   % h = w * dy/dx * ax/ay; 
   %

   pos = zeros(2,2);

   pos(1,2) = axesize(4);
   pos(1,1) = axesize(4) * dx/dy * aspect(2)/aspect(1);
   pos(2,1) = axesize(3);
   pos(2,2) = axesize(3) * dy/dx * aspect(1)/aspect(2);

   ii = find( sum( ( pos <= axesize([1 1],[3 4]) ) , 2 ) == 2 );

   pos = pos(ii(1),:);
 
   axesize([1 2]) = axesize([1 2])+axesize([3 4])/2-pos/2;
   axesize([3 4]) = pos;

end

 upp = [ dx  dy  0 ] ./ axesize([3 4 4]) ; % UnitPerPixel

%-----------------------------------
% Index for Cube

  ix = [ 1 2 2 1 1 2 2 1 ];
  iy = [ 1 1 2 2 1 1 2 2 ];

  x = [ 0  axesize(3) ];
  y = [ 0  axesize(4) ];
  
  x = x(ix);
  y = y(iy);

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [upp,axesize,x,y] = ppunit3(axe,axesize);

% PPUNIT  Returns AxesUnitPerPixel and AxesPixelSize for 3D-Axes
%


  xl = get(axe,'xlim');
  yl = get(axe,'ylim');
  zl = get(axe,'zlim');

   v = get(axe,'view');

  dx = diff(xl);
  dy = diff(yl);
  dz = diff(zl);


%-----------------------------------
% Index for Cube

  ix = [ 1 2 2 1 1 2 2 1 ];
  iy = [ 1 1 2 2 1 1 2 2 ];
  iz = [ 1 1 1 1 2 2 2 2 ];


  % Index for Axes in ix, iy, iz
  %  i01(1,:) --> i01(2,:)
  %       x  y  z

  i01 = [ 1  1  1 
          2  4  5 ];

%-----------------------------------

mode = get(axe,'dataaspectratiomode');

if strcmp(mode,'manual')

 x = xl(ix);
 y = yl(iy);
 z = zl(iz);

 aspect = get(axe,'dataaspectratio');

else

 x = ix-1;  % 0 | 1
 y = iy-1;
 z = iz-1;

 aspect = [ 1  1  1 ];

end

[x,y] = xy_proj( x , y , z , v , aspect );

% Normalize
x = x - min(x);
y = y - min(y);

xr = max(x);  % Range
yr = max(y);  % Range


if strcmp(mode,'manual')

   %  w/dx*ax == h/dy*ay 
   %
   % w = h * dx/dy * ay/ax;
   % h = w * dy/dx * ax/ay; 
   %

   pos = zeros(2,2);

   pos(1,2) = axesize(4);
   pos(1,1) = axesize(4) * xr/yr;
   pos(2,1) = axesize(3);
   pos(2,2) = axesize(3) * yr/xr;

   ii = find( sum( ( pos <= axesize([1 1],[3 4]) ) , 2 ) == 2 );

   pos = pos(ii(1),:);

   axesize([1 2]) = axesize([1 2])+axesize([3 4])/2-pos/2;
   axesize([3 4]) = pos;

end

 % Cube in Pixel

 x = x/xr * axesize(3);
 y = y/yr * axesize(4);

 % AxesLength in Pixel

 d01 = sqrt( ( x(i01(1,:)) - x(i01(2,:)) ) .^2 + ...
             ( y(i01(1,:)) - y(i01(2,:)) ) .^2   );
       

 upp = ~( d01 < 1 ) .* [ dx  dy  dz ] ./ ( d01 + ( d01 < 1 ) );

 
%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [x,y,z] = xy_proj(x,y,z,v,aspect)

% function  [x,y,z] = xy_proj(X,Y,Z,View,DataAspectRatio)
%
% View = [ Azimut Elevation ]
%

si = num2cell(size(x));

az = v(1)*pi/180; 
el = v(2)*pi/180;

aspect = aspect(:)' / min(aspect(:));

x = x(:)';
y = y(:)';
z = z(:)';

T = [ cos(az)           sin(az)           0
     -sin(el)*sin(az)   sin(el)*cos(az)  cos(el)
      cos(el)*sin(az)  -cos(el)*cos(az)  sin(el) ];

T = T ./ aspect(ones(3,1),:);

xyz = T * [ x ; y ; z ];

x = reshape(xyz(1,:),si{:});
y = reshape(xyz(2,:),si{:});
z = reshape(xyz(3,:),si{:});
 

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [h,v,nr] = hndlvarg(v,n,typ,tag);

% HNDLVARG  Checks Input-varargin for Handle
%
%  [H,V,InputNr] = HANDLVARG( V , n , Type , Tag )
%
%  Search in V{n(1)} .. V{n(2)} for Handle of Type and Tag
%
%  If no Handle found and Type =  'figure' | 'axes', 
%    GCF | GCA is returned 
%
%  see also: CHKHNDL
%

Nin = nargin;

h  = [];
nr = 0;

if Nin == 0
   v = cell(1,0);
   return
end


if ~iscell(v)
   error('Input must be a CellArray.');
end

nv = prod(size(v));

if Nin < 2
   n = [];
end

if Nin < 3
   typ = 'figure';
end

chkin = { typ };

if Nin == 4
   chkin = cat( 2 , chkin , {tag} );
end

%------------------------------------------------------------

if isempty(n)
   n = [ 1  nv ];
else
   n = n(:)';
   n = n( 1 : max(1,size(n,2)) );
   n = min( n , nv );
   if size(n,2) == 1
      n = [ 1  n ];
   end
end


%-------------------------------------------------------------

ok = 0;

for nr = n(1) : n(2)

   ok = chkhndl(v{nr},chkin{:});

   if ok
      break
   end

end

%-------------------------------------------------------------

if ok

  h = v{nr};

  v(nr) = [];   % Reduce v

  return

end

nr = 0;

if ~any( strcmp( typ , { 'figure' 'axes' } ) )
   return
end

h = get( 0 , 'currentfigure' );
if isempty(h)
   return
end

if strcmp( typ , 'axes' )
   h = get( h , 'currentaxes' );
end

ok = chkhndl(h,chkin{:});

if ~ok
   h = [];
end



%***************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [ok,msg] = chkhndl(h,typ,tag);

% CHKHNDL(H,Type,Tag)  Checks, if H is a Handle of specified Type and Tag
%
%  Tag   CharArray or CellStringArray which Tags, the Handle has to be.
%          The Wildcard '*' can be used in the Strings for Tag
%

Nin = nargin;

ok  = 0;
msg = [];


if Nin == 0
   return
end

ok = ( isnumeric(h) & ( prod(size(h)) == 1 ) );

if ~ok
   return
end

ok = ishandle(h);
if ~ok | ( Nin < 2 )
   return
end

%-------------------------------------------------------------------------
% Check with Type

is_typ = chkstr(typ,0);

if is_typ

  ok = isempty(typ);
  if ~ok
      ok = strcmp( get(h,'type') , typ  );
  end

  if ~ok | ( Nin < 3 )
     return
  end

elseif ( Nin == 3 )

  msg = 'Input Type must be a String.';
  ok  = 0;
  return

else
 
  tag = typ;
  
end


%-------------------------------------------------------------------------
% Check Tag

[ok,tag] = chkcstr(tag,0);

if ~ok
    msg = 'Input Tag must be a CharArray or CellStringArray.';
    return
end   


%-------------------------------------------------------------------------
% Check with Tag

 t = get(h,'tag');

nt = size(t,2);

ok = strwcmp(t,tag);

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,cmp,wc,nn,cc] = strwcmp(str,cmp,wc,nn,cc)

% STRWCMP   Compare Strings, including WildCards
%
% OK = STRWCMP( STR , CMP , WildCard )
%
% STR: CharacterArray or CellStringArray to compare with Comp
%
% CMP: CharacterArray or CellStringArray with strings to compare with STR
%         strings can contains WildCards
%
% WildCard: specify WildCard to use, default: '*'
%
% OK : logical Array with same size of STR, contains 1 if 
%        any strings of CMP match string of STR
%
% Special Output:  [ok,cmp,wc,nn,cc] = STRWCMP( STR , CMP , [WildCard] )
%
%  to use in follwing statements: ok = STRWCMP( STR , cmp , wc , nn , cc );
%
%  which makes it a bit faster.
%
% see also: STRCMP, FINDSTR
%

Nin  = nargin;
Nout = nargout;

Msg = '';
 nl = char(10);

%***************************************************************
% Check Inputs

%---------------------------------------------------------------
% String

if Nin < 1
   str = cell(0,0);
end

if ~( iscell(str) & isempty(str) )
   [ok,str] = chkcstr(str,0);
   if ~ok
      Msg=[ Msg  nl(1:(end*(~isempty(Msg)))) ...
           'First Input must be a CharacterArray or CellStringArray.' ];
   end
end

%---------------------------------------------------------------
% CompareString

if Nin < 2
   cmp = cell(0,0);
end

if ~( iscell(cmp) & isempty(cmp) )
   [ok,cmp] = chkcstr(cmp,0);
   if ~ok
      Msg=[ Msg  nl(1:(end*(~isempty(Msg)))) ...
            'Second Input must be a CharacterArray or CellStringArray.' ];
   end
end

%---------------------------------------------------------------
% WildCard

if Nin < 3
   wc = '*';
elseif ~( ischar(wc) & ( prod(size(wc)) == 1 ) )
   Msg=[ Msg  nl(1:(end*(~isempty(Msg)))) ...
         'WildCard must be a Single Character.' ];
end
  
%---------------------------------------------------------------

if ~isempty(Msg)
   error(Msg)
end

%***************************************************************

si = size(str);

if ( isempty(str) | isempty(cmp) ) & ( Nout <= 1 ) 
   ok = zeros(si);
   return
end

cmp = cmp(:);

if any(strcmp(cmp,wc)) & ( Nout <= 1 )  % Wildcard only
   ok = ones(si);
   return
end

%***************************************************************
% Analyze CompareStrings

nc = size(cmp,1);

ok = ( Nin == 5 );

if ok
   ok = ( isequal(size(cc),[nc 1]) & iscell(cc) & ...
          isequal(size(nn),[nc 2]) & isnumeric(nn)    );
   if ok
      try
         ok = cat( 1 , cc{:} ); 
         ok = ( size(ok,2) == 3 );
      catch
         ok = 0;
      end
   end
end


%--------------------------------------------------
if ~ok
%--------------------------------------------------

  cc    = cell(nc,1);
  cc(:) = { zeros(0,3) };  % { [ Start End N  ] }

  nn = zeros(nc,2);        %   [ Ncmp  sum(N) ] 

  for ii = 1 : nc

    if ~isempty(cmp{ii})

       iwc = ( double(cmp{ii}) == double(wc) );

       if ~any( iwc )

          nn(ii,:) = size(cmp{ii},2);
          cc{ii}   = [ 1  nn(ii,:) ];

       else

         %--------------------------------------------
         % Remove Duplicate WildCards

         iwc = find( iwc );
         if ~( prod(size(iwc)) == 1 )
            jj = find( diff(iwc) == 1 );
            if ~isempty(jj)
               cmp{ii}(iwc(jj+1)) = [];
            end
         end

         %--------------------------------------------
         % Get Start End
     
         iwc = ( double(cmp{ii}) == double(wc) );
  
          n  = size(iwc,2);

          if ( n == 1 ) & ( iwc == 1 ) & ( Nout <= 1 ) % Wildcard only
             ok = ones(si);
             return
          end

          i0 = ~iwc(1);
          i1 = ~iwc(n);

         iwc = cat( 2 , ones(1,i0) , iwc , ones(1,i1) );

         iwc = find( iwc );

         iwc = iwc(:);

           n = size(iwc,1) - 1;

         cc{ii} = zeros(n,3);

         if n > 0      

            cc{ii}(:,[1 2]) = cat( 2 , iwc(1:n)+1 , iwc((1:n)+1)-1 ) - i0;

            cc{ii}(:,3) = cc{ii}(:,2) - cc{ii}(:,1) + 1;
 
         end

         nn(ii,:) = cat( 2 , size(cmp{ii},2) , sum(cc{ii}(:,3),1) );

       end

    end

  end

%--------------------------------------------------
end
%--------------------------------------------------

if ( Nout > 1 )

  if ( isempty(str) | isempty(cmp) )
     ok = zeros(si);
     return
  end

  if any(strcmp(cmp,wc))  % Wildcard only
     ok = ones(si);
     return
  end

end

%***************************************************************
% Compare

ok = zeros(si);

for ii = 1 : prod(si)
 
    s2 = size(str{ii},2);

    for jj = 1 : nc

        ok(ii) = ( ( s2 == 0 ) & ( nn(jj,1) == 0 ) );
       
        if ok(ii)
           break
        end
       
        if ( s2 >= nn(jj,2) ) & ~( nn(jj,1) == 0 )

           ok(ii) = compare(str{ii},cmp{jj},cc{jj});

           if ok(ii)
              break
           end
            
        end

    end

end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = compare( str , cmp , cc )

sc = size(cmp,2);

ok = 1;

for ii = 1 : size(cc,1) 

    s2 = size(str,2);

    ok = ( ok &  ( s2 >= cc(ii,3) ) );

    if ~ok
       break
    end

    ic  = ( cc(ii,1) : cc(ii,2) );

    i01 = ( cc(ii,[1 2]) == [ 1  sc ] );

    i0  = ( 1 : cc(ii,3) );
    i1  = s2 - cc(ii,3) + i0;

    i2  = cat( 2 , findstr( str , cmp(ic) ) , 0 );
    i2 = i2(1);
    i2 = ( i2 + cc(ii,3) - 1 ) * ~( i2 == 0 );

    i2  = s2       * (  i01(1) &  i01(2) & strcmp(str    ,cmp    ) ) + ...
          cc(ii,3) * (  i01(1) & ~i01(2) & strcmp(str(i0),cmp(ic)) ) + ...
          s2       * ( ~i01(1) &  i01(2) & strcmp(str(i1),cmp(ic)) ) + ...
          i2       * ( ~i01(1) & ~i01(2) );

    ok = ( ok & ~( i2 == 0 ) );

    if ~ok
       break
    end

    str = str( i2+1 : s2 );

end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,str] = chkcstr(str,opt)


% CHKCSTR  Checks Input for CellString, contains Strings !
%
%  [ok,str] = chkcstr(str,Option)
%
%  Option ~= 0 ==> CharacterArrays not allowed,
%
%   default: Option == 0   ==>  CharacterArrays --> CellString
%
 
if nargin < 2
   opt = 0;
end

if strcmp(class(str),'char') & isequal(opt,0)
   str = cellstr(str);
end

ok = iscellstr(str);
if ~ok
   return
end

try
  s = cat(2,str{:});
catch
  ok = 0;
  return
end
 
ok = ( strcmp(class(s),'char')  &  ( prod(size(s)) == size(s,2) ) );

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

