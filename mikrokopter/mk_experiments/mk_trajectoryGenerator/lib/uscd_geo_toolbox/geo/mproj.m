function [x,y,s] = mproj(fcn,x,y,c,p,r)

% MPROJ  Projection of Coordinates using MappingToolbox
%
% [ X , Y , MapStruct ] = MPROJ( FCN , Lon , Lat , Origin , [Parallels] )
%
% FCN    = Function to call of MappingToolbox/MAPPROJ
%          specials: 'azma'  --> 'eqaazim'
%                    'azmd'  --> 'eqdazim'
%                    'cona'  --> 'eqaconic'
%                    'cond'  --> 'eqdconic'
%                    'cyla'  --> 'eqacylin'
%                    'cyld'  --> 'eqdcylin'
%                    'lamb'  --> 'lambert'
%                    'merc'  --> 'mercator'
%
% Origin = [ OriginLon OriginLat [Rotation] ]
%
% Use the Imaginary Unit "i" as last Input for reverse Projection:
%
% [ Lon , Lat ] = MPROJ( FCN , X , Y , Origin , [Parallels] , i )
%
%
% MapStruct = 
%    projection    :  FCN
%    geoid         : [ 6378137 0.08181919131097 ]  % WGS84
%    aspect        : 'normal'
%    zone          : '31N'
%    origin        : Origin([2 1 [3]])
%    mapparallels  : Parallels
%    flatlimit     : [ -90  90]
%    flonlimit     : [-Inf Inf]
%    angleunits    : 'degrees'
%    scalefactor   : 1
%    falseeasting  : 0
%    falsenorthing : 0
%
% Instead of a FunctionString (FCN) the MapStruct can used as 1. Input.
%
%
% MPROJ calls:
%
% [ X , Y ] = FEVAL( FCN , MapStruct , Lat , Lon , 'none' , 'forward' )
%
% reverse Projection:
%
% [Lat,Lon] = FEVAL( FCN , MapStruct ,  X  ,  Y  , 'none' , 'inverse' )
%
%


Nin  = nargin;
Nout = nargout;

s = struct( ...
'projection'    , { '' }        , ...
'geoid'         , { [ 6378137 0.08181919131097 ] } , ...
'aspect'        , { 'normal' }   , ...
'zone'          , { '31N'    }   , ...
'origin'        , { zeros(1,2) } , ...
'mapparallels'  , { zeros(1,1) } , ...
'flatlimit'     , { [ -90  90] } , ...
'flonlimit'     , { [-Inf Inf] } , ...
'maplatlimit'   , { [] } , ...
'maplonlimit'   , { [] } , ...
'angleunits'    , { 'degrees' }  , ...
'scalefactor'   , { 1 } , ...
'falseeasting'  , { 0 } , ...
'falsenorthing' , { 0 }          );

%*********************************************************
% Check Inputs

if Nin < 1
   error('Not enough InputArguments.');
end

if Nin == 2
   error('Input Y is missing.');
end
   
%-----------------------------------------------------

r = 0;

if Nin < 4
   c = [];
elseif isequal(c,i) & ( Nin == 4 )
   r = 1; 
   c = [];
end

if Nin < 5
   p = [];
elseif isequal(p,i)
   r = 1; 
   p = [];
elseif Nin == 6
   r = isequal(r,i);
end

%-----------------------------------------------------

msg = cell(0,1);

%-----------------------------------------------------
% FCN | MapStruct

si = size(fcn);

is_fcn = (   ischar(fcn) & ~isempty(fcn) & ( prod(si) == si(2) ) );
is_str = ( isstruct(fcn) & ( prod(si) == 1 ) );

if is_fcn

   fcn = lower(fcn);

   if si(2) == 4 
      switch fcn
       case 'lamb'
         fcn = 'lambert';
       case 'merc'
         fcn = 'mercator'
       otherwise
         suf = '';
         switch fcn(1:3)
          case 'azm'
                suf = 'azim';       
          case 'con'
                suf = 'conic';
          case 'cyl'
                suf = 'cylin';
        end
        if ~isempty(suf)
            fcn = sprintf('eq%s%s',fcn(4),suf);
        end
      end
   end

   m = '';

   if exist(fcn,'file') == 2

      s.projection = fcn;

      % Get Defaults

      try
         d = feval(fcn,s);
      catch
         m = sprintf('Error call FEVAL(%s,MapStruct):\n%s',upper(fcn),lasterr);
      end

      if isempty(m)
         [m,s] = structcmp( s , d , 1-i );
         if ~isempty(m)
             m = sprintf('Invalid MapStructure by %s.\n%s',upper(fcn),m);
         end
      end
          
   else
      m = sprintf('Function "%s" not found.',fcn);
   end

   if ~isempty(m)
       msg = cat( 1 , msg , {m} );
   end

elseif is_str

    
    [m,s] = structcmp( s , fcn , 1 );

    if ~isempty(m)
        msg = cat(1,msg,{sprintf('Invalid MapStructure.\n%s',m)});
    end

else

    msg = cat(1,msg,{'First Input must be a String or a MapStructure.'});

end

%-----------------------------------------------------
% Equal Size of X & Y

if ~isequal(size(x),size(y))
    msg = cat(1,msg,{'Size of X and Y must be agree.'});
end

%-----------------------------------------------------
% Origin

if ~isempty(c)
    sc = prod(size(c));
    if ~( isnumeric(c) & ( sc >= 2 ) )
       msg = cat(1,msg,{'Origin be numeric with min. 2 Elements: [ Lon Lat ].'});
    else
       c = c( 1 : min(3,sc) );
       c = c(:)';
       if sc == 2
          c = [ c  0 ];
       end
       c([1 2]) = c([2 1]);     % [ Lat  Lon ]
       s.origin = c;
    end
elseif is_str
    c  = s.origin;
    sc = prod(size(c));
    c  = c( 1 : min(3,sc) );
    c  = c(:)';
    if sc == 2
       c = [ c  0 ];
    end
    s.origin = c;
elseif is_fcn & ~r
    if ~isempty(y)
        s.origin(1) = meannan(y(:));
    end
    if ~isempty(x)
        s.origin(2) = meannan(x(:));
    end
end

%-----------------------------------------------------
% Parallel

if isempty(p)
   if ~is_str
       s.mapparallels = s.origin(1);
   end
else
   s.mapparallels = p;
end
 
%-----------------------------------------------------

if ~isempty(msg)
    error(sprintf('%s\n',msg{:}));
end

%*********************************************************
% Zone for UTM

if isequal(fcn,'utm')

   xm = s.origin(2);
   ym = s.origin(1);

   ok = all(isfinite([xm ym]));
   if ok
      ok = ( ( -80 <= ym ) &( ym <= 84 ) );
   end

   if ok

      xm = xm - 360 * floor(xm/360);
      xm = xm - 360 * ( xm > 180 );

      s.zone  = utmzone(ym,xm);
      s.geoid = utmgeoid(s.zone);

      [yl,xl] = utmzone(s.zone);

      s.origin([1 2]) = [ 0  mean(xl) ];

      s.falseeasting  = 500000;

   end

end

%*********************************************************
% SavedPoints

spt = struct( 'trimmed' , {[]} , ...
              'clipped' , {[]}        );

%*********************************************************
% Call Function

fcn = s.projection;

try

  if r

     [y,x] = feval(fcn,s,x,y,'none','inverse',spt);

  else

     [x,y] = feval(fcn,s,y,x,'none','forward',spt);

  end

catch

  error(sprintf('Error call FEVAL(%s):\n%s',upper(fcn),lasterr));

end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function x = meannan(x,dim)

% MEANNAN    Column mean with missing data.
%
%  Y = MEANNAN(X) returns the mean of each column of X as a row vector
%  where missing data values are encoded as NaNs. For vectors, MEANNAN(X)
%  returns the mean value of the elements in X.
%  Y = MEANNAN(X,DIM) averages over dimension DIM (only version 5 or higher)

%  This code has been suggested by Douglas M. Schwarz (schwarz@kodak.com) 
%  in the news group comp.soft-sys.matlab.
%
%  C. Mertens, IfM Kiel
%  $Revision: 1.1 $ $Date: 1995/03/08 14:27:31 $
%
% added ~isempty	G.Krahmann, IfM Kiel, Oct 1995
% added compatibility for version 5.	G.Krahmann, LODYC Paris, Jul 1997

% check for version

 
v=version;
if any(strcmp(v(1),{ '3' '4'}))

    dim = 0;

else

  if nargin==1
    dim = find( size(x) > 1 );
    if isempty(dim)
      dim = 1;
    else
      dim = dim(1);
    end
   end

end



      ok   = ~isnan(x);

  if dim == 0

      bad  = find( ~ok );
    x(bad) = 0*bad; 

    ok     = sum(ok);
    x      = sum(x) ./ ( ok + ( ok == 0 ) );

      bad  = find(ok==0);
    x(bad) = NaN*bad;

  else

    x(find(ok==0)) = 0;
 
    ok     = sum(ok,dim);
    x      = sum(x,dim) ./ ( ok + ( ok == 0 ) );

    x(find(ok==0)) = NaN;

  end


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,v0] = structcmp(v0,v1,exact,nr);

% STRUCTCMP  Compares 2 Structures
%
% [Msg,V] = STRUCTCMP( V0 , V1 , Exact )
%
%  Check if Fieldnames of V1 contained in V0 and
%    SubStructures too.
%  Checks Class DimensionLength and min. Size of Values. 
%
%  Returns ErrorMessages in Msg and the updated
%   Structure V = ( V0 with Values of V1 )
%
%  Exact == -1   No other FieldNames of V1 are allowed
%  Exact ==  1   All FieldNames of V0 must be matched
%  Exact ==  2   FieldNames of V0 and V1 must be equal
%
%  Exact = Exact + i   No Check of min. Size
%  Exact = Exact - i   No Check of min. Size and DimensionLength
%
%  default:  Exact == 0 only existing fields will checked
%
%  see also:  CLASSCHK, CHECK_CNF
%

Msg = '';

if nargin < 2
  Msg = 'Not enough Input Arguments.';
end

if nargin < 3
   exact = 0;
end

if nargin < 4
   nr = 0;
end

nl = char(10);

bl = char( 32 * ones(1,3*nr) );

check_ndim = ~( imag(exact) < 0 );
check_size = ~( imag(exact) > 0 );

check_size = ( check_size & check_ndim );

%----------------------------------------------------
% Check Class

c0 = class(v0);

[ok,v1] = classchk(v1,c0);

if ~ok & isempty(v1)
    [v1,m] = defval(c0,[0 0],v0);
    ok = isempty(m);
end

if ~ok
   Msg = [ bl  'Variable must be of Class "'  c0  '".' ];
   return
end

%----------------------------------------------------
% Check DimensionLength and Size

if check_ndim

   n0 = ndims(v0);
   n1 = ndims(v1);

   ok = ( n0 == n1 );

   if ~ok
      nd  = sprintf('%.0f',n0);
      Msg = [ bl 'Variable must have DimensionLength ' nd '.' ];
      return
   end

   if check_size

      s0 = size(v0);
      s1 = size(v1);

      ok = ( all( s0 <= s1 )  |  strcmp(c0,'char') );

      if ~ok

         si = sprintf('%.0f x ',s0);
         si = si( 1 : end-3 );

         Msg = [ bl 'Variable must have Size: [ ' si ' ].' ];

         return

      end

    end

end

%----------------------------------------------------
   
if ~strcmp(c0,'struct')
   v0 = v1;
   return
end

if ~( ( prod(size(v0)) == 1 ) & ( prod(size(v1)) == 1 ) )
    Msg = [ bl 'Structures must be single.' ];
    return
end

%----------------------------------------------------

f0 = fieldnames(v0);
f1 = fieldnames(v1);

n0 = prod(size(f0));
n1 = prod(size(f1));

ok0 = zeros(n0,1);
ok1 = zeros(n1,1);

for ii = 1 : n1

    ok1(ii) = any( strcmp( f0 , f1{ii} ) );

    if ok1(ii)

       jj = find( strcmp( f0 , f1{ii} ) );

       jj = jj(1);

       ok0(jj) = 1;

       [Msg1,p] = structcmp( getfield(v0,f0{jj}) , getfield(v1,f1{ii}) , ...
                         exact , nr+1 );

       if isempty(Msg1)
 
           v0 = setfield( v0 , f0{jj} , p );

       else

          Msg = [ Msg  nl(1:(end*(~isempty(Msg)))) bl ...
                  'Invalid SubStructure for "' f1{ii} '".' nl  Msg1 ];

       end

    elseif ( isequal(exact,2)  |  isequal(exact,-1) )

          Msg = [ Msg  nl(1:(end*(~isempty(Msg)))) bl ...
                  'Invalid FieldName "' f1{ii} '".' ];

    end

end

exact = real(exact);

if ~( isequal(exact,2)  | isequal(exact,1) )  |  all(ok0)
   return
end

jj = find(~ok0);

if ~isempty(jj)

  app = 's';
  app = app( 1 : ( end * ( prod(size(jj)) > 1 ) ) );

  Msg = [ Msg nl(1:(end*(~isempty(Msg)))) bl ...
          'Didn''t found Field'  app  ': "'  strhcat(f0(jj),'", "') '".' ];

end
 
%****************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,v] = classchk(v,cl);

% CLASSCHK  Check if Input is of given Class
%
% [ Ok , V ] = CLASSCHK( V , Class )
%
% Returns ONE, if V is of given Class or Class is empty.
%
% Special Class's are:
%
% 'numeric'   True if ISNUMERIC
% 'string'    True for ONE-Row CharacterArray or ONE-Element CellArray of it
% 'cellstr'   True of CellArray or CharacterArray of Strings
% 'object'    True if ISOBJECT
%
% Special Handlings for Matlab-defined Class's:
%
% 'double'    True for ISNUMERIC
% 'single'    True for ISNUMERIC
% 'sparse'    True for ISNUMERIC
% 'logical'   True for ISNUMERIC  with Values [ 0 | 1 ]
% 'int<B>'    True for ISNUMERIC  with Values [ -2^(B-1) .. 2^(B-1)-1 ]
% 'uint<B>'   True for ISNUMERIC  with Values [  0       .. 2^(B)-1   ]
% 'char'      True for ISCELLSTR
% 'cell'      True for CharacterArray
% 'struct'    True for Object
% ClassName   True if  Object is of given Class
% @ClassName  True if  Object inherits of given Class, using ISA
%



Nin  = nargin;
Nout = nargout;

if Nin < 1
   ok = 0;
   v  = ones(0,0);
   return
end

ok = 1;

if Nin < 2
   cl = '';
elseif  ~chkstr(cl,0)
   error('Class must be a String.')
end

%****************************************************

if isempty(cl)
   return
end

cv = class(v);

ok = isequal(cl,cv);

if ok
   return
end

%****************************************************

out = ( Nout > 1 );

switch cl

%----------------------------------------------------
case 'numeric'

      ok = isnumeric(v);

      cl = '';

%----------------------------------------------------
case 'string'

      ok = chkstr(v,0);

      if ~ok & strcmp(cv,'cell') & ( prod(size(v)) == 1 )
          ok = chkstr(v{1},0);
          if ok
             v = v{1};
          end
      end

      cl = '';

%----------------------------------------------------
case 'cellstr'

      [ok,v] = chkcstr(v,0);

      cl = '';

%----------------------------------------------------
case 'object'

      ok = isobject(v);

      cl = '';

%----------------------------------------------------
case { 'int8' 'int16' 'int32' 'uint8' 'uint16' 'uint32' }

      ok = isnumeric(v);

      if ok & ( ~strcmp(cv,'logical') | out )

         if strcmp(cv,'double')
            w = v;
         else
            w = double(v);
         end

         if ~strcmp(cv,'logical')

            s  = ~( cl(1) == 'u' );               % TRUE for Signed Integer

            b = eval( cl( (5-s) : end ) ) - s;    % Bytes

            l = 2^b * [ (0-s)  1 ] + [ 0  -1 ];   % Intervall

            ok = all( ( l(1) <= w(:) ) & ( w(:) <= l(2) ) );

         end

         if ok & out
            v = w;
         end

      end

%----------------------------------------------------
case { 'double'  'single'  'logical'  'sparse' }

      ok = ( isnumeric(v) | strcmp(cv,'logical') );

      if ok & ( strcmp(cl,'logical') | out )

         if strcmp(cv,'double')
            w = v;
         else
            w = double(v);
         end

         if strcmp(cl,'logical')
            ok = ( all( w(:) == 0 ) | all( w(:) == 1 ) );
         end

         if ok & out
            v = w;
         end

      end

%----------------------------------------------------
case 'char'

      ok = iscellstr(v);

%----------------------------------------------------
case 'cell'

      ok = strcmp(cv,'char');

      cl = 'cellstr';

%----------------------------------------------------
case 'struct'

      ok = isobject(v);
  
%----------------------------------------------------
otherwise

      ok = strcmp(cl(1),'@');

      if ok
         cl  = cl(2:end);
         ok  = isa(v,cl);
      end

end

%----------------------------------------------------

if ~( ok & out ) | isempty(cl)
   return
end

%----------------------------------------------------
% Convert Value

try
  v = feval(cl,v);
end

%****************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [val,msg] = defval(typ,si,val);

% DEFVAL   Returns DefaultValue for Type and Size
%
% [val,msg] = DEFVAL(typ,si,val);
%
% The 3. Input "val" is used, if the "typ" 'struct' is requested,
%   to create a Structure with new Size from existing Structure "val".
%
% typ: 'numeric'
%      'cellstr'
%      'string'
%      'object'
%      <class>
%      ''
%

Nin = nargin;

if Nin < 1
   typ = '';
end

if Nin < 2
   si = [ 0  0 ];
end

if Nin < 3
   val = [];
end


msg = '';


si(find(isnan(si))) = 0;
 

switch typ

  case { 'char'  'string' }
 
          val = char( 32 * ones(si) );

  case { 'double'  'single'  'numeric'  '' }

          val = NaN * ones(si);

  case {  'int8'     'int16'    'int32' ...
         'uint8'    'uint16'   'uint32' ...
         'logical'  'sparse'        }

          val = feval( typ ,  zeros(si) );

  case 'cell'

          val = cell(si);

  case 'cellstr'

          val = cell(si);

          if ~any( si == 0 )
              val(:) = { '' };
          end

  case 'struct'

       if ~strcmp(class(val),'struct');

           val = struct( 'NaN' , { [] } );

       end

       ns = prod(size(si));

      ind = cell(ns,1);

      for ii = 1 : ns
          ind{ii} = ones(1,si(ii));
      end
 
      val = val(ind{:});

  otherwise

     try

        val = feval(typ,si);

     catch

        msg = lasterr;

     end

end

%****************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%****************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

