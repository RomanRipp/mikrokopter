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
% 'cellvec'   True of CellArray or Numeric Array of Vectors
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
case 'cellvec'

      [ok,v] = chkcvec(v,0);

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
            w = round(v);
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

%****************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,v,n] = chkcvec(v,opt)

% CHKCVEC  Checks Input for CellArray of numeric Vectors
%
%  [ok,C,N] = chkcvec(V,Option)
%
%  V CellArray of Vectors or 2-dimensional numeric Array
%
%  Option ~= 0 ==> Numeric Arrays not allowed,
%
%   default: Option == 0   ==>  Numeric Arrays --> CellArray
%
%  C CellArray of Vectors (single Row or single Column)
%  N Numeric Array with a Vector of C in each Row
%
%
 
if nargin < 2
   opt = 0;
end

n = [];

if isnumeric(v)
   n = v;
   v = {};
   ok = ( isequal(opt,0) & ( ndims(v) == 2 ) );
   if ok & ~isempty(n)
      m = size(n,1);
      v = cell(m,1);
      for ii = 1 : m
          v(ii) = { n(ii,:) };
      end
   end
   return
end

ok = iscell(v);
if ~ok
   return
end

for d = [ 2  1 ]
    ok = 1;
    try
      n = cat(d,v{:});
    catch
      ok = 0;
    end
    if ok
       s = size(n); p = prod(s);
       ok = ( ( p == s(d) ) & isnumeric(n) );
    end
    if ok
       break
    end
end

if ~ok
    return
end

m = prod(size(v));

n = zeros( m , min(1,p) );

if isempty(n)
   return
end

for ii = 1 : m 
    w = v{ii};
    if isempty(w)
       n(ii,:) = NaN;
    else
       w = w(:)';
       s = size(n,2);
       t = size(w,2);
       if s < t
          n = cat( 2 , n , NaN * zeros(m,t-s) );
       end
       n(ii,:) = NaN;
       n(ii,1:t) = w;
    end
end 
