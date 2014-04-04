function [Msg,v0] = structcmp(v0,v1,exact,nr);

% STRUCTCMP  Compares 2 Structures
%
% [Msg,V] = STRUCTCMP( V0 , V1 , Exact )
%
%  Check if Fieldnames of V1 contained in V0 and
%    SubStructures too.
%  Checks Class DimensionLength and min. Size of Values. 
%
%  Works only for single Structures and SubStructures.
%
%  Returns ErrorMessages in Msg and the updated
%   Structure V = ( V0 with Values of V1 )
%
%  Exact == -1   No other FieldNames as in V0 are allowed
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

check = imag(exact);

check_ndim = ~( check < 0 );
check_size = ~( check > 0 );

check_size = ( check_size & check_ndim );

check = real(exact);

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

    elseif ( isequal(check,2)  |  isequal(check,-1) )

          Msg = [ Msg  nl(1:(end*(~isempty(Msg)))) bl ...
                  'Invalid FieldName "' f1{ii} '".' ];

    end

end

if ~( isequal(check,2)  | isequal(check,1) )  |  all(ok0)
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

