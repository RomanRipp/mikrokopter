function [Msg,MsgC,cnf,ok,cfg] = check_cnf(cnf,ini,varargin);

% CHECK_CNF   Checks ConfigurationStructure
%
% [Msg,MsgCNF,CNF,INI_OK,V] = CHECK_CNF( CNF , INI , V )
%
% INI = { FieldName  [Type]  [Size]  [Optional] }
%
% Checks, if all FieldNames defined in 1. Column of INI 
%  exist in the Structure CNF.
%
% A StructureValue which match the FieldNames is Ok,
%   if it match the Type and the Size or any of the OptionalValues.
%  
%  Type = { 'none'  'numeric'  'string'  
%           'cellstr'  'cellvec'  'object'   <ClassName> }
%
%  Optional =  CellArray of Values, which can be equal to the StructureValue
%                In case of a single element structure as element of Optional
%                a single StructureValue is checked for equal fieldnames only.
%
%  Specials for empty Type or 'none':
%
%        empty Optional: any StructureValue matching Size ok
%     nonempty Optional: only OptValues ok
%
%  The DimensionNumber of the StructureValues must be equal
%    then the the Numbers of Values in Size.
%  The DimensionLength of the StructureValues must be equal or larger
%    then the Value in Size, larger DimensionLengths will cutted.
%  Single StructureValues will multiplied to Size.
%
%  Specials for empty Size, negative or NaN-Values in Size:
%
%            empty Size:  any Size Ok
%         negative Value: all Lengths smaller of absolute Value  Ok 
%              NaN-Value: all Lengths in this Dimensions Ok
%
%  Special Optional Values for ColorNames: '#colspec#' | '#colspecX#'
%
%   #colspec#    Checks for Matlab's ColorSpecifier
%   #colspecX#   Checks for Matlab's ColorSpecifier or XRGB-Color
%                 XRGB required!
%
%--------------------------------------------------------------------------
%
% If an StructureInput V is given, the Value of CNF will replaced
%    with the Values of V before check.
%
%  V can be a Structure, a CellArray, defines a Structure using struct(V{:}),
%    or multiple FieldName-Value-Inputs.  
%
% INI_OK : -1  FieldName of INI doesn''t exist
%           0  Invalid Value
%           1  Value ok
%
%--------------------------------------------------------------------------
%
% see also:   CLASSCHK, DEFVAL, STRUCTCHK, STRUCTCMP, COLSPEC, XRGB
%

Msg  = '';
MsgC = '';

Msg0 = 'CHECK_CNF: ';

Nin  = nargin;

clspc = { '#colspec#' '#colspecX#' };  % Special Optional for COLSPEC

%-----------------------------------------------

if Nin < 2
  Msg = [ Msg0  'Not enough Input Arguments.' ];
  return
end

if isempty(ini) | isempty(cnf)
   return
end

Msg = cell(0,1);

if ~( isstruct(cnf)  &  ( prod(size(cnf)) == 1 ) ) 
  Msg = cat( 1 , Msg , {'CNF must be a Structure with Length 1.'} );
end

%--------------------------------------------------------------------
% INI

s2 = size(ini,2);

ok = ( iscell(ini) & ~isempty(ini) & ( ndims(ini) == 2 ) );

if ~ok

   Msg = cat( 1 ,  Msg , {'INI must be a CellArray: { FieldName Type [Size] [Optional] }.'} );

else

   s2 = size(ini,2);

   if s2 > 1

      if     s2 > 4

             ini = ini(:,1:4);

      elseif s2 < 4

             ini = ini(:,[ (1:s2) ones(1,4-s2) ]);

             if s2 < 2, ini(:,3) = { 'none' };    end  % Type
             if s2 < 3, ini(:,3) = { [] };        end  % Size
             if s2 < 4, ini(:,4) = { cell(0,1) }; end  % Optional

      end

   end

   if ~chkcstr( ini(:,1:min(s2,2)) );
       Msg = cat( 1 ,  Msg , {'Elements of 1. and 2. Column of INI must be Strings (FieldNames).'} );
   elseif any(strcmp(ini(:,1),''));
       Msg = cat( 1 ,  Msg , {'Elements of 1. Column of INI must be nonempty Strings (FieldNames).'} );
   end

   try
      v  = cat(2,ini{:,3});
   catch
      v = zeros(2,1);
   end

   if ~( isnumeric(v) & ( ( prod(size(v)) == size(v,2) ) | isempty(v) ) );
       Msg = cat( 1 ,  Msg , {'Elements of 3. Column of INI (Size) must be RowVectors.'} );
   elseif ~isempty(v)
       if ~all( ( mod(v,1) == 0 ) | isnan(v) )
           Msg = cat( 1 ,  Msg , {'Elements of 3. Column of INI (Size) must be Integers.'} );
       end
   end

end

%--------------------------------------------------------------------
% Structure

if Nin < 3

   cfg = [];

else

    % Check for CellArrays
    if Nin > 3
       for ii = 2 : 2 : Nin-2
           if iscell(varargin{ii}) & ~( prod(size(varargin{ii})) == 1 )
              varargin{ii} = {varargin{ii}};
           end
       end
    end

   [MsgV,cfg] = structchk(varargin);

   if ~( ( isempty(MsgV) & ( prod(size(cfg)) == 1 ) ) | isempty(cfg) )
         Msg = cat( 1 , Msg , ...
                    {'3. Input must define a ONE-Element Structure.'} , {MsgV} );
   end

end

%--------------------------------------------------------------------

if ~isempty(Msg)
    Msg = sprintf('%s\n',Msg{:});
    Msg = sprintf('%sInvalid Inputs.\n%s',Msg0,Msg);
    ok  = NaN * zeros(size(ini,1),1);
    return
end

Msg = '';

%*******************************************************
% Use Values of Structure

field = fieldnames(cnf);

lfld = lower(field);

if ~isempty(cfg)
    fld = fieldnames(cfg);
    for ff = fld(:)'
        ok = strcmp(field,ff{1});
        if ~any(ok)
            ok = strcmp(lfld,lower(ff{1}));
        end
        if any(ok)
           ok  = find(ok);
           ok  = ok(1);
           cnf = setfield( cnf , field{ok} , getfield(cfg,ff{1}) );
        end
    end
end

%*******************************************************

nl =  sprintf('\n      ');
          
ni = size(ini,1);

%------------------------------------------------
% Check for FieldNames of CNF

ok = zeros(ni,1);

for ii = 1 : ni
    ok(ii) = any( strcmp( field , ini{ii,1} ) );
end

if all(ok)
    MsgF = '';
else
    MsgF = sprintf( 'Missing FieldNames:%s%s' , nl , ...
                    strhcat( ini(find(~ok),1) , ', ' , 6 , nl ) );
end

ok = ok - 1 * ( ok == 0 );

if ~any( ok == 1 )
    MsgC = MsgF;
    return
elseif size(ini,2) < 2
    return
end

fok = find( ok == 1 );

%------------------------------------------------

vok     = ones(ni,1);

MsgC    = cell(ni,1);
MsgC(:) = {''};

for ii = fok(:)'

     m = cell(0,1);

  %-----------------------------------------------

  fld = ini{ii,1};
  typ = ini{ii,2};
  siz = ini{ii,3};
  opt = ini{ii,4};

  if ~iscell(opt)
      opt = { opt };
  end

  tn = ( strcmp(typ,'none') | isempty(typ) ); % Any Type allowed

  o0 = isempty(opt);

  o1 = ( tn & ~o0 );           % Only OptValues allowed

  s0 = isempty(siz);           % Any Size allowed

  s0 = ( s0 | o1 );

  sn = ~s0;

  if sn
     isn =  isnan(siz);
     s0  =    all(isn);
     if ~s0 & any(isn)
         jj  = size(siz,2) - sum(cumprod(double(isn(end:-1:1))));
         jj  = max(2,jj);
         siz = siz(1:jj);
         isn = isn(1:jj);
     end
     sn = any( isn | ( siz < 0 ) );  % Any Size allowed if empty VAL
  end


  val = getfield( cnf , fld );
  svl = size(val);

  v0 = isempty(val);           % Any Type allowed if NaN's in SIZ

  %-----------------------------------------------
  if sn & v0   % DefaultValue
  %-----------------------------------------------

      siz(find(isn)) = 0;

      if tn
         val = zeros(max(siz,0));
      else
         val = defval(typ,max(siz,0),val);
      end

  %-----------------------------------------------
  elseif ~( tn & s0 )   % Check Type & Size
  %-----------------------------------------------

      %------------------------------------------
      % Check Type

      if ~tn
          [tok,val] = classchk(val,typ);
          if ~tok
              m = cat( 1 , m , {sprintf('%sValue must be of Type "%s".',nl,typ)} );
          else
              svl = size(val);
          end
      end

      %------------------------------------------
      % Check Size

      if ~s0

          % Check for Vector

          sok = ( ( max(siz) == 1 ) & ( sum(isn) == 1 ) );

          sok = ( sok & ( prod(svl) == max(svl) ) );

          if sok
            
             val = val(:);

             jj = find(isn);

             val = shiftdim(val,-(jj-1));

          else

             nd = size(siz,2);

             sok = ( ndims(val) == nd );

             chk = sok;
             if chk
                chk = ( ( siz < 0 )  &  ( svl <= abs(siz) ) );
                chk = ~all( chk | isn | ( svl == siz ) );
             end
   
             if chk

                ind = cell(nd,1);
                ind(:) = {':'};

                for jj = find(~isn)
                    sv = svl(jj);
                    sz = siz(jj); 
                    s1 = ( sv == 1 );
                    sok = (  sok  &  ( s1 | ( sv >= sz ) )  );
                    if sz < 0
                       ind{jj} = ( 1 : min(sv,abs(sz)) );
                    elseif s1
                       ind{jj} = ones(1,sz);
                    else
                       ind{jj} = ( 1 : sz );
                    end
                end  

                if sok
                   val = val(ind{:});
                end
    
             end

          end 

          if ~sok

              siz = sprintf('%.0f x ',siz);
              siz = siz( 1 : end-3 );

              m = cat( 1 , m , {sprintf('%sValue must have Size: [ %s ].',nl,siz)} );

          end

      end

  %-----------------------------------------------
  end
  %-----------------------------------------------

  %-----------------------------------------------
  % Check OptValues

  if ( ~isempty(m) & ~o0 ) | o1

     if o1
        m = 'fail';
     end

     val = getfield(cnf,fld);

     for vv = opt(:)'

         vs = chkstr(vv{1},1);
         vc = vs;
         if vc
            vc = any(strcmp(vv{1},clspc));
         end
 
         if vc

            try
                cl = colspec(val);
                if isempty(cl) & any(vv{1}=='X')
                   cl = colspec(val,1);
                end
                if ~isempty(cl) & ~tn
                    [cok,cl] = classchk(cl,typ);
                    if cok
                       val = cl;
                    end
                end
            catch
                cl = [];
            end

            kk = ~isempty(cl);

         else

             kk = isequal(val,vv{1});
             if ~kk & chkstr(val,1) & vs
                 kk = isequal(lower(val),lower(vv{1}));
                 if kk
                    val = vv{1};
                 end
             elseif isstruct(val)   & ( prod(size(val))   == 1 ) & ...
                    isstruct(vv{1}) & ( prod(size(vv{1})) == 1 )
                 kk = isempty(structcmp(vv{1},val,1-i));
             end

         end

         if kk
            m = '';
            break
         end

     end

     if o1 & ~isempty(m)
        m = {sprintf('%sValue must be any of:%s%s.',nl,nl,var2mstr(opt,-1))};
     end

  end

  %-----------------------------------------------

  vok(ii) = isempty(m);

  if vok(ii)
      cnf = setfield( cnf , fld , val );
  else
      m = sprintf('%s',m{:});
      MsgC{ii} = sprintf('Invalid Data for Field "%s".%s',fld,m);
  end

end
% ii 

ok(fok) = vok;

ii = ~vok;

if any(ii)
   ii = find(ii);
   MsgC = sprintf('%s\n',MsgC{ii});
else
   MsgC = '';
end

if ~isempty(MsgF)
    MsgC = sprintf('%s\n%s',MsgF,MsgC);
end

%****************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [Msg,V] = structchk(V,txt);

% STRUCTCHK( V )   Checks if V is an StructArray,
%
%  or CellArray with Field-Value-Pairs and build it.
%
%  [ Msg , Structure ] = STRUCTCHK( V , ActionText )
%
%  returns the Error's in Msg
%    or the Stucture build from V
%

Msg = '';

if isstruct(V)
   return
end

if nargin == 2
   if ischar(txt) & ~isempty(txt)
      txt = sprintf(' after %s-Action',txt);
   else
      txt = '';
   end
else
   txt = '';
end


n = prod(size(V));

%---------------------------------------------
if     n == 0
%---------------------------------------------

   Msg = sprintf('Empty Input%s.',txt);

%---------------------------------------------
elseif n == 1
%---------------------------------------------

  % Single Input must be StructArray

  if iscell(V)
     V = V{1};
  end

  if ~isstruct(V)
    Msg = sprintf('Single Input%s must be a StructArray.',txt);
  end


%---------------------------------------------
else
%---------------------------------------------

  % FieldName-Value-Pairs ==>  build StructArray

  if mod(n,2) ~= 0

     Msg = sprintf('Number of Multiple Input%s must be even (Field-Value-Pairs).',txt);

  elseif ~chkcstr(V(1:2:(n-1)));  % FieldNames

      Msg = sprintf('Multiple Input%s must be FieldName-Value-Pairs.',txt);

  else

     V = reshape(V,2,n/2);

      % Check for duplicate FieldNames, take last

      [c,si] = sort(V(1,:));
      c = double(char(c));
      c = ( sum( ( diff(c,1,1) == 0 ) , 2 ) == size(c,2) );
      if any(c)
         c = find(c);     % Leave last
         V(:,si(c)) = [];
      end

     try 
        V = struct(V{:});
     catch
        Msg = sprintf('Error call STRUCT from Input%s.\n%s',txt,lasterr);
     end

  end

end

%****************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
%      'cellvec'
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

  case 'cellvec'

          val = cell(si);

          if ~any( si == 0 )
              val(:) = { [] };
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

function str = var2mstr(val,bl)

% VAR2MSTR   Converts Variables to String with Matlab expression
%
% String = VAR2MSTR( V )
%
%  Converts 2-dimensional MatlabVariables 
%   to their String-Expression, using by EVAL 
%             
% Variables could be Numeric, Char, Cell- or StructArrays.
%
% String = VAR2MSTR( V , N )
%
%  real(N) = Number of Blanks around the Seperator
%  imag(N) = Number of Blanks around a numeric
%
%  real(N) < 0  ==> No "," - Separator along 2. Dimension
%
% default: N = 1
%
 
if nargin < 2
   bl = [];
end

if isempty(bl)
   bl = 1 + 0*i;
end

sp =     real(bl(1));   % Blank arround Seperator
nb = abs(imag(bl(1)));  % Blank arround Numbers

nb = char(32*ones(1,nb));

   % Usage to convert Number x to char:
   % vk = floor(log(abs(x))/log(10))+1;
   % vk = vk*(vk>0);
   % sprintf(sprintf(kform,nk+vk),x);

 % Seperator for Dimensions
   sep0 = char(32*ones(1,abs(sp)));

   sep1 = [ sep0  ';'  sep0 ];

   sep0 = char(32*ones(1,abs(sp)));
   if sp < 0
      sep2 = sep0;
   else
      sep2 = [ sep0  ','  sep0 ];
   end

 % Make 2-Dimensional
   if isobject(val)
      si = size( struct(val) );
   else
      si = size(val);
   end

  val = reshape(val,si(1),prod(si(2:end)));

   si = size(val);
   si = si * (~isempty(val));

try

   str = get_str(val,si,sep0,sep1,sep2,nb,bl);

catch

   str = -1;

end

if isequal(str,-1)

   % Build String from Size and class

   cl = class(val);

   ff = [ '[]' ; '{}' ];

   ff = ff(1+strcmp(cl,'cell'),:);

   str = sprintf('%.0fx',si);

   str = cat(2,ff(1),str(1:end-1),' ',cl,ff(2));

end

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = get_str(val,si,sep0,sep1,sep2,nb,bl)

 str = -1;

 %**************************************************
 % Char
 if ischar(val) 

   if si(1) <= 1
      str = [ ''''  val '''' ];
   else
      str = [ '['  sep0 ]; 
      for ii = 1 : si(1)
        str = [ str '''' val(ii,:) '''' ...
                sep1(1:(end*(ii~=si(1))))  ];
      end
      str = [ str  sep0  ']' ];
   end 

 %**************************************************
 % Numeric
 elseif isnumeric(val)

   if ~strcmp(class(val),'double')
       val = double(val);
   end

   clear eps
   nk    = floor(abs(log(eps)/log(10)))-1 ;
   kform = [ nb  '%%.%.0fg'  nb ];

   if prod(si) == 0
      str = '[]';
   elseif prod(si) == 1
      str = numstr(val,nb,nk,kform);
   else
      str = [ '['  sep0 ]; 
      for ii = 1 : si(1)
        for jj = 1 : si(2)
            str = [ str  numstr(val(ii,jj),nb,nk,kform)  sep2(1:(end*(jj~=si(2)))) ];
        end
        str = [ str  sep1(1:(end*(ii~=si(1))))  ];
      end
      str = [ str  sep0  ']' ];
   end 

 %**************************************************
 % CellArray
 elseif iscell(val)

   if prod(si) == 0

         str = '{}';

   else

      str = [ '{'  sep0 ] ; 
      for ii = 1 : si(1)
        for jj = 1 : si(2)
          str = [ str  var2mstr(val{ii,jj},bl)           ...
                     sep2(1:(end*(jj~=si(2))))     ];
        end
        str = [ str  sep1(1:(end*(ii~=si(1))))  ];
      end
      str = [ str  sep0  '}' ];
 
   end

 %**************************************************
 % StructArray
 elseif isstruct(val)

    fnames = fieldnames(val);
    fnames = fnames(:);
    nf     = size(fnames,1);
   
    fsep = ',';

    str = [ 'struct(' ];
  
    for ff = 1 : size(fnames,1);

      str = [ str  var2mstr(fnames{ff},bl) fsep '{' ]; 

      for ii = 1 : si(1)
        for jj = 1 : si(2)
          str = [ str  var2mstr(getfield(val,{ii,jj},fnames{ff})) ...
                     sep2(1:(end*(jj~=si(2))))     ];
        end
        str = [ str  sep1(1:(end*(ii~=si(1))))  ];
      end

      str = [ str '}' fsep(1:(end*(ff~=nf)))  ];

    end

    str = [ str ')' ];

 %**************************************************
 % StructArray
 elseif isobject(val)

    str = disp(val);

 end

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = numstr(val,nb,nk,kform)

vi = imag(val);

if ~( vi == 0 )

    vr = real(val);

    sr = '';
    if ~( vr == 0 )
        sr = numstr( vr  , '' , nk , kform );
    end

    add_plus = ( ~isempty(sr) & ( ( vi > 0 ) | isnan(vi) ) );

    if abs(vi) == 1
       si = 'i';
       if vi < 0
           si = [ '-' si ];
       elseif add_plus
           si = [ '+' si ];
       end
    else
       si = numstr( vi , '' , nk , kform );
       if add_plus
           si = [ '+' si ];
       end
       si = [ si '*i' ];
    end

    str = [ nb sr si nb ];

    return

end

if val == 0
   str = '0';
elseif isnan(val)
   str = 'NaN';
elseif isinf(val)
   str = 'Inf';
   if sign(val) == -1
      str = [ '-' str ];
   end
elseif abs(mod(val,pi)) < 1e3*eps
   val = round(val/pi);
   str = 'pi';
   if ~( val == 1 )
       if val == -1
          str = '-pi';
       else
          str = [ numstr(val,'',nk,kform) '*pi' ];
       end
   end
else
   vk = floor(log(abs(val))/log(10))+1;
   vk = vk*(vk>0);
  str = sprintf(sprintf(kform,nk+vk),val);
end
        
str = [ nb  str  nb ];

%****************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [n,c] = colspec(v,x)

% COLSPEC   Returns Matlab's or XRGB color specifier 
%
% Get List of Colors:
%
% [ Name , RGB ] = COLSPEC      List of ColorSpecifier
% 
% [ Name , RGB ] = COLSPEC(1)   List of ColorSpecifier and XRGB-Colors
% 
% Find a Color:
%
%   RGB  = COLSPEC( Name )       Color corresponds with Name, Matlab's only
%   RGB  = COLSPEC( Name , 1 )   Color corresponds with Name, incl. XRGB-ColorNames
%
%   Name = COLSPEC( RGB )        Find Name corresponds with Color, Matlab's only
%
%   Name = COLSPEC( RGB , 1+R*i) Find Name corresponds with Color, incl. XRGB
%
%   R is the ToleranceRadius to search for a XRGB-Color, default: 0.03
%
% see also: XRGB
%

Nin = nargin;

nv = 0;

if Nin == 0
   nv = 1;        % Not V
   v  = [];
   x  = 0;
elseif Nin < 2
   nv = ( isnumeric(v) & ( prod(size(v)) == 1 ) );
   if nv
      x = v;
   else
      x = 0;
   end
end

rad = imag(x);       % ToleranceRadius

x   = isequal(real(x),1);  % Use XRGB

%------------------------------------------------------

n = {   'r'        [ 1  0  0 ]
        'g'        [ 0  1  0 ]
        'b'        [ 0  0  1 ]
        'c'        [ 0  1  1 ]
        'm'        [ 1  0  1 ]
        'y'        [ 1  1  0 ]
        'k'        [ 0  0  0 ]
        'w'        [ 1  1  1 ]
        'red'      [ 1  0  0 ]
        'green'    [ 0  1  0 ]
        'blue'     [ 0  0  1 ]
        'cyan'     [ 0  1  1 ]
        'magenta'  [ 1  0  1 ]
        'yellow'   [ 1  1  0 ]
        'black'    [ 0  0  0 ]
        'white'    [ 1  1  1 ]   };

c = cat(1,n{:,2});
n = n(:,1);

%------------------------------------------------------

if nv
   m = size(n,1) / 2;         % Unique Colors
   c = c(1:m,:);
   n = n(1:m);
   if x
      try
         [nx,cx] = xrgb(1);      % Unique Colors
         n  = cat( 1 , n , nx );
         c  = cat( 1 , c , cx );
      catch
         warn(sprintf('Can''t get ColorList by XRGB.\n%s',lasterr))
      end
   end
   return
end

%------------------------------------------------------

s = size(v);
p = prod(s);

%----------------------------------------------
if isnumeric(v) & isequal(s,[1 3])
%----------------------------------------------

   d = c - v(ones(1,size(c,1)),:);
   d = all(d==0,2);

   if  any(d)

       d = find(d);
       d = d(1);
       n = n{d};
       c = c(d,:);

   elseif x

       try

          rad = rad + 0.03 * ( rad == 0 );

          [n,c] = xrgb(1);  % Unique Colors

           m = ones(size(c,1),1);

           d = ( c - v(m,:) );

           d = sum( d.^2 , 2 );

           ok = ( d <= rad^2 );

           if any(ok)
              nk =  sum(ok);
              ok = find(ok);
              if nk > 1
                 [d,ii] = min(d(ok));
                  ok = ok(ii);
              end
              n = n{ok};
              c = c(ok,:);
           else
              n = '';
              c = [];
           end

        catch

           warn(sprintf('Can''t get Color by XRGB.\n%s',lasterr))
           
           n = '';
           c = [];

        end

   else
       n = '';
       c = [];
   end

%----------------------------------------------
elseif ischar(v) & ~isempty(v) & ( p == s(2) )
%----------------------------------------------

   v = lower(v);

   for ii = 1 : size(n,1)

       nn = n{ii}( 1 : min(p,size(n{ii},2)) );

       ok = strcmp( nn , v );      

       if ok
          break
       end

   end


   if  ok
       n = n{ii};
       c = c(ii,:);
   elseif x
       try
          [c,h,n] = xrgb(v);
       catch
          warn(sprintf('Can''t get ColorList by XRGB.\n%s',lasterr))
          n = '';
          c = [];
       end
   else
       n = '';
       c = [];
   end

   v = c;
   c = n;
   n = v;

%----------------------------------------------
else
%----------------------------------------------

   error('Input must be an RGB-Tripel or String.');

%----------------------------------------------
end
%----------------------------------------------

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function warn(wrn)

% Display Warning-Message with Beep

ww = warnstat;

if strcmp(ww,'off') | isempty(wrn)
   return
end

warning('on');
  
fprintf(1,'\n%s',char(7));

warning(wrn);

fprintf(1,'\n%s',char(7));

warning(ww);


%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ww = warnstat

% WARNSTAT  Returns global WarningStatus
%
%  WARNSTAT returns the Status of WARNING
%
% Matlab R<13   WARNING
% Matlab R>12   WARNING for Identifier ALL
%

ww = warning;

if isstruct(ww)   % New Matlab R>12 Syntax
   try
      id = strcmp({ww.identifier},'all');
      if any(id)
         id = find(id);
         ww = ww(id(1)).state;
      else
         ww = '';
      end
   catch
      ww = '';
   end
elseif ~chkstr(ww)
   ww = '';
end

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%****************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  str = strhcat(str,del,n,nl)

% STRHCAT  Concatenates Strings into ONE
%
% STRHCAT( StringArray , Delimiter )
%   Forms one long String from the Strings in the
%   StringArray, delimited with the delimiter.
%   The EndDelimiter will be removed.
%
% STRHCAT( StringArray , Delimiter , N , NewLine )
%   Build a  NewLine after each N-th String.
%   default: N = 10;  NewLine = char(10);
%
% Example:  
%         >> strhcat({'apples' 'pies' 'prunes'},', ')
%    
%         ans =
%
%         apples, pies, prunes
%
%         >> strhcat({'apples';'pies';'prunes'},', ',2)
%    
%         ans =
%
%         apples, pies
%         prunes
%



Nin = nargin;

if Nin < 4
 nl = char(10);
end
if Nin < 3
 n = [];
end
if Nin < 2
 del = char((32*ones(1,3)));
end


if isempty(str)
 str = '';
 return
end


if ischar(str)
  str = cellstr(str);
end

str = str(:);

if isempty(n)
   n = size(str,1) + 1;
end

str( : , 2 ) = { del };

str(n:n:size(str,1),2) = { nl };

str(    size(str,1),2) = { '' };


str = permute( str , [ 2  1 ] );

str = cat(2,str{:});


