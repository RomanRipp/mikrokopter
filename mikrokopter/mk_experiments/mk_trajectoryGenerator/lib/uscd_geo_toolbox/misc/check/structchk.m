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

  elseif ~chkcstr(V(1:2:end-1));  % FieldNames

      Msg = sprintf('Multiple Input%s must be FieldName-Value-Pairs.',txt);

  else

     try 
        V = struct(V{:});
     catch
        Msg = sprintf('Error call STRUCT from Input%s.\n%s',txt,lasterr);
     end

  end

end

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
   n = size(str,1);
   if n == 1
      str = strrep(str,char(32),char(1));
   end
   str = cellstr(str);
   if n == 1
      str = strrep(str,char(1),char(32));
   end
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

