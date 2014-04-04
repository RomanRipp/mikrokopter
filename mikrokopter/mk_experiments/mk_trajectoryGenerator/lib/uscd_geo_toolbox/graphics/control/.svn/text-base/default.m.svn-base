function [s,v] = default(typ,varargin)

% DEFAULT   Set default properties of root
%
% [State,Values] = DEFAULT( Type , Property , Value , ... )
%
%  calls: set( 0 , 'Default<Type><Property>' , Value , ... )
%

Msg = '';
 nl = char(10);

Nin = nargin;

  s = cell(0,0);
  v = cell(0,0);

%*********************************************

if Nin < 1
   error('Input Type is missing');
end

if ~chkstr(typ,0)
   Msg = 'Type must be a String.';
else
   typ = cat( 2 , 'default' , typ );
end

%---------------------------------------------

if Nin == 1

   if ~isempty(Msg)
      error(Msg)
   end

   try
     s = get(0,typ);
     v = set(0,typ);
   catch
     error(lasterr);
   end

   return

end

%*********************************************

p = varargin;
n = prod(size(p));

%---------------------------------------------

if n == 1

   p = p{1};

   if ~chkstr(p,0)
      Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
               'Property must be a String.' ];
   end

   if ~isempty(Msg)
      error(Msg)
   end

   try
     s = get( 0 , cat(2,typ,p) );
     v = set( 0 , cat(2,typ,p) );
   catch
     error(lasterr);
   end

   return

end

%---------------------------------------------

ok = ( mod(n,2) == 0 );
if ok
   p  = reshape( p , 2 , n/2 );
   ok = chkcstr( p(1,:) , 1 );
end

if ~ok
   Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
            'Following Inputs must be Property-Value-Pairs.' nl ...
            'Properties must be Strings.' ];
end

if ~isempty(Msg)
   error(Msg)
end

%*********************************************

p(1,:) = permute( cellstr( cat( 2 , typ(ones(n/2,1),:) , char(p(1,:)) ) ) , [ 2  1 ] );

try
  set( 0 , p{:} );
  s = get( 0 , p(1,:) );
  v = set( 0 , typ );
catch
  error(lasterr);
end


%************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


