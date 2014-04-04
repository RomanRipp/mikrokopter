function assign(n,v)

% ASSIGN  Assign variable in current workspace
%
% ASSIGN('name',V) assigns the variable 'name' in the
%    current workspace the value V
%
% see also: ASSIGNIN, EVAL, EVALIN
%

Nin = nargin;
if Nin == 0
   return
end

if ~( strcmp( class(n) , 'char' )     & ...
       ( prod(size(n)) == size(n,2) ) & ~isempty(n) )
    error('Name of Variable must be a String.');
end

if ~all( ( ( '0' <= n ) & ( n <= '9' ) ) | ...
         ( ( 'A' <= n ) & ( n <= 'Z' ) ) | ...
         ( ( 'a' <= n ) & ( n <= 'z' ) ) | ...
           ( '_' == n )    );
    error('Invalid Characters in Name of Variable.');
end


if Nin == 1
   v = [];
end

assignin('caller',n,v);



