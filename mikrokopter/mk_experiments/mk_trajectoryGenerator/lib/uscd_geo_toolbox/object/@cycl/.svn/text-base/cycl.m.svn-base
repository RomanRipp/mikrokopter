function obj = cycl(v)

% CYCL  Constructor for CYCL-Object
%
% A CYCL-Object allows to use cyclic subscripts for Arrays.
% 
%   A = cycl(ARRAY)  defines an Array as CYL-Object.
%
% The Array can be a numeric, char, cell or struct array,
%  or any other user defined object with array-syntax.
%
% Subscripts, which exceeds the Array-Dimension,
%  will transformed cyclic into the dimension.
% 
% The error: ???  Index exceeds matrix dimensions.
%  in subscripted reference or assignment will never occur.
%
% The Syntax  A{:}  or  A.VALUE  returns the original Array,
%  to use it for calculations etc. with other functions,
%  and allows to redefine the Array of the CYCL-Object.
% 
% Extractig multiple values from Cell-Arrays ( {}-Syntax ) 
%  or Structures doesn't work!
%
% examples:
% 
%  A = cycl([3 7 9]);   % create CYL-Object A from Array
% 
%  A(1:7)               % cyclic subscripted reference
%
%% ans =
%
%%     3     7     9     3     7     9     3
%
%  A{:} + 1             % PLUS-Function with CYCL-Object
%
%% ans = 
%
%%     4     8    10    % numeric array
% 
%  sum(A{:})            % SUM-Function with CYCL-Object
%
%% ans =
%
%%    19
%
%  A{:} = A{:} + 1      % redefine the Array
%
%% A =
%
%%     4     8    10    % CYCL-Object
%
%  A{:} = A(1:5) + 1    % redefine the Array
%
%% A =
%
%%     5     9    11     5     9
%
%  A(1:7) = ( 1 : 7 )   % cyclic subscripted assignment
%
%% A =
%
%%     6     7     3     4     5
%
%  A(1:7) = [1 2 3 4 5 6 7 8]  
%
%% ??? Error using ==> cycl/subsasgn
%% Subscripted assignment dimension mismatch.
%
%


%***********************************************

c = 'cycl';   % Class
p = '';       % ParentClass


f = { 'VALUE' };  % FieldNames


%***********************************************

Nin = nargin;

%-----------------------------------------------
% Check Input

if ( Nin == 0 );
   v = [];
end


%***********************************************
% Check if Input is a CYCL-Object

if strcmp( class(v) , c )

   obj = v;

   return
 
end

%***********************************************
% Create new Object

   %-----------------------------------------------
   % Build Structure from f
 
   f      = f(:,[1 1])';
   f(2,:) = { {[]} }; 
   f      = struct(f{:});

   cl = clock;

   if ~isempty(p)

      obj = feval(p);

      obj = class(f,c,obj);

   else

      obj = class(f,c);

   end

   obj = setfield(obj,'VALUE',v);



