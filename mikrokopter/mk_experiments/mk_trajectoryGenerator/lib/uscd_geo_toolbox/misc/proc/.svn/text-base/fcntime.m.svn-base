function [mdt,dt,str] = fcntime(varargin);

% FCNTIME    Returns Time to evaluate a Function
%
% EvaluationTime = FCNTIME( FunctionName , varargin )
%
% [ EvaluationTime , TotalTime ] = FCNTIME( LoopNumber , FunctionName , varargin )
%
%   returns the TotalTime used for [LoopNumber] Loops.
%
% A imaginary part of LoopNumber is the Number of OutputArguments
%  for the FunctionCall: FCNTIME( LoopNumber + OutputNumber * i , ... )
% 
% [ EvaluationTime , TotalTime , String ] = FCNTIME( ... )
%
%   returns a InfoCellString with EvaluationTime and TotalTime.
%
%  see also:  LOOPTIME, FEVAL
%


Nin  = nargin;
Nout = nargout;

if Nin < 1
   error('Not enough Input Arguments.');
end

n = 1;

%-----------------------------------------------

v  = varargin;
nv = prod(size(v));

if ( isnumeric(v{1}) & ( prod(size(v{1})) == 1 ) );
   n  = v{1};
   if nv == 1
      error('Input Function is missing.');
   end
   v  = v(2:nv);
end

%-----------------------------------------------

msg = '';
nl  = char(10);

m = imag(n);
n = real(n);

if ~( ( mod(n,1) == 0 ) & ( n >= 0 ) );
   msg = 'LoopNumber must be a positive Integer.';
end

if ~( ( mod(m,1) == 0 ) & ( m >= 0 ) );
   msg = [ msg  nl(1:(end*(~isempty(msg)))) ...
           'OutputNumber must be a positive Integer.' ];
end

n = max(1,n);

if ~chkstr(v{1},1);
   msg = [ msg  nl(1:(end*(~isempty(msg)))) ...
           'Function must be a String'          ];
end

if ~isempty(msg)
   error(msg);
end

%-----------------------------------------------

t0 = looptime;

if m == 0

   for ii = 1 : n
       feval(v{:});
   end

else

   m = cell(1,m);

   for ii = 1 : n
       [m{:}] = feval(v{:});
   end

end

   
[t1,dt,mdt] = looptime(n,n,t0,[]);

if ( Nout > 0 ) & ( Nout < 3 )
   return
end

%-----------------------------------------------
% InfoString

str = cell(2,1);

app = char( 's' * ones(1,(n>1)) );

str{1} = sprintf('EvaluationTime: %8.3f sec',mdt);
str{2} = sprintf('Time for %.0f Loop%s: %8.3f sec',n,app,dt);

ns = size(char(str),2);

for ii = 1 : size(str,1)
    str{ii} = cat( 2 , char(32*ones(1,ns-size(str{ii},2))) , str{ii} );
end

if ( Nout > 0 )
   return
end

%-----------------------------------------------
% Display

bl = char(32*ones(1,3));

fprintf(1,nl);
fprintf(1,bl);
fprintf(1,str{1});
fprintf(1,nl);

if n > 1
   fprintf(1,bl);
   fprintf(1,str{2});
   fprintf(1,nl);
end

fprintf(1,nl);

clear mdt

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

