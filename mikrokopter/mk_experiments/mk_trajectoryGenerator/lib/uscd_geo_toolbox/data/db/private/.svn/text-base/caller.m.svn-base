function [c,sm,lm] = caller(dep,sep);

% CALLER   returns the History of the calling functions
%
%  C = CALLER;   
%
%  C = { FunctionName  LineNumer SubFunction };  3-Column-CellArray
%
%  CALLER, executed from the Workspace or a function which is
%   executed directly from the Workspace, returns an empty CellArray. 
%
%  [ C , ShortMsg , LongMsg ] = CALLER( Depth , Seperator )
%
%  returns a Messages as String:
%
%   ShortMsg = 'FCN1.SUB1:FCN2.SUB2:FCN3.SUB3: ... '
%    LongMsg = 'fcn1(lnr1).sub1:fcn2(lnr2).sub2: ... '
%
%  the String Separator will used instead of ':'
%
%  Depth gives the Number to recurse the DirectoryTree to Display, 
%   default: 0
%
%  see  also:  DBSTACK
%

Nin  = nargin;
Nout = nargout;

c  = cell(0,3);   % { Name LineNr SubFcn }
sm = '';
lm = '';

dp =  0 ;  % Default Depth
sp = ':';  % Default Seperator

s = dbstack;

if prod(size(s)) < 2
   return
end

%--------------------------------------------
% Check Depth

if Nin < 1
   dep = [];
end

ok = ( isnumeric(dep) & ~isempty(dep) &  ( prod(size(dep)) == 1 ) );
if ok
   ok = isfinite(dep);
end

if ~ok
   dep = dp;
end

%--------------------------------------------

s = s( 2 : end );  % without "caller.m" 

n = prod(size(s));

c = cell( n , 3 );  % without FunctionName

c(:,3) = {''};

fs = filesep;

for ii = 1 : n

    cc = sepname(s(ii).name,dep,fs,0);

    c{ii,2} = s(ii).line;

    j1 = find( cc == '(' );
    j2 = find( cc == ')' );
    if ~( isempty(j1) | isempty(j2) ) 
        j2 = j2(end);
        j1 = j1( find( j1 < j2 ) );
        if ~isempty(j1)
            c{ii,3} = cc( j1+1 : j2-1 );
            cc      = cc( 1 : j1-1 );
        end
    end

    j1 = find( cc == '.' );
    j2 = find( cc == ' ' );

    if ~isempty(j1)
        j1 = max(j1);
    end
    if ~isempty(j2)
        j2 = min(j2);
    end

    if ~( isempty(j1) & isempty(j2) )
        cc = cc( 1 : min([j1,j2])-1 );
    end

    c{ii,1} = cc;

end

if Nout < 2
   c = c(2:n,:);  % without FunctionName
   return
end

%--------------------------------------------
% Check Separator

if Nin < 2
   sep = '';
end

if ~( ischar(sep) & ~isempty(sep) &  ( prod(size(sep)) == size(sep,2) ) )
   sep = sp;
end

ns = size(sep,2);

%--------------------------------------------
% ShortMessage

sm = upper(c(:,1));

is_sub = find(~strcmp(c(:,3),''));
if ~isempty(is_sub)
    is_sub = is_sub(:)';
    for ii = is_sub
        sm{ii} = sprintf('%s.%s',sm{ii},upper(c{ii,3}));
    end
end

sm = sprintf( ['%s' sep] , sm{n:-1:1} );

sm = sm( 1 : size(sm,2)-ns );

if Nout < 3
   c = c(2:n,:);  % without FunctionName
   return
end

%--------------------------------------------
% LongMessage

lm = permute(c(n:-1:1,:),[2 1]);

lm = sprintf( ['%s[%.0f].%s' sep] , lm{:} );

lm = lm( 1 : end-ns );
lm = lm( 1 : end-(lm(end)=='.') );

jj = findstr(lm,['.' sep]);
if ~isempty(jj)
    lm(jj) = [];
end

c = c(2:n,:);  % without FunctionName


%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = sepname(name,dep,sep,mode);

% SEPNAME  separates Name
%
% String = SEPNAME( Name , Depth , Seperator ,Mode )
%
% Depth gives the number of Recursions
% 
% Use Depth = NaN  to get all parts of Name in String.
%  In this case String is a CellStringArray.
%
% Mode =  1  Recursion from Start
%        -1  Recursion from End
%         0  Recursion from End, 
%             special Handling for Class- and Private-Directories
%
% Defaults:  Depth     = 0
%            Seperator = FILESEP
%            Mode      = 0    ( 1  if Depth = NaN )
%


Nin = nargin;

if Nin < 1
   name = '';
end

if Nin < 2
   dep = 0;
end

dep_nan = isnan(dep);

if Nin < 3
   sep = filesep;
end

if Nin < 4
   mode = dep_nan;
end

%********************************************

if dep_nan
   str = cell(1,0);
else
   str = char(zeros(1,0));
end

if isempty(name)
   return
end

if ~chkstr(name)
   error('Name must be a String.')
end

if ~( chkstr(sep,1) & ( prod(size(sep)) == 1 ) )
   error('Seperator must be a single Character.')
end

n = size(name,2);

%---------------------------------------------
% Find Seperator in Name

is = ( double(name) == double(sep) );

if all(is)
   if dep_nan
      str    = cell(1,n-1);
      str(:) = { '' };
   end      
   return
end

%---------------------------------------------

i0 = ~is(1);
i1 = ~is(n);

is = cat( 2 , ones(1,i0) , is , ones(1,i1) );

is = find( is );

is = is(:);

ni = size(is,1) - 1;

if ni == 0 
   return
end
     
%---------------------------------------------
% [ Start  End ]

ind = ( 1 : ni ) ;

is  = cat( 2 , is(ind)+1 , is(ind+1)-1 ) - i0;

%---------------------------------------------
% Take care for duplicate Seperators

if ~dep_nan | ( double(sep) == char(32) )

   is( find( is(:,1) > is(:,2) ) , : ) = [];

end

%---------------------------------------------

ni = size(is,1);

if ni == 0
   return
end

%---------------------------------------------
if dep_nan
   
   ind = [ 1  ni ];

  flip = ~( mode == 1 );

   ind = ind( [ 1  2 ] + [ 1 -1 ]*flip );

   ind = ( ind(1) : 1-2*flip : ind(2) );

   is = is(ind,:);

   str = cell(1,ni);
  
   for ii = 1 : ni

       str{ii} = name( is(ii,1) : is(ii,2) );

   end

   return

end

%---------------------------------------------

ii = ni - 1 * ( ni > 1 );

nn = name( is(ii,1) : is(ii,2) );

ic = strcmp( nn(1) , '@'       );
ip = strcmp( nn    , 'private' );

id = 1 * ic + 2 * ip;

dep = dep + 1 + id * ( mode == 0 ) * ( ni > 1 );

dep = min( dep , ni );

%---------------------------------------------

is(1,1) = is(1,1) - 1 * ( is(1,1) > 1 );  % Start incl. Seperator

ind = ( 1 : ni-1 );

is(ind,2) = is(ind,2) + 1;                % End incl. Seperator

is(:,2) = is(:,2) - is(:,1) + 1;          % Length

%---------------------------------------------

ind = ( 1 : dep ) + ( ni - dep ) * (~( mode == 1 ));

is  = is(ind,:);

ind = grp2ind(is(:,1),is(:,2));

str = name(ind);

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ii = grp2ind(i0,l);

% GRP2IND  Built IndexVector from StartIndex and Length
%
% Index = GRP2IND( StartIndex , GroupLength )
%

if isempty(i0);
   ii = [];
   return
end

si = size(i0);

if ( sum( si > 1 ) > 1 )
   error('StartIndex must be a Vector.');
end

i0 = i0(:);
l  =  l(:);

if ~isequal(size(i0,1),size(l,1))
   error('Size of StartIndex and GroupLenght must be the same.');
end

n = size(l,1);

ii = ones(sum(l),1);
jj = cumsum( cat(1,1,l) , 1 );

ii(jj(1:n)) = i0;

if n > 1
   ii(jj(2:n)) = ii(jj(2:n))-(i0(1:n-1,1)+l(1:n-1)-1);
end

ii = cumsum(ii,1);

if n == 1
   return
end

jj = find( si > 1 );

if jj == 1
   return
end
  
perm = cat( 2 , (1:jj-1)+1 , 1 , (jj+1:size(si,2)) );

ii = permute(ii,perm);

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

