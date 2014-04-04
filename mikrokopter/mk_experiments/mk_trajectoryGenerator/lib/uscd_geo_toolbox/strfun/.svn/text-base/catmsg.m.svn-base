function [msg,nl] = catmsg(msg,str,n);

% CATMSG   Expands Message
%
% [ NewMsg , Separator ] = CATMSG( OldMsg , String , BlankNumber )
%
%  OldMsg       Original MessageString
%  String       CharacterArray of CellStringArray of Strings to append to OldMsg, 
%                 separated with NewLine and Blanks
%  BlankNumber  Number of Blanks, added to NewLine to separate Strings.
%                 Separator = cat( 2 , char(10) , char(32*ones(1,BlankNumber)) )
%
%  If the Value for BlankNumber has an Imaginary Part, the UpperCase Name of the 
%   calling M-File will added to OldMsg: 
%    OldMsg = cat( 2 , MFILENAME , ': ' , OldMsg )
%   The BlankNumber will expand by the Length of [ MFILENAME ': ' ]
%
% example:  
%
%  txt = { 'Size of X and Y must be agree.'
%          'Invalid Method.'                }
%
%  msg = catmsg( 'Invalid Inputs', txt )
%  msg = catmsg( 'Invalid Inputs', txt , 8 )
%  msg = catmsg( 'Invalid Inputs', txt , 1+i )
%

Nin = nargin;

err = '';
 nl = char(10);

%*********************************************************************
% Check Inputs

%---------------------------------------------------------------
% Message

if Nin < 1
   msg = char(zeros(1,0));
   return
end

if isempty(msg)
   msg = char(zeros(1,0));
else
   if ~chkstr(msg,0)
      err = [ err  nl(1:(end*(~isempty(err)))) ...
              'Message must be a String.' ];
   end
end

%---------------------------------------------------------------
% String

if Nin < 2
   str = cell(0,0);
end

if ~( iscell(str) & isempty(str) )
   [ok,str] = chkcstr(str,0);
   if ~ok
      err = [ err  nl(1:(end*(~isempty(err)))) ...
             'String must be a CharacterArray or CellStringArray.' ];
   end
end

%---------------------------------------------------------------
% BlankNumber

if Nin < 3
   n = 0;
elseif ~( isnumeric(n) & ( prod(size(n)) == 1 ) )
      err = [ err  nl(1:(end*(~isempty(err)))) ...
              'BlankNumber must be a single Number.' ];
end
 
%---------------------------------------------------------------

if ~isempty(err)
   error(err)
end

%*********************************************************************

%---------------------------------------------------------------
% Prepend CALLER

if ( imag(n) == 0 )

   c = char(zeros(1,0));

else

   c = dbstack;

   if prod(size(c)) < 2
      c = 'Matlab';
   else
      c = sepname(c(2).name,0,filesep,0);  % Calling M-FileName
      c = upper(c);
     jj = find( double(c) == double('.') );
     if ~isempty(jj)
        c = c(1:max(jj)-1); 
     end
   end

   c = cat( 2 , c , ': ' );

   n = real(n) + 0 * size(c,2); % !!!!

end

%---------------------------------------------------------------

nl = cat( 2 , nl , char( 32 * ones(1,ceil(n))) );

if ~isempty(str)
   str = strrep(str,char(10),nl);
   str = strhcat(str,nl);
end

if isempty(str) & isempty(msg)
   return
end

inl = ~( isempty(msg) | isempty(str) );

msg = cat( 2 , c , msg ,  nl(1:(end*inl)) , str );

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

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


