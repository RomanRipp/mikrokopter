function [f,p,n] = whichfile(varargin)

% WHICHFILE   Find Files in Directory or Matlab's SearchPath by Suffix
%
% [ FileName , DirName , FullName ] = WHICHFILE( SUF1 , SUF2 , ... )
%
% List Files in Matlab's SearchPath with the Syntax: *SUF#
%
% [ FileName , DirName , FullName ] = WHICHFILE( {PathList} , SUF1 , SUF2 , ... )
%
% List Files in the Directories of PathList with the Syntax: *SUF#
%
% WHICHFILE( {PathList} , ... , '-r' ) works recursive
%
% WHICHFILE( {PathList} , ... , '-R' ) works recursive, follow Links
%
% The Inputs for FileSuffixes SUF# must be Strings.
%
% FileName, DirName and Fullname are CellArray's of Strings.
%

Nout = nargout;

f = cell(0,1);
p = cell(0,1);
n = cell(0,1);

if nargin == 0
   return
end

pfad = {};
if chkcstr(varargin{1},1);
   if nargin == 1
      return
   end
   pfad = varargin{1};
   varargin = varargin(2:end);
end

[ok,v] = chkcstr(varargin);

if ~ok
   error('Inputs must be Strings');
end

rec = NaN;

if  isempty(pfad)
    pfad = cat(2,sepname(matlabpath,NaN,pathsep),{'.'});
else
    ii  = ( strcmp(v,'-r') | strcmp(v,'-R') );
    if any(ii)
       if all(ii)
          return
       end
       rec = v(find(ii));
       rec = rec{end}(2);
       rec = ( rec == upper(rec) ); % True for follow Links
       rec = ( rec | ~isunix );     %  !!!
       v   = v(find(~ii));
    end
    pfad = pfad(:)';
end

nv = prod(size(v));

fs = filesep;

ii = strcmp(pfad,'');    % Empty PathName, use current
if any(ii)
   ni =  sum(ii);
   ii = find(ii);
   if any( strcmp(pfad,'.') | strcmp(pfad,['.' fs]) )
      pfad(ii) = [];
   elseif ni > 1
      pfad(ii(1)) = {'.'};
      pfad(ii(2:ni)) = [];
   else
      pfad(ii) = {'.'};
   end
end

for pp = pfad

    if ~( pp{1}(end) == fs );
          pp{1} = cat( 2 , pp{1} , fs );
    end

    for ii = 1 : nv

        d = dir(cat(2,pp{1},'*',v{ii}));

        if ~isempty(d)
            ok = ~cat(1,d.isdir); 
            if any(ok)
               nk = sum(ok);
               ok = find(ok);
                f = cat( 1 , f , {d(ok).name}' );
                p = cat( 1 , p , pp(ones(nk,1)) );
            end
        end
 
    end

    if ~isnan(rec)

        d = dir(pp{1});

        if ~isempty(d)

            isd = cat(2,d.isdir);
            dnm = { d.name };

            isd = ( isd & ~( strcmp(dnm,'.') | strcmp(dnm,'..') ) );

            if any(isd)

               isd = find(isd);

               dnm = dnm(isd);

               for d = dnm(:)'

                   pf = cat( 2 , pp{1} , d{1} );

                   ok = rec;

                   if ~ok
                       l  =  ls('-ld',['"' pf '"']);  % Without FileSep at End !!!
                       ok = ~isequal(l(1),'l');
                   end

                   if ok
                      [f1,p1] = whichfile({pf},varargin{:});
                       f = cat( 1 , f , f1 );
                       p = cat( 1 , p , p1 );
                   end

               end

           end

        end

    end

end   

if Nout == 0

   s = strhcat(permute(cat(2,p,f),[2 1]),'',2);

   fprintf(1,'\n%s\n\n',s);

   clear f

elseif Nout == 3

   n = p;
   for ii = 1 : size(n,1)
       n{ii} = [ n{ii} f{ii} ];
   end

end

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

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

