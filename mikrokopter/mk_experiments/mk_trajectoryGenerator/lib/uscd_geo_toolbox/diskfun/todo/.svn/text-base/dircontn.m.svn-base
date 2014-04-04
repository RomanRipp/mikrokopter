function [c,l,t] = dircont(pfad,varargin)

% DIRCONT  Returns Contents of Directory
%
% [ C , LISTE , TEXT ] = DIRCONT( Pfad , Mode , WildCard )
%
% gives Contents of Directory "Pfad"
%
%  C    = { IsDir  Date  Bytes  Name  DateNum };
%
%  LISTE = CellString(C)
%
%  TEXT  = STRHCAT( LISTE,'',1);
%
%  Mode == {'n'} , sort by Name            of c(:,4), default
%          {'a'}           Alphabetical
%          {'c'}           Character
%          {'l'}           Letter
%
%  Mode == {'t'} , sort by Time == DateNum of c{:,5}
%          {'d'}           DateTime
%          {'z'}           Zeit
%
%  Mode == {'s'} , sort by Size            of c(:,3)
%          {'b'}           Bytes
%
% If WildCard ends with FileSep, only Directory will returned.
%


Nin  = nargin;
Nout = nargout;

c = cell(0,5);  % { IsDir Date  Bytes  Name DateNum}
l = cell(0,1);
t = '';

%------------------------------------------------
% Get Inputs

[msg,pfad,mode,wc,recurse,follink] = checkin(pfad,varargin{:});

if ~isempty(msg)
   error(msg);
end

%------------------------------------------------
% Get DirectoryContents

fs = filesep;

is_win = strcmp(upper(computer),'PCWIN');

d = get_cont(pfad,fs,mode,recurse,follink,is_win,wc,0);


N = size(d,1);


if N == 0
  return
end

%------------------------------------------------

% Format for Bytes
form = '%8.0f';  

% Seperator: [ sep0 Date sep1 Bytes sep2 Name ];

sep0 = ' ';
sep1 = '  ';
sep2 = '  ';


l = cell(N,1);
c = cell(N,5);

for ii = 1 : N
 
    c{ii,1} = d(ii).isdir * ( 1 - 2 * strcmp(d(ii).name,'..') );
    c{ii,2} = d(ii).date;
    c{ii,3} = d(ii).bytes;
    c{ii,4} = d(ii).name;
 
    if isempty(d(ii).date)
       c{ii,5} = NaN;
    else
       c{ii,5} = datenum(ger2eng(d(ii).date));  
    end

end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if is_win

   c(:,4) = lower( c(:,4) );

end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


nd = size(char(c{:,2}),2);

for ii = 1 : N

     fs_end = size(fs,2) * ( c{ii,1}                      &  ...
                            ~strcmp(c{ii,4}(end) ,  fs  ) &  ...
                            ~strcmp(c{ii,4}      , '..' )        );

    bytes  = sprintf( form , c{ii,3});
    if d(ii).bytes == 0
       bytes = char( 32 * ones( 1 , size(bytes,2) ) );
    end

    if isempty(c{ii,2})
       date = char( 32 * ones(1,nd) );
    else
       date = cat( 2 , c{ii,2} , char( 32 * ones(1,nd-size(c{ii,2},2)) ) );
    end

    l{ii}  = cat( 2 , sep0 , date                ,  ...
                      sep1 , bytes               , ...
                      sep2 , c{ii,4} , fs(1:fs_end)        );

 
end


c0 = cell(0,5);
l0 = cell(0,1);

bad = strcmp( c(:,4) , '.' );

if get_dir

   bad = ( bad  |  ~cat(1,c{:,1}) );

else

   jj = find( strcmp( c(:,4) , '..' ) );

   if ~isempty(jj);

     c0 = c(jj,:);
     l0 = l(jj);

     c(jj,:) = [];
     l(jj)   = [];

   end

end



bad = find(bad);

c(bad,:) = [];
l(bad,:) = [];

% Sort alphabetical

is_dir = cat(1,c{:,1});

ii = find( is_dir);
jj = find(~is_dir);

switch mode
  case { 't' 'd' 'z' }

    % via DateNum

    [h,s_ii] = sort( (-1)*cat(1,c{ii,5}) );
    [h,s_jj] = sort( (-1)*cat(1,c{jj,5}) );

  case { 'n' 'a' 'c' 'l' }

    % via Name

    [h,s_ii] = sort( c(ii,4) );
    [h,s_jj] = sort( c(jj,4) );

  case { 's' 'b' }

    % via Size

    [h,s_ii] = sort( c(ii,4) );  % Directories via Name
    [h,s_jj] = sort( (-1)*cat(1,c{jj,3}) );

  otherwise

       s_ii  = ( 1 : prod(size(ii)) );
       s_jj  = ( 1 : prod(size(jj)) );

end

c = cat( 1 , c0 , c(ii(s_ii),:) , c(jj(s_jj),:) );
l = cat( 1 , l0 , l(ii(s_ii),:) , l(jj(s_jj),:) );


if Nout == 3
   t = strhcat(l,'',1);
end

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  d = get_cont(pfad,fs,mode,recurse,follink,is_win,wc,z);


if isempty(wc)

   d = dir(pfad);

   return

end

fs = fs(1:(end*(~strcmp(pfad(end),fs))));

d = struct( 'name' , {} , ...
            'date' , {} , ...
           'bytes' , {} , ...
           'isdir' , {}        );

d = d(ones(0,1));


for ii = 1 : size(wc,1)

    d1 = dir( cat(2,pfad,fs(1:(end*(~isempty(wc{ii,1})))),wc{ii,1}) );

    ok = ~wc{ii,2} + cat( 1 , d1.isdir ) + recurse;

    d = cat( 1 , d , d1(find(ok)) );

end

n = { d.name };

if ~( z == 0 )
   bad = find( strcmp( n , '..' ) );
   n(bad) = [];
   d(bad) = [];
end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if is_win

   n = lower(n);

end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

[n,si] = sort(n);

n = find( sum( diff(double(n),1,1) , 2 ) == size(n,2) ) + 1;

d(n) = [];

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function dir2cell(d,pre,mode,is_win,fs,dok)


%------------------------------------------------

% Format for Bytes
form = '%8.0f';  

% Separator: [ sep0 Date sep1 Bytes sep2 Name ];

sep0 = ' ';
sep1 = '  ';
sep2 = '  ';


l = cell(N,1);
c = cell(N,5);

for ii = 1 : N
 
    c{ii,1} = d(ii).isdir * ( 1 - 2 * strcmp(d(ii).name,'..') );
    c{ii,2} = d(ii).date;
    c{ii,3} = d(ii).bytes;
    c{ii,4} = cat(2,pre,d(ii).name);
 
    if isempty(d(ii).date)
       c{ii,5} = NaN;
    else
       c{ii,5} = datenum(ger2eng(d(ii).date));  
    end

end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if is_win

   c(:,4) = lower( c(:,4) );

end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


nd = size(char(c{:,2}),2);

for ii = 1 : N

     fs_end = size(fs,2) * ( c{ii,1}                      &  ...
                            ~strcmp(c{ii,4}(end) ,  fs  ) &  ...
                            ~strcmp(c{ii,4}      , '..' )        );

    bytes  = sprintf( form , c{ii,3});
    if d(ii).bytes == 0
       bytes = char( 32 * ones( 1 , size(bytes,2) ) );
    end

    if isempty(c{ii,2})
       date = char( 32 * ones(1,nd) );
    else
       date = cat( 2 , c{ii,2} , char( 32 * ones(1,nd-size(c{ii,2},2)) ) );
    end

    l{ii}  = cat( 2 , sep0 , date                ,  ...
                      sep1 , bytes               , ...
                      sep2 , c{ii,4} , fs(1:fs_end)        );

 
end


c0 = cell(0,5);
l0 = cell(0,1);

bad = strcmp( c(:,4) , '.' );

if ~dok

   bad = ( bad  |  ~cat(1,c{:,1}) );

else

   jj = find( strcmp( c(:,4) , '..' ) );

   if ~isempty(jj);

     c0 = c(jj,:);
     l0 = l(jj);

     c(jj,:) = [];
     l(jj)   = [];

   end

end



bad = find(bad);

c(bad,:) = [];
l(bad,:) = [];

% Sort alphabetical

is_dir = cat(1,c{:,1});

ii = find( is_dir);
jj = find(~is_dir);

switch mode
  case { 't' 'd' 'z' }

    % via DateNum

    [h,s_ii] = sort( (-1)*cat(1,c{ii,5}) );
    [h,s_jj] = sort( (-1)*cat(1,c{jj,5}) );

  case { 'n' 'a' 'c' 'l' }

    % via Name

    [h,s_ii] = sort( c(ii,4) );
    [h,s_jj] = sort( c(jj,4) );

  case { 's' 'b' }

    % via Size

    [h,s_ii] = sort( c(ii,4) );  % Directories via Name
    [h,s_jj] = sort( (-1)*cat(1,c{jj,3}) );

  otherwise

       s_ii  = ( 1 : prod(size(ii)) );
       s_jj  = ( 1 : prod(size(jj)) );

end

c = cat( 1 , c0 , c(ii(s_ii),:) , c(jj(s_jj),:) );
l = cat( 1 , l0 , l(ii(s_ii),:) , l(jj(s_jj),:) );

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,pfad,mode,wc,recurse,follink] = checkin(pfad,varargin);


% CHECKIN  Check of Inputs 


msg = '';
nl  = char(10);

fs  = filesep;

Nin = nargin;


%*********************************************************
% Defaults

recurse = 0;
follink = 0;

mode = 'a';

wc = cell(0,2);

%*********************************************************
% Get and Check Inputs

%---------------------------------------------------------
% Check Pfad and String

if Nin < 1
   pfad = cd;
elseif isempty(pfad)
   pfad = fs;
end

if ~( ischar(pfad)  &  ( size(pfad,2) == prod(size(pfad)) ) )
  msg = 'Directory must be a String.';
end

if ~isempty(pfad) 
   pfad = pfad( 1 : ( end - any(strcmp(pfad(end),{ '.'  '*' })) ) );
end

if isempty(pfad)
   pfad = cd;
end

%---------------------------------------------------------

if Nin < 2
   return
end


%---------------------------------------------------------
% Get Options from VarArg

for ii = 1 : Nin-1

   val = varargin{ii};

   nv = prod(size(val));

   recurse = ( recurse | isequal(val,'-r') );
   follink = ( follink | isequal(val,'-l') );
      
   if iscellstr(val)

      if ~isempty(val{1})
        mode = lower(val{1}(1));
      end

   elseif ( ischar(val)  &  ( nv > 0 )  & ...
            ( nv == size(val,2) )     )

      switch val

        case '-r'

           recurse = 1;

        case '-l'

           follink = 1;

        otherwise
 
           dd = strcmp(wc(nv),fs);

           wc = wc = cat( 1 , wc , { val( 1 : nv-dd )  dd };
      end

   end

end


%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [y,ndx] = sortchar(x,col)

% SORTCHAR Sort rows in ascending order.
%
%   SORTCHAR(X) sorts the rows of the matrix X in ascending order as a
%   group.  For ASCII strings, this is the familiar dictionary sort.
%   When X is complex, the elements are sorted by ABS(X). Complex
%   matches are further sorted by ANGLE(X).
%
%   SORTCHAR(X,COL) sorts the matrix based on the columns specified in
%   the vector COL.  For example, SORTCHAR(X,[2 3]) sorts the rows of X
%   by the second and third columns of X.
%
%   [Y,I] = SORTCHAR(X) also returns an index matrix I such that Y = X(I,:).
%
%   See also SORT.

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 1.11 $  $Date: 1998/04/02 18:05:26 $

if nargin<1, error('Not enough input arguments.'); end

if ndims(x)>2, error('X must be a 2-D matrix.'); end

[m,n] = size(x);

if nargin<2, col = 1:n; end

% Sort back to front
if m>0
  ndx = (1:m)';
else
  ndx = [];
end

x = double(x);

for i=length(col):-1:1,
  [v,ind] = sort(x(ndx,col(i)));
  ndx = ndx(ind);
end

y = char( x(ndx,:) );

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = ger2eng(str)

% GER2ENG   Converts German MonthNotation to English
%
% EnglishString = GER2ENG( GermanString );
%
%  mär --> mar
%  mai --> may
%  okt --> oct
%  dez --> dec
%
% In other case you get problems to use DATEVEC, DATENUM
%

if nargin < 0
   str = '';
   return
end

if ~( iscellstr(str) | ischar(str) );
   error('Input must be a CharArray or CellStringArray.');
end


is_char = ischar(str);

if is_char
   str = cellstr(str);
end

rep = { ...
  'mär'  'mar'
  'mai'  'may'
  'okt'  'oct'
  'dez'  'dec'  };


for ii = 1 : prod(size(str))

    s = str{ii};

    if ~isempty(s)  

        s = lower(s); 

        for jj = 1 : size(rep,1)

            kk = findstr( s , rep{jj,1} );

            if ~isempty(kk)

               n = size(rep{jj,1},2);

               for ll = kk(:)'

                   ind = ll + (1:n) - 1;     

                    nn = find( ~( double( s(ind) ) == double( str{ii}(ind) ) ) );

                    str{ii}(ind) = rep{jj,2};

                    str{ii}(ind(nn)) = upper(str{ii}(ind(nn)));

               end
               % ll == kk(findstr)
            end
            % ~isempty(findstr) 
        end
        % jj = rep
    end
    % ~isempty(str{ii})
end
%ii

if is_char
   str = char(str);
end
