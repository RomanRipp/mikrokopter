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


fs = filesep;


c = cell(0,5);  % { IsDir Date  Bytes  Name DateNum}
l = cell(0,1);
t = '';


if nargin < 1
  pfad = cd;
else
  if ~( ischar(pfad)  &  ( prod(size(pfad)) == size(pfad,2) )  )
    error('Input must be a String.');
  end
  if isempty(pfad)
     pfad = fs;
  end
end
 
mode = 'a';
wc   = '';


for ii = 1 : Nin-1

   val = varargin{ii};

   if iscellstr(val)

      if ~isempty(val{1})
        mode = lower(val{1}(1));
      end

   elseif ( ischar(val)  &  ~isempty(val) & ...
            ( prod(size(val)) == size(val,2) )     )

      wc = val;

   end

end

get_dir = 0;
if ~isempty(wc)
   get_dir = strcmp(wc(end),fs);
   if get_dir
      wc = wc(1:end-1);
   end
end



if ~isempty(wc)
  pfad = cat( 2 , pfad , fs(1:(end*(~strcmp(pfad(end),fs)))) , wc );
end

% Format for Bytes
form = '%8.0f';  

% Seperator: [ sep0 Date sep1 Bytes sep2 Name ];

sep0 = ' ';
sep1 = '  ';
sep2 = '  ';

%------------------------------------------------
% Start

d = dir(pfad);

N = size(d,1);


if N == 0
  return
end


sep0 = ' ';
sep1 = '   ';
sep2 = '   ';

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
if strcmp(upper(computer),'PCWIN')

   c(:,4) = lower( c(:,4) );

end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


nd = size(char(c{:,2}),2);

for ii = 1 : N*( Nout >= 2 )

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


if get_dir

   bad = ( strcmp( c(:,4) , '.' )  |  ~cat(1,c{:,1}) );

else

   jj = find( strcmp( c(:,4) , '..' ) );

   if ~isempty(jj);

     c0 = c(jj,:);
     l0 = l(jj);

     c(jj,:) = [];
     l(jj)   = [];

   end

   bad = strcmp( c(:,4) , '.' );

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
