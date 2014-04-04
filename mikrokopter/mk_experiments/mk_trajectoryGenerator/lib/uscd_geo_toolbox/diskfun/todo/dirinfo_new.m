function [inf,str] = dirinfo(varargin)

% DIRINFO  Returns Info about Contents of Directory
%
% INF = DIRINFO( Pfad , [options] )
%
% gives DirectoryInfo CellArray INF:
%
%  INF    = { ID  Date  Bytes  Name  ...
%             SubINF  SubSize  SubNumbers       }   %  7 Columns
%
%  with: ID  =  1    for Directory
%               0    for Files
%              -1    for Links
%
%        SubNumbers = [ NumberOfDirectories  NumberOfFiles ]
%  
%        SubINF     =  CellArray INF  for  Directory
%               
% 
% Options:  -r  recurse trough SubDirectories
%           -l  give CellstringArray of List in 8 Column of D
%
%           -Ffilename  writes LIST into File with name filename
%                        using DIRLIST
%
%             Same Result like Option '-Ffilename':
%                 >> info = dirinfo('.','-r','-l');
%                 >> dirlist(info{8},filename);
%
% The Option "-l" returns a CellStringArray of DirectoryContents
%   in the 8. Element of INF, with the String-Element-Columns:
%   { ID  Date  Byte  Name  [NDir NFiles] }
%
%   with: ID  = 'D'   for Directory
%               'F'   for Files
%               'L'   for Links
%
% [ INF , LISTE ] = DIRINFO( ... )
%
%  Returns a String, build from CellArray, described above.
%


VarArg = varargin;
VarArg = VarArg(:);

Nin  = size(VarArg,1);
Nout = nargout;


nl = char(10);

%---------------------------------------
% Search for Options

listout  = ( Nout == 2 );

recurse  = 0;
listmode = 0;

file     = '';

ismode = zeros(Nin,1);

for ii = 1 : Nin

  mode = VarArg{ii};

  if ischar('mode')    &  ...
     size(mode,1) == 1 & ...
     size(mode,2) >= 2    

    if mode(1) == '-'

      if ( ( mode(2) == 'F'   ) & ...
           ( size(mode,2) > 2 )       );
        file = mode(3:size(mode,2));
        file(find(abs(file)==32)) = [];
        ismode(ii) = ~isempty(file);
      end

        recurse = ( recurse  | ( mode(2) == 'r' )  );
       listmode = ( listmode | ( mode(2) == 'l' )  );     

       ismode(ii) = ( ismode(ii)         | ...
                      ( mode(2) == 'r' ) | ...
                      ( mode(2) == 'l' )       );

       if ~ismode(ii)
         fprintf([nl 'DIRINFO: Invalid Option: '  mode nl ]);
       end

    end 
    % '-'
  end
  % ischar
end

filemode = ~isempty(file);

list = ( listmode  |  listout );



VarArg(find(ismode)) = [];
VarArg               = VarArg(:);
Nin                    = size(VarArg,1);

%---------------------------------------
% Search for Pfad

pfad = '';
if Nin >= 1
  if ischar(VarArg{1})
    pfad = VarArg{1};
  else
    fprintf([nl 'DIRINFO: Pfad must be a String.'  nl ]);
  end
end


if isempty(pfad)
 pfad = cd;
end

if any(strcmp(pfad,{ '.'  '*' }))
  pfad = cd;
end

fs = filesep;

if ~strcmp(pfad(size(pfad,2)),fs)
    pfad = [ pfad fs ];
end

np = size(pfad,2);


%------------------------------------------------
% Prepare Output's
 
% M Columns in CellArray:
%  { IsDir  Date  Byte  Name  ...
%    SubCell  SubBytes  [ NDir NFiles ]  }

  M = 7;
  
%  { IsDir  Date  Byte  Name  SubCell  SubBytes NSub List }


inf   = cell(0,M+listmode);
liste = '';

% Format for Bytes
bform = '%10.0f';  

% Format for [ NDir NFiles ]  
mform = { '#%.0f(%0.f)'  '#%.0f' };  % { InclLink    NoLink }
nform = sprintf('%%s   %%s%s',nl);   %   'DirForm  FileForm'

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ini,cnf] = get_info([pfad name{ii} fs],dpt,verb,rec,sm,fs);

ni = 9;  % [ mode byte day dpt NsubDir NsubF SNsubD  SNsubF Sbyte ]
nc = 2;  % { name date }

ini = zeros(0,ni);
cnf =  cell(0,nc);

dpt = dpt + 1;

%------------------------------------------------
% Start

if ~verb
    fprintf(1,'%s',pfad);
end

[ok,name,byte,day,name,date,num] = read_dir(pfad,sm);

if ~verb
    if ~ok
        fprintf(1,' ... can''t read Directory');
    end
    fprintf(1,'\n')
end

n = size(mode,1);

if dpt == 1
   ini = [ 1 0 
if ok
   sum = [ sum(byte)  num ];
else
   sum = zeros(1,3);
end

add = ( dpt == 1 );

%--------------------------------------------------
if n > 0

   ind = ( 1 : n ) + add;

   ini = zeros(n+add,ni); % [ mode byte day dpt NsubDir NsubF SNsubD  SNsubF Sbyte ]
  
   ini(ind,1) = mode;
   ini(ind,2) = byte;
   ini(ind,3) = day;
   ini(ind,4) = dpt;

   cnf = cell(n+add,nc);  % { name date }
   cnf(ind,1) = 
inf = { mode name byte day date num sum };

%--------------------------------------------------
% Recurse

if rec & ( n > 0 )

   id = ( mode == 1 );

   if any(isd)

      nd = sum(isd);
      id = find(isd) ;

      in =  cell(nd,1);
      cn =  cell(nd,1);
      nn = zeros(nd,1);

      ok = zeros(n,1);

      for ii = 1 : nd
          [in{ii},cn{ii}] = get_info([pfad name{ii} fs],dpt,verb,rec,sm,fs);
           nn(ii) = ~isempty(in{ii});
      end  

   sum = sum + cat(1,sub{:,3});
  
   ni = sum(ok);
   jj = find(ok);

   ii = ones( n+ni , 1 );
   jj = jj + ( 1 : ni )';

   ii = cumsum(ii,1);

 %  x = x(ii,:);
 %  x(jj,:) = y;


end


%-------------------------------------------------
% Build FileList

id = ['LFD'];  % Identifer
if strcmp(pfad,'/d1/project/matlab/instruments/'),keyboard,end
if list
%  Dir0/1  Date  Byte  Name 

   mm      = cat(1,real(nn),imag(nn));
   mm(1,:) = sum(mm,1);                 % Total

   Ncell = (1+recurse)*(mm(1,1)+1) + mm(1,2) + 2;
   Nstr  = (1+recurse)*(mm(1,1)+1) + mm(1,2) + 2;

   str      = cell(Nstr,5);
   str(:)   = {''};
   str(:,1) = {' '};
   str(:,5) = {nl};

   nf = sprintf( nform , mform{[1 1]+(mm(2,:) == 0)} );
   ni = ones(2,2);
   ni(2,:) = ~( mm(2,:) == 0 );
   ni = find(ni);

   str(1,:) = [ {'D'}                  dact(2)  ...
                  {sprintf(bform,si)}  dact(4)  ...
                  {sprintf(nf,mm(ni))}              ];


   ll = 1;

   for ii = 1 : size(inf,1)

       ll = ll + 1;

       if isd(ii)  & recurse

          str(ll+(1:size(inf{ii,8},1)),:) = inf{ii,8};

          ll = ll+size(inf{ii,8},1);

          inf{ii,8} = '';

       else

          str(ll,[1 2 3 4]) = [ {id(md(ii)+2)}              ...
                                   inf(ii,2)                  ...
                                  {sprintf(bform,inf{ii,3})}  ...
                                  {cat(2,inf{ii,4},fs(1:is_dir(ii)))}  ];

       end


   end
   % ii


   bb = find( strcmp(' ',str(:,1)) );
   if ~isempty(bb)
       jj = find( diff(bb,1,1) == 1 );
       if ~isempty(jj)
           str( bb(jj)+1 , : )  = [];
       end
   end

   ok = zeros(size(str,1),1);
   for bb = ' LDF'
       ok = ( ok | strcmp( str(:,1) , bb ) );
   end

   str = str(find(ok),:);

end


inf = [ dact(1:4) { inf(:,1:M)  si  nn  str }  ];


if listout | filemode

 VarIn = { file };

 str = dirlist(str,VarIn{1:filemode});

 if ~listout
     str = '';
 end

end

if ~listmode
    inf = inf(:,1:M);
end

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,mode,byte,day,name,date,num] = read_dir(pfad,sm,dpt);

mode = zeros(0,1);
byte = zeros(0,1);
day  = zeros(0,1);
name =  cell(0,1);
date =  cell(0,1);
num  = zeros(1,2);

%-----------------------------------------------------------------

d = dir(pfad);

ok = ~isempty(d)
if ~ok
    num = NaN * num;
    return
end

%-----------------------------------------------------------------

mode = cat(1,d.isdir);
name =      {d.name};

add = ( dpt == 1 );

i1 = strcmp(name,'.');
i0 = strcmp(name,'..');

if dpt == 1
   mode = mode - 2 * i1.*mode;
end

ok = ~( ( mode == 1 ) & ( i0  | i1 );

if ~any(ok)
    mode = zeros(0,1);
    name =  cell(0,1);
    return
end

mode = mode(ok);
name = name(ok);

byte = cat(1,d(ok).bytes);
date =      {d(ok).date};

day = datenum(ger2eng(char(date)));

isd = mode;   % TRUE for Directory

%-----------------------------------------------------------------
% UNIX: Check for Links

cmp     = computer;
if ~any( strcmp( cmp([1 2]) , { 'PC' 'MA' } ) );

    n  = size(mode,1);

    ok = zeros(n,1);    % TRUE for NOT a Link

    for ii = 1 : n
     ok(ii) = ~( unix(sprintf('test -L "%s"',[pfad name{ii}])) == 0 );
    end

    byte = byte .* ok;

    mode = min( mode , 2*ok-1 );

end

%-----------------------------------------------------------------
% Number of: [  Dir+i*DirLink  File+i*FileLink ]

isl = ( mode == -1 );
ild = ( isl &  isd );   % Link to Directory
ilf = ( isl & ~isd );   % Link to File

num = ( sum(mode==1)+i*sum(ild  sum(mode==0)+i*sum(ilf) );

%-----------------------------------------------------------------
% Sort

switch sm

   case 'a'

        [h,si] = sort(name);

   case 's'

        [h,si] = sort(-byte);

   case 't'

        [h,si] = sort(-day);

end

[h,sj] = sort(isd(si));    % Directory last

si = si(sj);

mode = mode(si);
byte = byte(si);
day  =  day(si);
name = name(si);
date = date(si);

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
