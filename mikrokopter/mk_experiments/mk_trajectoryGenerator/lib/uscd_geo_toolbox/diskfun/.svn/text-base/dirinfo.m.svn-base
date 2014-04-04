function [out,liste] = dirinfo(varargin)

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


out   = cell(0,M+listmode);
liste = '';

% Format for Bytes
bform = '%10.0f';  

% Format for [ NDir NFiles ]  
mform = { '#%.0f(%0.f)'  '#%.0f' };  % { InclLink    NoLink }
nform = '%s   %s';                   %   'DirForm  FileForm'


%------------------------------------------------
% Start

fprintf(strrep(pfad,'\','\\'));

d = dir(pfad);

N = size(d,1);


if N == 0
  fprintf([' ... can''t read Directory ' nl ]);
  return
end

fprintf(nl)

out = cell(N,M+1);
out(:,5) = { cell(0,M) };
out(:,6) = { 0  };
out(:,7) = { [ 0  0 ] };
out(:,8) = { { '' '' '' '' '' } };

cmp     = computer;
is_unix = ~any( strcmp( cmp([1 2]) , { 'PC' 'MA' } ) );

is_dir = cat(1,d.isdir);

[name,si] = sort({d.name});

[is_dir,sj] = sort(is_dir(si));

for ii = 1 : N
    jj = si(sj(ii));
    out{ii,1} = d(jj).isdir;
    out{ii,2} = d(jj).date; 
    out{ii,3} = d(jj).bytes;
    out{ii,4} = [ pfad  d(jj).name ];
    if is_unix
       isl = ( unix(sprintf('test -L "%s"',out{ii,4})) == 0 );
       out{ii,1} = min( out{ii,1} , 1-2*isl );
       out{ii,3} = out{ii,3} * ( 1 - isl );
    end
end


isact = find( is_dir & strcmp(out(:,4),[pfad '.' ]) );
ispre = find( is_dir & strcmp(out(:,4),[pfad '..']) );

dact      = out(isact,:);
dact{1,4} = pfad;

 
out([isact,ispre],:) = [];

is_dir([isact,ispre]) = [];

nn = dact{1,7};  % [ 0  0 ]
si = dact{1,6};  %  0

%--------------------------------------------------
% Recurse

if size(out,1) > 0

   VarIn = { '-r'  '-l' };

   md = cat(1,out{:,1});

   isd = ( md == 1 );
   isf = ( md == 0 );

   isl = ( md == -1 );
   ild = ( isl &  is_dir );
   ilf = ( isl & ~is_dir );

   if recurse

      for ii = find(isd')
          out(ii,1:M+list) = dirinfo(out{ii,4},VarIn{1:1+list});
      end  

   end   
  
   si = sum( ( cat(1,out{:,3}) + cat(1,out{:,6}) ) , 1 );
   
   nn = [ sum(isd)+i*sum(ild) sum(isf)+i*sum(ilf) ] + sum(cat(1,out{:,7}));

end


%-------------------------------------------------
% Build FileList

id = ['LFD'];  % Identifer

if list
%  Dir0/1  Date  Byte  Name 

   mm      = cat(1,real(nn),imag(nn));
   mm(1,:) = sum(mm,1);                 % Total

   Nlist = (1+recurse)*(mm(1,1)+1) + mm(1,2) + 2;

   liste      = cell(Nlist,5);
   liste(:)   = {''};
   liste(:,1) = {' '};

   nf = sprintf( nform , mform{[1 1]+(mm(2,:) == 0)} );
   ni = ones(2,2);
   ni(2,:) = ~( mm(2,:) == 0 );
   ni = find(ni);

   liste(1,:) = [ {'D'}                dact(2)  ...
                  {sprintf(bform,si)}  dact(4)  ...
                  {sprintf(nf,mm(ni))}              ];


   ll = 1;

   for ii = 1 : size(out,1)

       ll = ll + 1;

       if isd(ii)  & recurse

          liste(ll+(1:size(out{ii,8},1)),:) = out{ii,8};

          ll = ll+size(out{ii,8},1);

          out{ii,8} = '';

       else

          liste(ll,[1 2 3 4]) = [ {id(md(ii)+2)}              ...
                                   out(ii,2)                  ...
                                  {sprintf(bform,out{ii,3})}  ...
                                  {cat(2,out{ii,4},fs(1:is_dir(ii)))}  ];

       end


   end
   % ii


   bb = find( strcmp(' ',liste(:,1)) );
   if ~isempty(bb)
       jj = find( diff(bb,1,1) == 1 );
       if ~isempty(jj)
           liste( bb(jj)+1 , : )  = [];
       end
   end

   ok = zeros(size(liste,1),1);
   for bb = ' LDF'
       ok = ( ok | strcmp( liste(:,1) , bb ) );
   end

   liste = liste(find(ok),:);

end


out = [ dact(1:4) { out(:,1:M)  si  nn  liste }  ];


if listout | filemode

 VarIn = { file };

 liste = dirlist(liste,VarIn{1:filemode});

 if ~listout
   liste = '';
 end

end

if ~listmode
  out = out(:,1:M);
end

