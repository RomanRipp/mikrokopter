function c = dirhist(pfad,mode)

% DIRHIST  returns Directory-History
%
% History = DIRHIST( Pfad , SortMode )
%
%  for SortMode see DIRCONT
%

 c = cell(0,2);

fs = filesep;


%--------------------------------------------

if nargin < 1
  pfad = cd;
else
  if ~( ischar(pfad) & ( prod(size(pfad)) == size(pfad,2) ) )
    error('Input must be a String.');
  end
  if isempty(pfad)
    pfad = fs;
  end
end
 
if nargin < 2
  mode = { 'a' };  % Alphabetical
end

%--------------------------------------------
% Try to get absolute Path

if isempty(pfad)
  pfad = cd;
end


p0 = cd;
if isempty(p0)
   p0 = fs;
end

p0 = cat( 2, p0 , fs(1:(end*(~strcmp(p0(end),fs)))) , pfad );


d0 = dir(p0);

d1 = dir(pfad);

if isequal( d0 , d1 )

   pfad = p0;

end

if isempty( d1 )
  c = get_disk;
  return
end

pfad = cat( 2 , pfad , fs(1:(end*(~strcmp(pfad(end),fs)))) );


%---------------------------------------------------------
% Remove duplicate FileSep


if size(pfad,2) > 1
   is_fs = ( double(pfad) == double(fs) );
   bad = find( ( diff(is_fs) == 0 ) & ( is_fs(2:end) == 1 ) );
   pfad(bad) = [];
end

%---------------------------------------------------------
% Remove '.'

pfad = strrep( pfad , [ fs '.' fs ] , fs );

%---------------------------------------------------------
% Remove '..'

bc = [ fs '..' fs ];
sb = size(bc,2);

ok = 0;

while ~ok

  is_bc = findstr( pfad , bc );

  ok = ( isempty( is_bc )  | ( size(pfad,2) < sb ) );

  if ~ok

     is_bc = is_bc(1);

     is_fs = find( double(pfad) == double(fs) );

     n = sum( is_fs < is_bc );

     pfad = cat( 2 , pfad(1:is_fs(n+(n==0))) , pfad(is_bc+sb:end) );

  end

end 

%---------------------------------------------------------
% History

is_fs = find( double(pfad) == double(fs) );

nfs = prod(size(is_fs));

%---------------------------------------------------------
% Contents of Directory

cp = dircont( pfad , {mode} , fs );

bad = find( strcmp( cp(:,4) , '.' )  |  strcmp( cp(:,4) , '..' ) );

cp(bad,:) = [];

np = size(cp,1);

c = cell( np+nfs , 2 );


%---------------------------------------------------------

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if strcmp(upper(computer),'PCWIN')

   pfad = lower( pfad );

end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%---------------------------------------------------------
% Directories on same Level

for ii = 1 : np

  c{ii,1} = cat( 2 , pfad , cp{ii,4} );

  c{ii,2} = 1;

end

%---------------------------------------------------------
% Parent Directories

for ii = 1 : nfs

  iic = nfs - ii + 1 + np;

  c{iic,1} = pfad( 1 : ( is_fs(ii) - 1 * ( is_fs(ii) > 1 ) ) );

  c{iic,2} = -1 * ( ii < nfs ) ;

end


c = cat( 1 , c , get_disk(c{end,1}) );


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function c = get_disk(pfad)

% GET_DISK returns Drives of PCWIN-System

c = cell(0,2);

if ~strcmp( upper(computer) , 'PCWIN' )

  return

end

if nargin < 1
   pfad = '';
end


fs = filesep;

np = size(pfad,2);

if ~isempty(pfad)
   d0 = dir(pfad);
end

disk = 'cdefgijhklmnop';
ds   = ':';       % DiskSep

nd = size(disk,2);

c = cell(nd,2);

c(:,2) = { 0 };

for ii = 1 : nd

    c{ii,1} = cat(2,disk(ii),ds);
    
    try
      c{ii,2} = -2*(~isempty(dir(c{ii,1})));
    end
  
    nc = size(c{ii,1});

    if ~isempty(pfad) 

       ok = ( np >= nc );
       if ok
          ok = strcmp( lower(pfad(1,1:nc)) , c{ii,1} ); 
       end
 
       if ~ok
           ok = isequal( d0 , dir(c{ii,1}) );
       end

       c{ii,2} = c{ii,2} * (~ok);

    end

end

c = c ( find( cat(1,c{:,2}) == -2 ) , : );

