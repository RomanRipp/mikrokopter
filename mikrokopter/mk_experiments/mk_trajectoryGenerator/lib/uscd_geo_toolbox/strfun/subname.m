function [n,p,ok] = subname(p,def);

% SUBNAME   Split Directory- or FileName in Root and Sub  
%
% [ SubName , RootName , UseDefault ] = SUBNAME( DirName , DefaultSubName )
%
%
% see also: FILEPARTS
%

  Nin = nargin;

  n = '';

  ok = 1;

  if Nin < 1
     p = '';
  end

  if Nin < 2
     def = '';
  end

  if ~( ischar( p )  &  ( prod(size( p )) == size( p ,2) ) &  ...
        ischar(def)  &  ( prod(size(def)) == size(def,2) )         );
     error('Inputs must be Strings.');
  end

  fs = filesep;

   n = def;

  if isempty(p)
     return
  end

  %-------------------------------------
  % Remove Duplicate FileSep

  ii = find( double(fs) == double(p) );

  if ~isempty(ii)

       jj = find( diff(ii) == 1 );

       p(ii(jj+1)) = [];

  end

  if strcmp(p,fs)
     return
  end

  p = p( 1 : ( end - strcmp(p(end),fs) ) );

  %--------------------------------------------
  % Split

  ok = 0;

  ii = find( double(fs) == double(p) );

  if isempty(ii)
     ii = 0;
  else
     ii = max(ii);
  end

  n = p(ii+1:end);
  p = p(1:ii);
