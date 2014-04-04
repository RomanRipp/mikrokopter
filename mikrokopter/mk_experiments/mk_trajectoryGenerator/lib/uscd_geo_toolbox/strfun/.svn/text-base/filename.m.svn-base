function [msg,pfad,name] = filename(file,ext)


% FILENAME  Get's Full FileNmae of File
%
% [ Msg , FullName , Name ] = FILENAME( File , Extension )
%
%  If the 2. Input "Extension" is given, for a "File.Extension" will tried.
%


msg  = '';
msg0 = 'FILENAME: ';

fullname = '';
    name = '';

      fs = '.';


%---------------------------------------------------------------

if nargin == 0
  return
end

if isempty(file)
   return
end

%---------------------------------------------------------------
% Check File

if ~( ischar(file)  &  ( prod(size(file)) == size(file,2) ) )
  msg = [ msg0 'Input FileName must be a String.'];
  return
end

%---------------------------------------------------------------
% Check Extension

if nargin < 2
   ext = char(ones(1,0));
elseif isempty(ext)
   ext = char(ones(1,0));
end

if ~( ischar(ext)  &  ( prod(size(ext)) == size(ext,2) ) )
  msg = [ msg0 'Input Extension must be a String.'];
  return
end

%---------------------------------------------------------------

fullname = GetWhichFileName(file);

if isequal(fullname,0)
   fullname = '';
   msg = [ msg0 'A Directory exist: ' file ];
   return
end


if ~isempty(ext)

   ext = ext( (1+strcmp(ext(1),fs)) : end );

   new = isempty(fullname);

   if new

      % Check File with Extension
  
      file = file( 1 : (end-1*strcmp(file(end),fs)) );

   else

      new = ( strcmp(ext,'m') & ( exist(fullname) == 6 ) )

      if new

        % Save before P-Code

        jj = find( double(fullname) == double(fs) );

        if ~isempty(jj)
           file = fullname( 1 : (max(jj)-1) ) ;
        else
           file = fullname;
        end

      end

   end   

   if new

     file = cat( 2 , file , fs , ext );

     fullname = GetWhichFileName(file);

     if isequal(fullname,0)
        fullname = '';
        return
     end

   end

end

%----------------------------------------------------------

if nargout > 2

    fs = filesep;

    if ~isempty(fullname)
        file = fullname; 
    end

    jj = find( double(file) == double(file) );

    if isempty(jj)
       name = file;
       return
    end

    jj = max(jj);

    p0 = cd;

    try
      cd(file(1:jj))
      p1 = cd;
      if isempty(p1) & strcmp(p1,fs)
         name = fullname( (max(jj)+1) : end ) ;
      end 
    catch
         name = fullname; 
    end

end


%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [FullWhichPfadName,FullWhichFileName] = GetWhichFileName(NoWhichFileName)

FullWhichPfadName = '';
FullWhichFileName = which(NoWhichFileName);

bad = { 'variable'  'built-in' };

if any(strcmp(FullWhichFileName,bad)) & ~any(strcmp(NoWhichFileName,bad))

   FullWhichFileName = '';

   return

end

if ~isempty(FullWhichFileName)

   return

end

%----------------------------------------------

d = dir(NoWhichFileName);

if isempty(d)
   return
end


bad = ( prod(size(d)) > 1 );
if ~bad
    bad = d.isdir;
end

if bad
   FullWhichFileName = 0;
   return
end

fs = filesep;

jj = find( double(NoWhichFileName) == double(fs) );

if isempty(jj)
   return
end

jj = max(jj);

p0 = cd;

try
  cd(NoWhichFileName(1:jj));
  FullWhichFileName = cd;
catch
  return
end

cd(p0);

if isempty(FullWhichFileName)
   FullWhichFileName = fs;
else
   FullWhichFileName = cat( 2 , FullWhichFileName , ...
                                fs(1:(end*(~strcmp(FullWhichFileName(end),fs)))) );
end

FullWhichFileName = cat( 2 , FullWhichFileName , ...
                             NoWhichFileName( (jj+1) : end )  );


