function file = man2txt(page,file)

% MAN2TXT  converts UNIX-MAN-Pages to TextFiles
%
%  FILE = MAN2TXT( PAGE )
%
% calls UNIX: man <Page>
%
%   and writes the TextFile:  FILE = '<PAGE>.man.txt'
%
% MAN2TXT( PAGE , FILE ) specifies the TextFile to write.
%


Nin = nargin;

fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
fcn = upper(fcn);

%**********************************************************
% Check Inputs

if Nin == 0
   error(sprintf('%s: Input PAGE is missing.',fcn));
end

if ~( ischar(page)  & ~isempty(page) & ...
      ( prod(size(page)) == size(page,2) ) )
    error(sprintf('%s: Input PAGE must be a String.',fcn))
end

if Nin < 2
   file = sprintf('%s.man.txt',page);
elseif ~( ischar(file)  & ~isempty(file) & ...
          ( prod(size(file)) == size(file,2) ) )
    error(sprintf('%s: Input FILE must be a String.',fcn))
end

%**********************************************************
% Call MAN

command = sprintf('man %s',page);

fprintf(1,'\nUNIX: %s ',command);

[s,w] = unix(command);

if ~( s == 0 )
    error(sprintf('\n%s: Error call UNIX.\n%s',fcn,w));
end

fprintf(1,'\n');

%**********************************************************
% Remove Bold- and UnderLine Characters

w = strrep(w,char(['_' 8]),'');   % UnderLine

ii = find( w == 8 );

w([ii ii+1]) = [];     % Bold


%**********************************************************
% Write File

fprintf(1,'\nWRITE: %s ',file);

fid = fopen(file,'wt');

if fid == -1
   error(sprintf('\n%s: Can''t open File "%s" for writing.',fcn,file));
end

fwrite(fid,w,'char');

fclose(fid);

fprintf(1,'\n\n');
%**********************************************************

if nargout == 0
   clear file
end
