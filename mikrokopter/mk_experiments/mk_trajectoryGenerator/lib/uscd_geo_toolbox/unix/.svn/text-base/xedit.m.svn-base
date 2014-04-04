function xedit(file)

% XEDIT  Open's XEDIT to editing a file
%
% XEDIT( FileName )
%


if nargin == 0
  file = '';
end

if ~( ischar(file)  &  ( prod(size(file)) == size(file,2) ) )
  error('XEDIT: Input FileName must be a String.')
end

if ~isempty(file)

  file0 = which(file);
  if isempty(file0) 
     file0 = file;
     if ~( exist(file,'file') == 2 )
        jj = findstr(file0,'.');
        if isempty(jj)
           file0 = [ file0 '.m' ];
        end
     end
  end

  if exist(file0) == 6
  % P-Code

    [p,n] = fileparts(file0);

    file1 = [ fullfile(p,n) '.m' ];
   
    if exist(file1) == 2
      file0 = file1;
    end

  end

end

jj = [ 0  findstr(file0,'/') ];
name = file0(jj(end)+1:end);

command = ['xedit -geometry 500x600-10+10  -title "'  name '" ' file0 ' & '];

fprintf([ char(10) 'UNIX: '  command char([10 10]) ]);

%%% msg = lib_kim(0);

[s,w] = unix(command);

%%% msg = lib_kim(1);

if ~isempty(w)
  fprintf([ char(10) w char([10 10]) ]);
end
