function medit(file)

% MEDIT  Open's MCEDIT to editing a file
%
% MEDIT( FileName )
%

file0 = which(file);
if isempty(file0)
  file0 = file;
  jj = findstr(file0,'.');
  if isempty(jj)
     file0 = [ file0 '.m' ];
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

 
jj = [ 0  findstr(file0,'/') ];
name = file0(jj(end)+1:end);


geometry = [ 80  40 ];

fid = fopen(file0,'r');

if fid ~= -1

  bb = fread(fid,'char');

  fclose(fid);

  bb = cat(1,10,bb,10);
  ww = max(diff(find(bb==10)))+10;

  geometry(1) = min(120,max(80,ww));

end

geometry = sprintf('%.0fx%.0f',geometry);

command = ['xterm -geometry '  geometry  '-10+10  -title "'  name '"'   ...
           ' -exec mcedit -c -x '  file0  ' & '];

fprintf([ char(10) 'UNIX: '  command char([10 10]) ]);

[s,w] = unix(command);

if ~isempty(w)
  fprintf([ char(10) w char([10 10]) ]);
end
