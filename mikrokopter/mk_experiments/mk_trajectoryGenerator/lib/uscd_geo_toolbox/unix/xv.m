function xv(file)

% XV  Open's XV Image Viewer
%
% XV( FileName )
%
    

if nargin == 0

  file = '';

end

if ~( ischar(file)  &  ( prod(size(file)) == size(file,2) ) )
  error('XV: Input FileName must be a String.')
end

if ~isempty(file)

  file0 = file;

  d = dir(file);

  if isempty(d)

    file = cat(2,file,'*');

    d = dir(file);

  end


  if isempty(d)

    fprintf([char(10) 'XV: File ' file0 ' not found.' char([10 10]) ]);

    return

  end

end

command = ['xv ' file  ' &' ];

fprintf([ char(10) 'UNIX: '  command char([10 10]) ]);

%%% msg = lib_kim(0);

[s,w] = unix(command);

%%% msg = lib_kim(1);

if ~isempty(w)
  fprintf([ char(10) w char([10 10]) ]);
end
