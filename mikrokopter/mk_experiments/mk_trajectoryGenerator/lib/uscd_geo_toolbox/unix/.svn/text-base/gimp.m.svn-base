function gimp(file)

% GIMP  Open's GIMP Image Manipulating Tool
%
% GIMP( FileName )
%
    

if nargin == 0

  file = '';

end

if ~isempty(file)

  file0 = file;

  d = dir(file);

  if isempty(d)

    file = cat(2,file,'*');

    d = dir(file);

  end


  if isempty(d)

    fprintf([char(10) 'GIMP: File ' file0 ' not found.' char([10 10]) ]);

    return
 
 end

end


command = ['gimp ' file  ' &' ];

fprintf([ char(10) 'UNIX: '  command char([10 10]) ]);

%%% msg = lib_kim(0);

[s,w] = unix(command);

%%% msg = lib_kim(1);

if ~isempty(w)
  fprintf([ char(10) w char([10 10]) ]);
end
