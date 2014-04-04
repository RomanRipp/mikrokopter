function cp_list(file,dst);

% CP_LIST  Copy Files from ListFile using CP2FIND
%
% CP_LIST( ListFile , Destination )
%

if ~( exist(dst,'dir') == 7 )
    error('Directory doesn''t exist.')
end

[msg,file] = readfile(file);

if ~isempty(msg)
    error(msg);
end

if ~isempty(file)
    file(find(strcmp(file,''))) = [];
end

if isempty(file)
   return
end

if ~isunix
    error('UNIX-Enviroment required.');
end

for ff = file(:)'

    [p,n,e] = fileparts(ff{1});

    if ~isempty(n)
        name = which(cat(2,n,e));
        if ~isempty(name)
            cp2find(name,dst);
        end
    end

end
