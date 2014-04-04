function [p,field,f] = methods(obj);

% METHODS   Returns Methods of Object
%
% MTD = Method( Object )
%
% MTD = { Class  {Methods} }  % 2-Column-CellArray
%
% [ MTD , FieldNames , Index ] = METHODS( Object )
%
% Returns Index, that you get corressponding: cat(1,MTD{:,2}) <--> MTD(Index,1) 
%

Nout = nargout;

[p,field] = history(obj);

if Nout >= 2
   field = cat(1,field{:});
end

%------------------------------------
% Get Methods

n = size(p,1);

m = cell(n,1);
f = cell(n,1);

for ii = 1 : n

     mm = methods(p{ii});

     % Remove Creator-Method
 
     mm(find(strcmp(mm,p{ii}))) = [];

     m{ii} = mm;
     
     f{ii} = ii*ones(size(mm,1),1);

end

%------------------------------------
% Sort Methods, remove duplicate

m = cat(1,m{:});
f = cat(1,f{:});

if isempty(m)
   p = cell(0,2);
   return
end

[h,si] = sort(m);

h = char(h);

bad = find( sum( diff(h,1,1) == 0 , 2 ) == size(h,2) );

m(si(bad)) = [];
f(si(bad)) = [];

%------------------------------------
% Sort Methods, remove duplicate

ok = zeros(n,1);
p  = p(:,[1 1]);

for ii = 1 : n
    
    kk = find( f == ii );

    ok(ii) = ~isempty(kk);

    if ok(ii)
       p{ii,2} = m(kk);
    end

end

p = p(find(ok),:);

if ~( Nout == 0 )
    return
end

%**************************************+
% Display

fs = filesep;
cs = '@';
nl = char(10);

np = size(p,1);

str =  cell(np,1);

for ii = 1 : np

     txt = p{ii,2}(:,[1 1]);
     txt(:,1) = {cat(2,cs,p{ii,1})};

     txt = permute(txt,[2 1]);

     txt = strhcat(txt,fs,2);

     str{ii} = cat(2,txt,nl);

end

str = strhcat(str,'',1);

fprintf(1,'\nMethods for %s-Object:\n',upper(class(obj)));
fprintf(1,'\n%s\n',str);

clear p field f

