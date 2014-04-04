function [ini,ok] = dbcheck(file);

% DBCHECK  Checks Files for Header- and DataVariables
%
% [ INI , Ok ] = DBCHECK( Files )
%
%  INI IdentiferStructure with unique Header- and DataVariables 
%       from Files
%
%  Ok  Matrix with Size of File, 
%      One if File is a DataBase-File, Zero otherwise
%
% see also: DBSCAN, DBREAD, DBFILE
%
 
Nin = nargin;

if Nin < 1
   file = [];
end

ini = dbident;
 ok = [];

if isempty(file)
   return
end

[ok,file] = chkcstr(file);

if ~ok
    error('FileNames must be Strings.');
end

%------------------------------------------------------

cnf = dbget;

cid = dbident(cnf.ColumnKeys);
cid = cid.ID;
   
%------------------------------------------------------

nf = prod(size(file));

ok = zeros(size(file));

dat = ini(ones(nf,1),1);

for ii = 1 : nf
      
    [msg,dat(ii)] = dbread(file{ii},ini,cid,cnf);

    ok(ii) = ~isempty(dat(ii).ID);

end

if ~any(ok)
    return
end

if ~all(ok)
    dat = dat(find(ok),:);
end

fld = fieldnames(ini);

for f = fld(:)'

    ini = setfield( ini , f{1} , eval(['cat(1,dat.' f{1} ')']) );

end

%--------------------------------------------------------
% Check for Duplicate Names

ini = chknam(ini,1);

%--------------------------------------------------------
% Sort by Type

[h,si] = sort(ini.Type);

for f = fld(:)'
    v   = getfield(ini,f{1});
    v   = v(si,:);
    ini = setfield(ini,f{1},v);
end


