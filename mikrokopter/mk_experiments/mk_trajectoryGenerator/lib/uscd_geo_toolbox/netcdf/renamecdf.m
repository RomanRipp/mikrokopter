function RENAMECDF(file,old,new)

% RENAME_CDF  Renames Variables in a NetCDF-File
%
% RENAMECDF( NetCDF_File , OldName , NewName )
%
% Renames in the NetCDF_File the Variables 
%   wich Names are OldName to NewName
%
% NetCDF_File can be a CharacterArray or CellArray of Strings.
%
% OldName and NewName can be CharacterArrays or CellArrays of Strings.
% OldName and NewName must have the same Number of Elements.
%
%---------------------------------------------------------
%
% RENAMECDF works with 5 basic MEXCDF-Commands:
%
%  fid    = ncmex('OPEN',CDF_File,'write');
%  status = ncmex('REDEF,fid);
%  status = ncmex('VARRENAME',fid,vid,NewName);
%  status = ncmex('ENDEF',fid);
%  status = ncmex('CLOSE',fid);
%
% to get the ID (vid) and the original Name of an Variable, 
%  following MEXCDF-Commands are used:
%
% [ndims,nvars] = ncmex('INQUIRE',fid);
%        name   = ncmex('VARINQ' ,fid,vid);
%
%---------------------------------------------------------
%
% see also: NCMEX
%


Narg = nargin;

if Narg < 3
 error('not enough InputArguments.')
end

%*******************************************************
% Check Inputs

msg = cell(0,1);

%---------------------------------------------------------------
% File

[ok,file] = chkcstr(file);
if ~ok
    msg = cat(1,msg,{'CDF_File must be a StringArray.'});
end

%---------------------------------------------------------------
% OldName

[ok,old] = chkcstr(old);
if ~ok
    msg = cat(1,msg,{'OldName must be a StringArray.'});
else
   old = old(:);
   if any(strcmp(old,''))
      msg = cat(1,msg,{'Empty Strings in OldName.'});
   end
   cc = char(sort(old));
   if any( sum(diff(cc,1,1)==0,2) == size(cc,2) )
      msg = cat(1,msg,{'Duplicate Strings in OldName.'});
   end
end

%---------------------------------------------------------------
% NewName

[ok,new] = chkcstr(new);
if ~ok
    msg = cat(1,msg,{'NewName must be a StringArray.'});
else
   new = new(:);
   if any(strcmp(new,''))
      msg = cat(1,msg,{'Empty Strings in NewName.'});
   end
   cc = char(sort(new));
   if any( sum(diff(cc,1,1)==0,2) == size(cc,2) )
      msg = cat(1,msg,{'Duplicate Strings in NewName.'});
   end
end

%---------------------------------------------------------------

if ~isempty(msg)
    error(cat(1,msg{:}));
end

%---------------------------------------------------------------

if ~( size(old,1) == size(new,1) )
     error('Number of Strings in OldName and NewName must be the same.')
end
   
%*******************************************************

file(find(strcmp(file,''))) = [];

if isempty(file)
   return
end

n = size(old,1);

fprintf(1,'\n');

for ff = file(:)'

    if ~( exist(ff{1},'file') == 2 )
        fprintf(1,'File doesn''t exist: %s',ff{1});
    else

        f = which(ff{1});
        if ~isempty(f)
            ff{1} = f;
        end

        fprintf(1,'Open %s ... ',ff{1});

        [fid,stat] = ncmex('open',ff{1},'write');

        if ~( ( fid > 0 ) & ( stat == 0 ) )
           fprintf(1,'error');
        else

           fprintf(1,'ok');

           status = ncmex('redef',fid);

           if status == -1
              fprintf(1,', Can''t redefine!');
           else

              [nd,nv] = ncmex('INQUIRE',fid);

              % Search for OldName

              vn = cell(nv,1);
              for ii = 1 : nv
                  vn{ii} = ncmex('VARINQ',fid,ii-1);
              end

              iv = zeros(n,1);
              for ii = 1 : n
                  jj = strcmp(vn,old{ii});
                  if any(jj)
                     iv(ii) = find(jj);
                  end
                  jj = strcmp(vn,new{ii});
                  if any(jj)
                     iv(ii) = -1;
                  end
              end
     
              if any( iv == 0 )
                 jj = find( iv == 0 );
                 str = sprintf('      %s\n',old{jj});
                 fprintf(1,'\n   Didn''t find Old Variables:\n%s',str);
              end

              if any( iv == -1 )
                 jj = find( iv == -1 );
                 str = sprintf('      %s\n',new{jj});
                 fprintf(1,'\n   Allready existing New Variables:\n%s',str);
              end

              if all( iv > 0 )
                 fprintf(1,', Rename ... ');
                 for ii = 1 : n
	             status = ncmex('VARRENAME',fid,iv(ii)-1,new{ii});
                     if status == -1
                        break
                     end
                 end
                 if status == -1
                    fprintf(1,'error');
                 else
                    fprintf(1,'ok');
                 end
              end

              ncmex('endef',fid);

           end % REDEF ok

           ncmex('close',fid);

        end % FID ok

    end % FILE exist

    fprintf(1,'\n');

end % FILE

fprintf(1,'\n');

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,str] = chkcstr(str,opt)

% CHKCSTR  Checks Input for CellString, contains Strings !
%
%  [ok,str] = chkcstr(str,Option)
%
%  Option ~= 0 ==> CharacterArrays not allowed,
%
%   default: Option == 0   ==>  CharacterArrays --> CellString
%
 
if nargin < 2
   opt = 0;
end

if strcmp(class(str),'char') & isequal(opt,0)
   n = size(str,1);
   if n == 1
      str = strrep(str,char(32),char(1));
   end
   str = cellstr(str);
   if n == 1
      str = strrep(str,char(1),char(32));
   end
end

ok = iscellstr(str);
if ~ok
   return
end

try
  s = cat(2,str{:});
catch
  ok = 0;
  return
end
 
ok = ( strcmp(class(s),'char')  &  ( prod(size(s)) == size(s,2) ) );


