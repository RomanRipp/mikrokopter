function [p,f] = history(obj);

% HISTORY  resturns the ClassHistory and FieldNames of an Object
%
% [ ClassHistory , FieldNameHistory ] = HISTORY( object )
%
%  calls:  @class/fieldnames.m
%          @class/subsref.m
%          @class/history.m
%
 
p = cell(1,1);
f = cell(1,1);

p{1} = class(obj);

[f{1},par,strct] = fieldnames(obj);

if ~isempty(par)

    s = struct( 'type' , { '.'  } , ...
                'subs' , { par  }      );

  obj = subsref( strct , s );

  [p1,f1] = history(obj);

  p = cat(1,p1,p);
  f = cat(1,f1,f);

end

if ~( nargout == 0 )
    return
end

%**************************************+
% Display

cs = '@';
bl = char(32*ones(1,3));
sp = char(32*ones(size(cs)));
nl = char(10);

np = size(p,1);

str =  cell(np,1);

for ii = 1 : np

    % @class
    pre = cat(2,cs,p{ii});
         
    txt = f{ii}(:,[1 1]);
    txt(:,1) = {bl};

    txt = permute(txt,[2 1]);

    txt = strhcat(txt,' ',2);

    str{ii} = cat(2,pre,nl,nl,txt,nl);

end

str = strhcat(str,'',1);

fprintf(1,'\nHistory for %s-Object:\n',upper(class(obj)));
fprintf(1,'\n%s\n',str);

clear p f

