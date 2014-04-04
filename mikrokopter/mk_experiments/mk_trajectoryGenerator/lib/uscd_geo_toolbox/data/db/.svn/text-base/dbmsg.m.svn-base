function dbmsg(msg,mode,cnf);

% DBMSG  display Messages in Command or LogFile
%
% DBMSG( Message , Mode )
%
% Mode == 1 writes the DateTime and FunctionStack
%           into the LogFile.
% 

Nin = nargin;

if Nin < 2
   mode = 0;
end

if Nin < 3
   cnf = [];
end


if isempty(cnf)
   cnf = dbget;
end

ism = ~isempty(msg);

if ism & strcmp(cnf.Verbose,'on')
   fprintf(1,'%s',msg);
end

if strcmp(cnf.LogMode,'off')
   return
end

[f,s] = fopen(cnf.LogID);

if isempty(f)
   return
end

if ( s(1) == 'a' )
   if isequal(mode,1)
       c = caller;
       if isempty(c)
          m = '';
       else
          c = upper(c(:,1));
          m = sprintf('%s:',c{:});
          m = m(1:end-1);
       end
       c = clock;
       c = datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6)));
      fprintf(cnf.LogID,'>> %s %s\n',c,m);
   end
   if ism
      fprintf(cnf.LogID,'%s',msg);
   end
end