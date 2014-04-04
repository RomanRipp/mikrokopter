function dbsave(ident,varargin)

%  'H2' , H1 , 'H2' , H2 , ... , 'D1' , D1 , 'D2' , D2 , ...
%
%  'H2' , H1 , 'H2' , H2 , ... , D1 , D2 , ...
%
%  'H2' , H1 , 'H2' , H2 , ... , Data , ...
%
%  HeadStruct , Data
%  HeadStruct , DataStruct
%
%  HeadStruct , 'D1' , D1 , 'D2' , D2 , ...
%  HeadStruct , D1 , D2 , ...
%
%  HeadDataStruct  
%
%  IdentStruct:   [ID], Name, [Format], [Length], [Type], Value
%
%  DI:   [ NData by NFiles ] | [ NData by 1 ]
%
%  Data: [ NData by NPar by NFiles ]
%

Nin = nargin;

if Nin < 2
   error('Not enough InputArguments.');
end

[file,msg,ok] = dbfile(ident,1);

if ~isempty(msg) 
    dbmsg(msg,1);
elseif isempty(file)
    msg = 'No valid FileNames.');
    return
end

[msg,v] = structchk(varargin);  % HeaderDataStruct

ok = isempty(msg);
if ~ok
    [msg,v] = structchk(varargin{1}); % HeaderStruct
    ok = isempty(msg);
    if ok
       varargin = varargin{2:end};
    end
end


if ( Nin == 1 )
    
   if ~isstruct(varargin{1})
       error('Single Input must be a Structure');
   end
 
%
% write data matrix
%

eval(['data = v',num2str(i + 1),';'])

fprintf(fp,format,data');


% close output file
fclose(fp);





