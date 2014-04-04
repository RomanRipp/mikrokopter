function dos2unix(varargin)

% DOS2UNIX  Converts CR/LF to LF
%
% dos2unix(Pfad,'-r','.ext')
% dos2unix(Pfad,'-r','.ext')
%

vin = varargin;

if ~chkcstr(vin)
    error('Inputs must be Strings.');
end

vin = vin(:);
Nin = size(vin,1);

recurse = 0;
for ii = 1 : Nin
    recurse = ( recurse | strcmp(vin{ii},'-r') );
end

ext = cell(1,0);
for ii = 2 : Nin
  if strcmp(vin{ii}(1),'.');
    ext = cat(2,ext,vin(ii));;
  end
end

if isempty(ext)
   ext = {'.m'};
end

Pfad = '';
if Nin >= 1
  if ~strcmp(vin{1},'-r');
    Pfad = vin{1};
  end
end

if any(strcmp(Pfad,{ '.'  '*' }))
   Pfad = '';
end

if isempty(Pfad)
   Pfad = cd;
end


%****************************************************************
% Check for File or Directory

fprintf(1,'%s',Pfad);

IsFile = ( exist(Pfad,'file') == 2 );
if IsFile
   f = which(Pfad);
   if ~isempty(f)
       Pfad = f;
   end
else
   IsDir = ( exist(Pfad,'dir') == 7 );
   if IsDir
      p0 = cd;
      try
         cd(Pfad);
         Pfad = cd;
      catch
         IsDir = 0;
      end
      cd(p0)
   end
   if IsDir
      fs = filesep;
      Pfad = cat( 2 , Pfad , fs(1:(end*(~strcmp(Pfad(end),fs)))) );
   else
      fprintf(1,' ... File or Directory doesn''t exist\n');
      return
   end
end

%****************************************************************
% Get DirectoryContents

d = dir(Pfad);

if isempty(d) 
  fprintf(1,' ... can''t read File or Directory\n');
  return
end

fprintf(1,'\n');

%****************************************************************
% Translate
is_win  = strcmp(upper(computer),'PCWIN');

is_file = find(~cat(1,d.isdir));

if is_win
   ext = lower(ext);
end

for ii = is_file(:)'

  name = d(ii).name;

  % Check for Extension

  [p,n,e] = fileparts(name);

  if is_win
     e = lower(e);
  end

  if any(strcmp(e,ext))

      if IsFile
         file = Pfad;
      else
         file = [ Pfad  name ];
      end

      fid = fopen(file,'r');
      if fid == -1
         fprintf(1,'   Error read "%s"\n',file);
      else
         bb = fread(fid,'char');
         fclose(fid);
         if any(bb==0)
            fprintf(1,'   Invalid Character in "%s"\n',file);
         else
            bb = strrep(char(bb(:)'),char([13 10]),char(10));
            bb = strrep(char(bb(:)'),'%','%%');
            bb = strrep(char(bb(:)'),'\','\\');
           fid = fopen(file,'wt');
           if fid == -1
              fprintf(1,'   Error open to write "%s"\n',file);
           else
              fprintf(fid,bb);
              fclose(fid);
           end
         end
      end 

  end

end

%******************************************************
% Recurse

if recurse

  is_dir = find(cat(1,d.isdir));

  for ii = is_dir(:)'

    if ~any(strcmp(d(ii).name,{ '.'  '..' }))

        dos2unix([Pfad d(ii).name],'-r',ext{:});

    end

  end

end      

%*******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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



   

