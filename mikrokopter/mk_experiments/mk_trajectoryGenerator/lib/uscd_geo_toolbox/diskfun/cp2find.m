function [f,ok] = cp2find(src,dst,Target);

% CP2FIND   Copy Files recursive to Destination
%
%-------------------------------------------------------------
%
% CP2FIND( Name , Destination );
%
%    Copy File with Name to same Files, found in Directory Destination
%
% CP2FIND uses the File with Name which is found in Matlab's SearchPath.
%
%-------------------------------------------------------------
%
% CP2FIND( File , Destination , Target );
%
%    Copy File to Directory of Targets, found in Destination
%
%-------------------------------------------------------------
%
% [Files,Ok] = CP2FIND( ... )
%
% Returns a list of Files which are found in Destination. 
% The corresponding Vector Ok is true if COPYFILE was successfull.
%
%-------------------------------------------------------------
%
% see also: WHICHFILE (required), UPDATE, COPYFILE
%

Nout = nargout;

f  = {};
ok = [];

src0 = src;

[pfd,name,ext] = fileparts(src);

if isempty(pfd)
   if exist(src,'file') == 2
      f = which(src);
      if ~isempty(f)
          src = f; 
         [pfd,name,ext] = fileparts(src);
      end
   else
      if ~strcmp(src(1),filesep)
          src = [ pwd filesep src ];
      end
  end
end

if ~( exist(src,'file') == 2 )
    d = dir(src);
    if ~( prod(size(d)) == 1 )
        error(['File not found: ' src0 ]);
    end
end

if ~( exist(dst,'dir') == 7 )
    error(['Destination not found: ' dst ]);
end

if strcmp(dst,'.')
   dst = pwd;
end

if strcmp(dst,'..')
   dst = pwd;
   dst = dst(1:max(findstr(dst,filesep))-1);
end

if nargin < 3
   Target = [ name ext ];
end


fprintf(1,'\nSearch in "%s" for "%s" ...',dst,Target);

msg = '';
try
   [n,p,f] = whichfile({dst},Target,'-e','-r');
catch
   msg = lasterr;
end

if ~isempty(msg)
    fprintf(1,'%s error\n%s\n',char(7),msg);
   if Nout == 0, clear('f'), end
    return
end

fprintf(1,' ok\n');

if isempty(f)
   fprintf(1,'%sNo targets found in "%s"\n\n',char(7),dst);
   if Nout == 0, clear('f'), end
   return
end

ii = ~strcmp(f,src);
if ~all(ii)
    if ~any(ii)
        ii = [];
    else
        ii = find(ii);
    end
    n = n(ii);
    p = p(ii);
    f = f(ii);
end

if isempty(f)
   fprintf(1,'%sNo other targets found in "%s"\n\n',char(7),dst);
   if Nout == 0, clear('f'), end
   return
end

n = size(f,1);

ok = zeros(n,1);

for ii = 1 : n

       fprintf(1,'Copy: %s %s',src,p{ii});

       [ok(ii),msg,id] = copyfile(src,p{ii});

       if ~ok(ii)

           fprintf(1,'%s error\n%s\n',char(7),msg);

       end

       fprintf(1,'\n');
     
end

if Nout == 0
   clear('f')
end
