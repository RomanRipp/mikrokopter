function msg = lib_kim(mode)

% LIB_KIM  Renames Links of Matlab-LIB-files
%
% required to access X using Matlab5 on SuSE 9, HOST kim
%
% BINDIR: /opt/matlab53/bin/lnx86
%
% Changes: libXt.so.6  <---> libXt.so.6.tmp
%          libX11.so.6 <---> libX11.so.6.tmp
%
% LIB_KIM( Mode )
%
% Mode:    0   --->
%          1  <---
%
% 

host = { 'kim' 'kimm' };
vers = '5';

bin = '/opt/matlab53/bin/lnx86';

lnk = { 'libX11.so.6' 
        'libXt.so.6'  };

suf = '.tmp';

n = size(lnk,1);

%**************************************************************

if nargin == 0
   error('Mode required.');
end

mode = isequal(mode,1);

%**************************************************************
% Check for HOST, VERSION, BINDIR, FILES

msg = cell(0,1);

if ~any(strcmp(getenv('HOSTNAME'),host))
    msg = cat(1,msg,{sprintf('HOSTNAME must be "%s"',host)});
end

v = version;

if ~( v(1) == vers )
    msg = cat(1,msg,{sprintf('VERSION must be "%s"',vers)});
end

if ~( exist(bin,'dir') == 7 )
    msg = cat(1,msg,{sprintf('Directory doesn''t exist: %s',bin)});
else
   for ii = 1 : n
       lnk{ii} = fullfile(bin,lnk{ii});
       f = cat( 2 , lnk{ii} , suf( 1 : end*mode ) );
       if ~( exist(f,'file') == 2 )
           msg = cat(1,msg,{sprintf('File doesn''t exist: %s',f)});
       end
   end
end

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    if nargout == 0
       warning(msg);
       clear msg
    end
    return
end

%**************************************************************

for ii = 1 : n

    src = cat( 2 , lnk{ii} , suf( 1 : end*mode ) );
    dst = cat( 2 , lnk{ii} , suf( 1 : end*(~mode) ) );

    cmd = sprintf('mv %s %s',src,dst);

    [s,w] = unix(cmd);

    if ~( s == 0 )
        msg = cat(1,msg,{sprintf('Error call UNIX: %s\n%s',cmd,w)});
    end

end

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    if nargout == 0
       warning(msg);
    end
else
    msg = '';
end

if nargout == 0
   clear msg
end

