function [s,w] = kim(varargin)

% KIM executes UnixCommands on HOST kim, using LIB_KIM and UNIX
%
% [S,W] = KIM( Command , [Parameter] )
%
% required to access X under Matlab5, SuSE 9
%
% see also: LIB_KIM, UNIX
%

Nout = nargout;

if ~chkcstr(varargin)
    error('Inputs must be Strings.');
end

cmd = sprintf('%s ',varargin{:});

%%% msg = lib_kim(0);

[s,w] = unix(cmd);

%%% msg = lib_kim(1);

if ( nargout < 2 ) & ~isempty(w)
   if s == 0
      fprintf(1,'%s\n',w);
   else
      warning(w);
   end
end

if nargout == 0
   clear s
end

%**************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%**************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function msg = lib_kim(mode)

% LIB_KIM  Renames Links of Matlab-LIB-files
%
% Nessecary to open an Window under X using Matlab5 on SuSE9
%
% Works on HOST kim, VERSION 5, BINDIR /opt/matlab53/bin/lnx86
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

host = 'kim';
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

if ~strcmp(getenv('HOSTNAME'),host)
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

