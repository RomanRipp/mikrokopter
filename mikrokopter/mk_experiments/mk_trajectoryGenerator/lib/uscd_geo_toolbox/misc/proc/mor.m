function [s,w] = more(varargin);

% Overload for MORE
%
% MORE   Control paged output in command window.
%    MORE OFF disables paging of the output in the MATLAB command window.
%    MORE ON enables paging of the output in the MATLAB command window.
%    MORE(N) specifies the size of the page to be N lines.
% 
%    When MORE is enabled and output is being paged, advance to the next
%    line of output by hitting the RETURN key;  get the next page of
%    output by hitting the spacebar. Press the "q" key to exit out
%    of displaying the current item.
%
% MORE FILE - file perusal filter for crt viewing
%
% SYNOPSIS
%     more [-dlfpcsu] [-num] [+/ pattern] [+ linenum] [file ...]
%
% DESCRIPTION
%     More is a filter for paging through text one screenful at a time.  This
%     version is especially primitve.  Users should realize that less(1) pro­
%     vides more(1) emulation and extensive enhancements.
%

Nin  = nargin;
Nout = nargout;

fcn = 'more';
alt = 'cat';

f = '';

if Nin > 0
   f  = varargin{Nin};
   ok = chkcstr(varargin);
   if ~ok & ( Nin == 1 )
       ok = ( isnumeric(f) & prod(size(f==1)) );
       if ok
          ok = ( ( mod(f,1) == 0 ) & ( f > 0 ) );
       end
    end
    if ~ok
        error('Inputs must be Strings or a single Integer larger ZERO.');
    end
end

ism = ( ( Nin == 0 ) | ( ( Nin == 1 ) & isnumeric(f) ) );
if ~ism & ( Nin == 1)
    f = lower(f);
    ism = any(strcmp(f,{'on' 'off' }));
end

if ism
   if Nin == 1
      builtin(fcn,f);
   end
   s = builtin(fcn);
   try, w = get(0,fcn); catch, w = ''; end
   if Nout == 0, clear s w, end
   return
end

f  = varargin{Nin};  % Again

if exist(f,'file') == 2
   ff = which(f);
   if ~isempty(ff), varargin{end} = ff; end
end

cmd = sprintf(' %s',varargin{:});

if Nout > 0
   fcn = alt;
end

cmd = sprintf('%s %s',fcn,cmd);

% Nout, cmd, return

if Nout == 0
   unix(cmd);
   clear s w
else
   [s,w] = unix(cmd);
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


