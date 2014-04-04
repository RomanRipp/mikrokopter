function ps(varargin)

% PS  Unix ps  [options] -u USER
%
% PS( options )
%
 
if nargin > 0
   if ~iscellstr(varargin)
      error('Options must be Strings');
   end
end

v = strhcat(varargin,' ');

%------------------------------------
% Get User

w = getenv('USER');

if isempty(w)

  command = 'whoami';

  [s,w] = unix(command);

  if ~( s == 0 )
    msg = ['Error call UNIX: ' command ];
    if ~isempty(w)
       msg = [ msg char(10) w ];
    end
    error(msg);
  end

end

w = w( 1 : ( end - 1*( double(w(end)) == 10 ) ) );

%-------------------------------------------------------

command = ['ps ' v  ' -u ' w  ' --format "%U%p  %P%t   %c  %a" '];

fprintf(1,'\nUNIX: %s\n\n',command);

[s,w] = unix(command);

if ~isempty(w)
  fprintf('\n%s\n\n',w);
end

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  str = strhcat(str,del,n,nl)

% STRHCAT  Concatenates Strings into ONE
%
% STRHCAT( StringArray , Delimiter )
%   Forms one long String from the Strings in the
%   StringArray, delimited with the delimiter.
%   The EndDelimiter will be removed.
%
% STRHCAT( StringArray , Delimiter , N , NewLine )
%   Build a  NewLine after each N-th String.
%   default: N = 10;  NewLine = char(10);
%
% Example:  
%         >> strhcat({'apples' 'pies' 'prunes'},', ')
%    
%         ans =
%
%         apples, pies, prunes
%
%         >> strhcat({'apples';'pies';'prunes'},', ',2)
%    
%         ans =
%
%         apples, pies
%         prunes
%



Nin = nargin;

if Nin < 4
 nl = char(10);
end
if Nin < 3
 n = 10;
end
if Nin < 2
 del = char((32*ones(1,3)));
end


if isempty(str)
 str = '';
 return
end


if ischar(str)
  str = cellstr(str);
end

str = str(:);

str( : , 2 ) = { del };

str(n:n:size(str,1),2) = {nl};

str(    size(str,1),2) = {''};


str = str';

str = cat(2,str{:});


