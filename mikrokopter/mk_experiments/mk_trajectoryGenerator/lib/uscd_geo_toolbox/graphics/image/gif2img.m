function msg = gif2img(gif,bgc,pfad);

% GIF2IMG   Extracts Frames from animated GIF's
%           Replace transparency colors with Background
%
% GIF2IMG( GIF_File , BackgroundColor )
%
% Requires UNIX-Command  CONVERT
%
% See also:  READXPM, WRTXPM
%
 
if nargin < 2
   bgc = [ 1  1  1 ];
end

try
   [n,bgc] = colspec(bgc);
catch
   bgc = [];
end

if isempty(bgc)
   error('Invalid BackGroundColor.');
end

f = which(gif);
if ~isempty(f)
    gif = f;
end

if ~( exist(gif,'file') == 2 )
    error(sprintf('File doesn''t  exist: %s',gif));
end

%---------------------------------------------------
% Temporary XPM-Files from CONVERT

[pfd,tmp] = fileparts(tempname);
tmp = fullfile(pfd,sprintf('gif_%s',tmp));
while ( exist(tmp,'dir') == 7 ) | ( exist(tmp,'file') == 2 )
      tmp = sprintf('%s0',tmp);
end
tmp = sprintf('%s.xpm',tmp);

rmf = sprintf('%s*',tmp);       % Files to Remove

%---------------------------------------------------
% Call CONVERT to extract Frames

cmd = sprintf('convert %s %s',gif,tmp);

[s,w] = unix(cmd);

if ~( s == 0 )
    msg = sprintf('Error call UNIX: %s',cmd);
    if ~isempty(w)
        msg = sprintf('%s\n%s',msg,w);
    end
    if ~isempty(dir(rmf))
        unix(sprintf('rm %s',rmf));
    end 
    error(msg);
end

%---------------------------------------------------

frm = sprintf('%s.%s',tmp,'%.0f');

[pfd,gif,ext] = fileparts(gif);

%%% ext = '.png'; %%% !!!

out = fullfile(pfd,'new',sprintf('%s%s%s',gif,'_%3.3d',ext));

z = 0;

msg = '';

while 1

    f  = sprintf(frm,z);
    ok = ( exist(f,'file') == 2 );

    if ~ok & ( z == 0 )
        f  = tmp;
        ok = ( exist(f,'file') == 2 );
    end

    if ~ok
        break
    end

    %----------------------------------------------------------
    % Read XPM

    try
       [msg,c,m] = readxpm(f);
    catch
       msg = lasterr;
    end
 
    if ~isempty(msg)
        msg = sprintf('Error call READXPM( %s )\n%s',f,msg);
        break
    end

    %----------------------------------------------------------
    % NaN --> BackGroundColor

    isn = isnan(c);
    if any(isn(:))
       isn = find(isn);
       ok  = ( m - bgc(ones(1,size(m,1)),:) );
       ok  = all( ok == 0 , 2 );
       if any(ok)
          ok = find(ok);
          c(isn) = ok(1);
       else
          m = cat( 1 , m , bgc );
          c(isn) = size(m,1);
       end
    end

    %----------------------------------------------------------
    % Write XPM

    try
       msg = wrtxpm(c,m,tmp);
    catch
       msg = lasterr;
    end

    if ~isempty(msg)
        msg = sprintf('Error call WRTXPM( %s )\n%s',tmp,msg);
        break
    end

    %----------------------------------------------------------
    % Convert XPM --> GIF

    f = sprintf(out,z+1);  % Start with ONE

    cmd = sprintf('convert %s %s',tmp,f);

    [s,w] = unix(cmd);

    try, delete(tmp), end

    if ~( s == 0 )
        msg = sprintf('Error call UNIX: %s',cmd);
        if ~isempty(w)
            msg = sprintf('%s\n%s',msg,w);
        end
        break
    end

    z = z + 1;

end

if exist(tmp,'file') == 2
   try, delete(tmp), end
end

if ~isempty(dir(rmf))
    unix(sprintf('rm %s',rmf));
end 

if ~isempty(msg)
    error(msg);
end
