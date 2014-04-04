function startup

% STARTUP   local Startup-File, using STARTMENU and STARTCNF
%
% If the EnviromentVariable MATLAB_PATH is defined,
% it's Value will added to Matlab's SearchPath first.
%
% If the EnviromentVariable MATLAB_WORK is defined,
%  the current Working Directory will changed 
%  to this Directory first.
%
% The EnviromentVariable MATLAB_START defines 
%  a Matlab Expression to execute first.
%
% If the first 2 Characters of MATLAB_START are 
%  equal to "@@", Matlab will terminated after
%  the call of MATLAB_START in any case.
%

%------------------------------------------

drawnow

more off

%------------------------------------------
% GraphicDefaults

setdefaults;

%------------------------------------------

mpath = getenv('MATLAB_PATH');
if ~isempty(mpath)
    fprintf(1,'Add MATLAB_PATH: %s\n',mpath);
    try 
       addpath(mpath);
    catch
       fprintf(1,'Error: %s\n',lasterr);
    end
end

mwork = getenv('MATLAB_WORK');
if ~isempty(mwork)
    fprintf(1,'Chdir MATLAB_WORK: %s\n',mwork);
    try 
       cd(mwork)
    catch
       fprintf(1,'Error: %s\n',lasterr);
    end
end

mstart = getenv('MATLAB_START');
if ~isempty(mstart)
    ext = ( size(mstart,2) >= 2 );
    if ext
       ext = strcmp(mstart(1:2),'@@');
       if ext
          mstart = mstart(3:end);
       end
    end
    fprintf(1,'Eval MATLAB_START: %s',mstart);
    if ext
       fprintf(1,', exit');
    end
    fprintf(1,'\n');
    eval(mstart,'fprintf(1,''Error: %s\n'',lasterr)')
    if ext
       exit
    end
end

%------------------------------------------
% StartupMenu, Sessions defined in STARTCNF.M

try
   cnf = startcnf;
catch
   cnf = [];
end

%**********************************************
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
hst = getenv('HOST');

switch hst

  case 'sirius'   % Thinkpad 600, Linux 7.3

     n = prod(size(cnf));
     ok = ones(n,1);
     for ii = 1 : n
         if ~isempty(cnf(ii).WorkDir)
             ok(ii) = ( exist(cnf(ii).WorkDir,'dir') == 7 );
         end
     end
     cnf = cnf(find(ok));


  case { 'kim' 'kimm' }  % Thinkpad A22M, Linux 9.0, Check for Matlab5
                         % Matlab should be started with the ShellScript from the end

     msg = lib_kim(0);   % LIB_KIM see blow !!!

end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%**********************************************

Msg = '';

if ~isempty(cnf)
    try
       [Msg,cnf] = startmenu(cnf);
    catch
       warning(sprintf('Error call STARTMENU.\n%s',lasterr));
       cnf = [];
    end
end

%------------------------------------------

if prod(size(cnf)) == 1
   session = cnf.Name;
else
   session = '';
end

pp = get(0,'defaultfigurepapertype');

fprintf(1,'\n %s STARTUP %s:   ******   %s   ******\n\n', ...
          session,upper(pp(1:2)),pwd);

if ~isempty(Msg)
   fprintf(1,'%s\n\n',Msg);
end


more on


%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function setdefaults


pp = 'a4';
if isunix
   [s,w] = unix('date +%z');   % Deviation to UTC
   if ( s == 0 ) & ~isempty(w)
      try
         if eval(w)/100 <= -5
            pp = 'us';
         end
      end
   end
end


ColorOrder1 = [  
   
         0.00         0.00         1.00
         0.00         1.00         0.00
         1.00         0.00         0.00  ];
        
ColorOrder2 = [
    
         0.00         1.00         1.00
         1.00         0.00         1.00
         0.95         0.95         0.00   ];
     
ColorOrder3 = [ 
    
         0.00         0.75         1.00
         1.00         0.75         0.00
         0.00         1.00         0.75   ];
            
            
ColorOrder = cat( 1 , ColorOrder1 , ColorOrder2 , ...
                      ColorOrder1*0.75 , ColorOrder3 , ColorOrder2*0.75 , ...
                      [ 0.75  0.75  0.75 ] );     


set( 0 , 'defaultfigurecolor'           , 'w' , ...
         'defaultfigurepapertype'       , [ pp 'letter' ] , ...
         'defaultfigureposition'        , [ 420 330 480 360 ] , ...
         'defaultfigureinverthardcopy'  , 'off' , ...
         'defaulttextcolor'             , 'k'        , ...
         'defaulttextinterpreter'       , 'none'     , ...
         'defaultaxescolor'             , 'w'        , ...
         'defaultaxescolororder'        , ColorOrder , ...
         'defaultaxesxcolor'            , 'k' , ...
         'defaultaxesycolor'            , 'k' , ...
         'defaultaxeszcolor'            , 'k' , ...
         'defaultaxesfontsize'          , 10  , ...
         'defaultimagecdatamapping'     , 'scaled' , ...
         'defaultpatchfacecolor'        , 'r' , ...
         'defaultpatchedgecolor'        , 'k' , ...
         'defaultsurfaceedgecolor'      , 'k' , ...
         'defaultlinecolor'             , 'k' , ...
         'hideundocumented'             , 'off' , ...
         'showhiddenhandles'            , 'on'              )
     
    
%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%**********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

% ## /usr/local/bin/matlab5
% #! /bin/bash
% 
% pfd="/opt/matlab53/bin/lnx86/"
% suf=".tmp"
%
% file="libXt.so.6 libX11.so.6"
%
% for ff in $file; do
%     test -e "$pfd$ff"
%     if [ $? -ne 0 ]; then
%        if test -L "$pfd$ff$suf"; then
%           echo "MATLAB5: $pfd$ff$suf -->  $pfd$ff"
%           mv $pfd$ff$suf $pfd$ff
%           wait $!
%        fi
%     fi
% done
%
% /opt/matlab53/bin/matlab -c /opt/matlab61/license/license.dat.solo1
%
