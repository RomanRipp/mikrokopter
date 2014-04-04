function [Msg,cnf] = startmenu(par,group,action,varargin)

% STARTMENU  Opens Startup-Menu and returns StartupConfiguration
%
% [Msg,CNF] = STARTMENU
% [Msg,CNF] = STARTMENU( ConfigurationStructure )
% [Msg,CNF] = STARTMENU( ConfigurationFile      )
%
% Returns the selected ConfigurationStructure CNF with the Fields:
%
%  Name          Name of Session 
%  Label         Label of Session (1 Word)  
%  Startup       Full Name of pecial STARTUP-Script
%  Option        Inputs for STARTUP-Script, CellArray !!!
%  WorkDir       Working Directory
%  EnhcDir       Additional Directories to add to SearchPath, 
%                  CellStringArray or CharacterArray
%  EnhcOpt       Option for additional Directories
%                   'begin' | 'end'
%
% The Directories of "EnhcDir" will added recursivly to Matlabs SearchPath,
%  before the STARTUP-Script will executed.
%
% Give explicit FileNames ("Startup") and DirectoryNames !
%
% The Input-ConfigurationStructure is by default returned from the function STARTCNF.
%
% A ConfigurationFile can be an M-FileFunction, which returnes the ConfigurationStructure
%  or an Script in MatlabSyntax, which defines the ConfigurationStructure in a Variable "c".
%
% After a successfull Startup, the selected ConfigurationStructure (CNF) will written
%  into a Script: '$MatlabHome/startmenu.ini'
% This will used as Item for the StartupMenu of the next MatlabSession.
%
% The Name of the Session will set to the Property "Tag" of the Root.
%   
%  >> get(0,'tag')  % returns the Label or Name of actual session
%
% see also:  STARTCNF,  RECPATH,  ADDPATH,  RMPATH
%


Nin = nargin;

nl = char(10);

Msg  = '';
Msg0 = 'STARTMENU: ';

cnf  = [];

%----------------------------------------------------------
% Basic check

%---------------------------------
if ( Nin == 0 )

  try

    cnf = startcnf;

  catch

    Msg = [ 'Error call STARTCNF.' nl lasterr ];

  end

  action = 'New';

%---------------------------------
elseif  isempty(par)  | isa(par,'struct')  | ...
       ( isa(par,'char') & ( prod(size(par)) == size(par,2) ) )

  cnf = par;
 
  action = 'New';

%---------------------------------
else

  ok = ( isnumeric(par)  &  ( prod(size(par)) == 1 ) );
  if ok
     ok = ishandle(par);
  end

  if ~ok

    Msg = 'First Input must be a STARTCNF-Structure, IniFile or an ObjectHandle.';

  end

end

if ~isempty(Msg)
   Msg = [ Msg0  Msg ];
   return
end

%----------------------------------------------------------
switch upper(action)

%**********************************************************
case 'NEW'

    HomePath = getenv('HOME');
    if isempty(HomePath)
       HomePath = fileparts( which('startup.m') );
    else
       HomePath = fullfile( HomePath , 'matlab' );
    end

    IniFile = fullfile( HomePath , 'startmenu.ini' );

    %-------------------------------------
    % Check Configuration

    [MsgC,cnf] = check_cnf(cnf);

    %  cnf.Name
    %      Label
    %      Startup
    %      Option
    %      WorkDir
    %      EnhcDir
    %      EnhcOpt
    %      StartDir
    %      StartFcn
    %      Text

    if ~isempty(MsgC)

       MsgC = [ Msg0 'Invalid ConfigurationStructure.' nl MsgC ];

    end


    %-------------------------------------
    % Check Configuration from IniFile

    if exist( IniFile , 'file' ) == 2

      [MsgI,ini] = check_cnf(IniFile);

      if ~isempty(MsgI)

       MsgI = [ Msg0 'Invalid IniFile.' nl MsgI ];

      end

    else

      MsgI = '';
      ini  = [];

    end

    %------------------------------------

    nc = prod(size(cnf));
    ni = prod(size(ini));

    Msg = cat( 2 , MsgC , nl(1:(end*(~isempty(MsgC)))) , MsgI );

    if ni+nc == 0
       Msg = [ Msg0 'No Configuration available.' ...
               nl(1:(end*(~isempty(Msg)))) Msg  ];
       return
    end

    %-------------------------------------
    if ( ni == 0 )
    % No valid Configuration in IniFile

       val = 1;

    else
    %-------------------------------------

       if ( nc == 0 ) 
       % No Input/StartCNF-Configuration
        
           cnf = ini;
           val =   1;

       else
       % Check IniConfiguration with Input/StartCNF-Configuration

           ini = ini(1);

           val = 0;

           for ii = 1 : nc
               if isequal( ini(1) , cnf(ii) )
                  val = ii;
                  break
               end
           end

           if val == 0 
           % Old Config not found in New

              ini.Name = cat( 2 , '~' , ini.Name );

              cnf = cat( 1 , ini , cnf(:) );
              val = 1;

           end

       end

    end

    %---------------------------------------
    % Create Figure

    [MsgF,ud,ok,ex] = gui_menu(layout(Msg),'Matlab Startup Menu',cnf,val);
  
    if ~isempty(MsgF)

       Msg = [ Msg0 'Can''t create Startup-Menu-Figure.' ...
               nl  MsgF  nl(1:(end*(~isempty(Msg))))     ...
                         nl(1:(end*(~isempty(Msg))))  Msg    ];

       cnf = cnf(ones(0,1));

       return

    end

    %---------------------------------------
    if ~ok
       cnf = [];
       if ex
          exit
       end
       return
    end

    %------------------------------------------------
    % Execute Startup

    cnf = ud.Config(ud.Value);

    try
      Msg = startmenu_evaluate_startupfcn( cnf );
    catch
      Msg = lasterr;
      cnf = cnf(ones(0,1));
    end

    if ~isempty(Msg)

       Msg = [ Msg0 'Error execute STARTUP.' nl Msg ];

    else

       %---------------------------------------------
       % Write actual Configuration

       wrt_cnf(cnf,IniFile);

       %---------------------------------------------
       % Set Title of xterm, for Linux and XTerm only

       if ~isempty(findstr(upper(computer),'LNX'))

           if strcmp( getenv('TERM') , 'xterm' )

              v = version;
              name = cat(2,'Matlab',v(1));
             

              if ~isempty(cnf.Label)
                 name = cat(2,name,' ',cnf.Label);
              end

              command = cat( 2 , '! echo -ne "\033]0;' , name , '\007"' );

              try, eval(command); end
            
           end

       end

       %---------------------------------------------
       % Set Tag of ROOT

       if isempty(cnf.Label)
          cnf.Label = cnf.Name;
       end

       set( 0 , 'tag' , cnf.Label );


    end

    %--------------------------------------------------


%**********************************************************
case 'LIST'

  [Msg,fig] = recpar(par,'Root');
 
  if ~isempty(Msg) 
     return
  end

  fig = fig(1);

  ud = get(fig,'userdata');

  pud = get( par , 'userdata' );

  ch = getfield( pud.Children.Control , group );
  h = getfield( ch , action );

  val = get( h , 'value' );

  if ud.Value == val
     return
  end

  ud.Value = val;
  
  set( pud.Children.Control.Par.List , 'value'      , 1 , ...
                                       'listboxtop' , 1 , ...
                                       'string'     , ud.Config(val).Text );
              
  set( fig , 'userdata' , ud );

                          
%**********************************************************
case { 'OK'  'QUIT'  'CANCEL' }

  [Msg,fig] = recpar(par,'Root');
 
  if ~isempty(Msg) 
     return
  end

  fig = fig(1);

  pud = get(par,'userdata');

  ch = getfield( pud.Children.Control , group );
  h = getfield( ch , action );

  set( h , 'userdata' , 1 );

  ud = get(fig,'userdata');

  set( fig , 'waitstatus' , ud.ResumeStatus );

end


%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,ud,ok,ex] = gui_menu(gui,titel,cnf,val)


    ud = struct( 'Config' , { cnf } , ...
                 'Value'  , { val }       );

    ok = 1;
    ex = 0;

    Msg = '';

    %***************************************
    % Check for non-exiting X

    if isunix
       if strcmp( get(0,'TerminalProtocol') , 'none' ) | ...
                ( get(0,'ScreenDepth')      ==   0   ) | ...
             all( get(0,'ScreenSize')       <  100   )

          [ud,ok,ex] = asc_menu(titel,ud);

          return

       end
    end
    %***************************************

    %---------------------------------------
    % Create Figure
  
    [Msg,fig] = make_gui(gui,titel);

    if ~isempty(Msg)

       return

    end

    %---------------------------------------
    % Initialize Figure
  
      ud = get(fig,'userdata');
     fud = get( ud.Children.Frame,'userdata');

      hs = fud.Children.Control.Select.List;   % SelectPopup
      hl = fud.Children.Control.Par.List;      % ParameterList

      hc = fud.Children.Control.Control;       % ControlButtons


      set( hs , 'string' , cellstr(str2mat(cnf.Name)) , ...
                'value'  , val                  );

      set( hl , 'string' , cnf(val).Text , ...
                'value'  , 1             , ...
            'listboxtop' , 1                    );

    %--------------------------------------------------
    % Activate

    ResumeStatus = 'inactive';

    ResumeCB     = get( hc.Cancel , 'callback' );
    ResumeStatus = 'inactive';

    ud.ResumeStatus = ResumeStatus;
    ud.Config       = cnf;
    ud.Value        = val;

    set( fig , 'userdata'         ,  ud        , ... 
               'handlevisibility' , 'callback' , ...
               'WindowStyle'      , 'modal'    , ...
               'waitstatus'       , 'waiting'  , ...
               'CloseRequestFcn'  , ResumeCB   , ...
               'visible'          , 'on'        );

    waitfor( fig , 'waitstatus' , ResumeStatus );

    %------------------------------------------------
    % Selection was done

    ud  = get(fig,'userdata');

     ok = get( hc.Ok   , 'userdata' );
     ex = get( hc.Quit , 'userdata' );

    delete(fig);

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ud,ok,ex] = asc_menu(titel,ud,cnt)


ok = 1;
ex = 0;

mx = 3;  % Max Number of (invalid) Trials

spc = char(32*ones(1,3));

if nargin < 3
   fprintf(1,' ****** %s ******\n\n',titel);
   cnt = 0;
end

nc = size(ud.Config,1);

%***************************************************
% Display Sessions

str = sprintf('Select a Session ...%s',spc);
spl = char(' '*ones(size(str)));

fprintf(1,' %s\n %s\n',str,spl);

for ii = 1 : size(ud.Config,1)
    if ii == ud.Value
       fprintf(1,'%s*%4.0f. %s\n',spc,ii,ud.Config(ii).Name);
    else
       fprintf(1,'%s%5.0f. %s\n',spc,ii,ud.Config(ii).Name);
    end
end

fprintf(1,'\n');
fprintf(1,'%s%s<C> Cancel\n',spc,spc);
fprintf(1,'%s%s<E> Exit\n',spc,spc);
fprintf(1,'\n');

%***************************************************
% Get Selection

inp = sprintf('%s Enter your choice, "-N%s" for Info, {%.0f}: ', ...
             char(7),char(186),ud.Value);

sel = input(inp,'s');

fprintf(1,'\n');

%---------------------------------------------------
% Check Selection

sel = rmblank(sel,2);

ok = isempty(sel);

inf = ~ok;

if inf
   inf = ( sel(1) == '-' );
   if inf
      ok = ( size(sel,2) == 1 );
      if ~ok
          sel = sel(2:end);
      end
   end
end

if ~ok
    ex = strcmp(upper(sel(1)),'E') - strcmp(upper(sel(1)),'C');
    ok = ( ex == 0 );
    if ok
       sel = eval(sel,'0');
       ok = ( isnumeric(sel) & ( prod(size(sel)) == 1 ) );
       if ok
          ok = ( ( 1 <= sel ) & ( sel <= nc ) );
          if ok
             ud.Value = sel;
          end
       end
    end
end

if ( ok & ~inf ) |  ~( ex == 0 )
   ex = ( ex == 1 );
   return
end

%***************************************************
% Info  |  Invalid Selection

%---------------------------------------------------
if ~inf % Invalid
%---------------------------------------------------

    cnt = cnt + 1;

    ret = ( cnt >= mx );

    if ret
       str = 'No trials left';
    else
       str = sprintf('%.0f of %.0f trials left',mx-cnt,mx);
    end

    str = sprintf('!!! Invalid Selection !!!   %s ...%s',str,spc);
    spl = char('-'*ones(size(str)));

    fprintf(1,'%s %s\n %s\n %s\n',char(7),spl,str,spl);

    if ret
       return
    end

%---------------------------------------------------
else  % Info
%---------------------------------------------------

   str = sprintf('Info %.0f. Session : %s%s',ud.Value,ud.Config(ud.Value).Name,spc);
   spl = char('-'*ones(size(str)));

   fprintf(1,' %s\n %s\n %s\n\n',spl,str,spl);

   txt      = ud.Config(ud.Value).Text(:,[1 1]);
   txt(:,1) = {spc};
   txt      = permute(txt,[2 1]);

   fprintf(1,'%s%s\n',txt{:});

   str = sprintf('Press any key to continue ...%s',spc);
   spl = char('-'*ones(size(str)));

   fprintf(1,'\n%s %s\n %s',char(7),spl,str);

   pause

   fprintf(1,'\n %s\n',spl);

%---------------------------------------------------
end
%---------------------------------------------------

fprintf(1,'\n');

%***************************************************

[ud,ok,ex] = asc_menu(titel,ud,cnt);


%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  Msg = startmenu_evaluate_startupfcn(cnf);

% START  executes Startup defined in Configuration
%
%  cnf.Name
%      Label
%      Startup
%      Option
%      WorkDir
%      EnhcDir
%      EnhcOpt
%      StartDir
%      StartFcn
%      Text

Msg = '';
nl  = char(10);


%-------------------------------------------------------------
% Add Optional Directories

for pp = cnf.EnhcDir(:)'

    MsgR = recpath(pp{1},['-' cnf.EnhcOpt] );

    if ~isempty(MsgR)
        Msg = [ Msg nl(1:(end*(~isempty(Msg)))) MsgR ];
    end

end

%-------------------------------------------------------------
% Change to StartDir or WorkDir

MsgD = '';
MsgW = '';

if ~isempty(cnf.StartDir)
    MsgD = change_dir(cnf.StartDir,'Startup');
    if ~isempty(MsgD)
      Msg = [ Msg nl(1:(end*(~isempty(Msg)))) MsgD ];
    end
elseif ~isempty(cnf.WorkDir)
    MsgW = change_dir(cnf.WorkDir,'Working');
    if ~isempty(MsgW)
      Msg = [ Msg nl(1:(end*(~isempty(Msg)))) MsgW ];
    end
end

%-------------------------------------------------------------
% Execute StartFcn( Option )

MsgS = '';

if ~isempty( cnf.StartFcn )  &  isempty(MsgD)

   try
     feval( cnf.StartFcn , cnf.Option{:} );
   catch
     MsgS = lasterr;
   end
 
   if ~isempty(MsgS)
      Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
              'Error call ' upper(cnf.StartFcn) '.'  nl MsgS  ];
   end

end

%-------------------------------------------------------------
% Change to WorkDir if not done before

if ~isempty(cnf.StartDir) & ~isempty(cnf.WorkDir)
    MsgW = change_dir(cnf.WorkDir,'Working');
    if ~isempty(MsgW)
      Msg = [ Msg nl(1:(end*(~isempty(Msg)))) MsgW ];
    end
end

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function Msg = change_dir(pfad,lab)

Msg = '';

if nargin < 1
   pfad = '';
end

if nargin < 2
   lab = '';
end

if isempty(pfad)
   return
end
 
try
   cd(pfad)
catch
   Msg = lasterr;
end

if ~isempty(Msg)
    Msg = sprintf('Can''t change to %sDirectory: \n%s',lab,pfad);
end

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function   [Msg,cnf] = check_cnf(c1);

% CHECK_CNF  Checks StartupConfiguration
%
%  cnf.Name
%      Label
%      Startup
%      Option
%      WorkDir
%      EnhcDir
%      EnhcOpt 
%
%      StartDir
%      StartFcn
%      Text
%

Msg = '';
nl  = char(10);

c0 = struct( 'Name'     , { 'No Name' } , ...
             'Label'    , { '' } , ...
             'Startup'  , { '' } , ...
             'Option'   , { {} } , ...
             'WorkDir'  , { '' } , ...
             'EnhcDir'  , { {} } , ...
             'EnhcOpt'  , { 'begin' } , ...
             'StartDir' , { '' } , ...
             'StartFcn' , { '' } , ...
             'Text'     , {{''}}        );


cnf = c0(ones(0,1));

if isempty(c1)
   return
end

%-----------------------------------------------------
% Check for IniFile

if ischar(c1) & ( prod(size(c1)) == size(c1,2) )

   if ~( exist(c1,'file') == 2 )
      Msg = [ 'File ' c1 ' does not exist.' ];
      return
   end

   %-------------------------------------------------
   % Check for Function

   ok = 0;

   p0 = cd;

   [p,f] = fileparts(c1);

   if ~strcmp(f,'startmenu');  % !!!!!!!

     try

       if ~isempty(p)
          cd(p);
       end

       c = feval(f);

     catch

        c = [];

     end

     cd(p0);

     ok = isa(c,'struct');

   end


   if ~ok
   %-------------------------------------------------
   % Check for Script

     fid = fopen(c1,'r');
   
     if fid == -1
        Msg = [ 'Can''t open File ' c1 ];
        return
     end
 
     bb = fread(fid,'char');

     fclose(fid);

     ok = 1;

     eval( char(bb(:)') , 'ok=0;' );

     if ok
        ok = ( exist('c','var') == 1 );
        if ok
           ok = isa(c,'struct');
        end
     end

   end 
   % ~ok == ~Function

   if ~ok
        Msg = [ 'ConfigurationFile must be a Function, returns a Structure,' nl ...
                ' or a Script, defines the Structure "c".' ];
        return
   end

   c1 = c;

end

%-----------------------------------------------------
% Check for Structure

if ~isa(c1,'struct')
   Msg = 'Input must be a ConfigurationFile or a Structure.';
   return
end

%-----------------------------------------------------
% c1 --> cnf

n = prod(size(c1));

cnf = c0(ones(n,1));

f0 = fieldnames(c0);
f1 = fieldnames(c1);

for ff = f1(:)'

    if any(strcmp(ff{1},f0))

       for ii = 1 : n

           cnf = setfield( cnf , { ii } , ff{1} , ...
                           getfield( c1 , { ii } , ff{1} )  );
                     
       end
 
    end

end

%-----------------------------------------------------
% Check Values in cnf

%  cnf.Name
%      Label
%      Startup
%      Option
%      WorkDir
%      EnhcDir
%      EnhcOpt 
%

ok = ones(n,1);

for ii = 1 : n

    msg = '';

    name = '';

    %---------------------------------------------
    % Check Values

    for ff = { 'Name'  'Label'  'Startup' 'Option'  'WorkDir'  'EnhcDir'  'EnhcOpt ' }

        mv = '';

        v = getfield( cnf , {ii} , ff{1} );
 
        %--------------------------------------------------------------
        if isempty(v)

           v = getfield(c0,ff{1});

        %--------------------------------------------------------------
        elseif any( strcmp( ff{1} , { 'Name'  'Label'  'Startup'  'WorkDir'  'EnhcOpt ' } ) )

           if ~ischar(v) & ( prod(size(v)) == size(v,2) )

              mv = [ 'Value for Field "' , ff{1} , '" must be a String.' ];

           else

              v = rmblank(v,2);

              switch ff{1}

                %--------------------------------------------------------
                case 'Name'

                   name = v;

                %--------------------------------------------------------
                case 'Startup'

                   [p,f] = fileparts(v);

                   if isempty(f)
                      mv = 'Invalid Filename for Field "Startup".';
                   else
                      v = cat( 2 , fullfile(p,f) , '.m' );
                   end
                  
                %--------------------------------------------------------
                case 'WorkDir'

                    if ~( exist(v,'dir') == 7 )
                       mv = 'Directory for Field "WorkDir" does not exist.';
                    end

                %--------------------------------------------------------
                case 'EnhcOpt '           

                    if isempty( findstr(v,'end') )
                       v = 'begin';
                    else
                       v = 'end';
                    end

              end
              % switch

           end
           % String ok

        %--------------------------------------------------------------
        elseif strcmp( ff{1} , 'Option' )

           if ~isa(v,'cell')
 
               v = { v };

           end

           v = v(:)';

        %--------------------------------------------------------------
        elseif strcmp( ff{1} , 'EnhcDir' )

           if ischar(v)
              v = cellstr(v);
           end

           if iscellstr(v)

               v = v(:);
              nv = size(v,1);
             vok = zeros(nv,1);

               for jj = 1 : nv

                   vok(jj) = ( exist(v{jj},'dir') == 7 );

               end

               v = v(find(vok));
 
           else
 
              mv = 'Value for Field "EnhcDir" must be a CharacterArray or CellStringArray.';   

           end

        end
        % ~isempty(v)

        if isempty(mv)
           cnf = setfield( cnf , {ii} , ff{1} , v );
        else
            msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , mv );
        end

  
   end
   % ff, Check Values

   %---------------------------------------------

   ok(ii) = isempty(msg);

   if ~ok(ii)

       if ~isempty(name)
           name = sprintf('Invalid Values for %.0f. Session: %s.',ii,name);
       else
           name = sprintf('Invalid Values for %.0f. Session.',ii);
       end

       Msg = cat( 2 , Msg , nl(1:(end*(~isempty(Msg)))) , ...
                  name , nl , msg , nl );

   else
   % Fill additional Fields
  
       %  cnf.Name
       %      Label
       %      Startup
       %      Option
       %      WorkDir 
       %      EnhcDir
       %      EnhcOpt 
       %
       %      StartDir
       %      StartFcn
       %      Text

       v = cnf(ii);

       %----------------------------------------------

       if ~isempty(v.Startup)

         [v.StartDir,v.StartFcn] = fileparts(v.Startup);

       end

       %----------------------------------------------
       % Text

       na = max( 1 , prod(size(v.EnhcDir)) ); 

       nt = 6 + na;

       ia = ( 1 : na ) + 5;
 
       v.Text    = cell(nt,1);
       v.Text(:) = {''};
     
       opt = rmblank( var2mstr( v.Option(:)' ) , 2 );
       opt = opt( 1+strcmp(opt(1),'{') : end-strcmp(opt(end),'}') );

       v.Text{ 1} = [ 'StartupDir: ' v.StartDir ]; 
       v.Text{ 2} = [ 'StartupFcn: ' v.StartFcn ]; 
       v.Text{ 4} = [ 'WorkingDir: ' v.WorkDir ];
       v.Text(ia) ={[ '            '            ]};
       v.Text{ 6} = [ 'EnhanceDir: '            ];
       v.Text{nt} = [ 'EnhanceOpt: '            ]; 

       if ~isempty(v.Option) & ~isempty(v.StartFcn)
          opt = rmblank( var2mstr( v.Option(:)' ) , 2 );
          opt = opt( 1+strcmp(opt(1),'{') : end-strcmp(opt(end),'}') );
          v.Text{2} = [ v.Text{2} '(' opt ')' ]; 
       end

       if ~isempty( v.EnhcDir )
          for jj = 1 : na
              v.Text{ia(jj)} = [ v.Text{ia(jj)}  v.EnhcDir{jj} ];
          end
          v.Text{nt} = [ v.Text{nt}  v.EnhcOpt ];
       end

       %----------------------------------------------

       cnf(ii) = v;
              
   end

end
%  ii = 1 : n

cnf = cnf(find(ok));

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  wrt_cnf(cnf,file);

% WRT_CNF  writes Configuration into File
%
%  cnf.Name
%      Label
%      Startup
%      Option
%      WorkDir
%      EnhcDir
%      EnhcOpt 
%      StartDir
%      StartFcn
%      StartOpt
%      Text


nl = char(10);

fid = fopen( file , 'wt' );

if isequal( fid , -1 )
   return
end

cl = round(clock);

fprintf( fid , '%% Matlab Startup Configuration' );
fprintf( fid , nl );
fprintf( fid , '%% %2.2d.%2.2d.%.0f %2.2d:%2.2d:%2.2d' , cl([3 2 1 4 5 6]) );
fprintf( fid , nl );
fprintf( fid , nl );
fprintf( fid , 'c.Name     = ''%s'';' , cnf.Name );
fprintf( fid , nl );
fprintf( fid , 'c.Label     = ''%s'';' , cnf.Label );
fprintf( fid , nl );
fprintf( fid , 'c.Startup  = ''%s'';' , cnf.Startup );
fprintf( fid , nl );
fprintf( fid , 'c.Option   = %s;' , var2mstr(cnf.Option(:)') );
fprintf( fid , nl );
fprintf( fid , 'c.WorkDir  = ''%s'';' , cnf.WorkDir );
fprintf( fid , nl );
fprintf( fid , 'c.EnhcDir  = %s;' , var2mstr(cnf.EnhcDir(:)') );
fprintf( fid , nl );
fprintf( fid , 'c.EnhcOpt  = ''%s'';' , cnf.EnhcOpt );
fprintf( fid , nl );

fclose(fid);


%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = var2mstr(val)

% VAR2MSTR   Converts Variables to String with Matlab expression
%
% String = VAR2MSTR( V )
%
%  Converts 2-dimensional MatlabVariables 
%   to their String-Expression, using by EVAL 
%             
% Variables could be Numeric, Char, Cell- or StructArrays.
%


 % Numeric Format

   clear eps
   nk    = ceil(abs(log(eps)/log(10))) + 1 ;
   kform = '%%.%.0fg';

   % Usage to convert Number x to char:
   % vk = floor(log(abs(x))/log(10))+1;
   % vk = vk*(vk>0);
   % sprintf(sprintf(kform,nk+vk),x);

 % Seperator for Dimensions
   sep1 = ';';
   sep2 = ',';

 % Make 2-Dimensional
   si = size(val);
  val = reshape(val,si(1),prod(si(2:end)));

   si = size(val);
   si = si * (~isempty(val));


 %**************************************************
 % Char
 if ischar(val) 
   if si(1) <= 1
      str = [ ''''  val '''' ];
   else
      str = '[' ; 
      for ii = 1 : si(1)
        str = [ str '''' val(ii,:) '''' ...
                sep1(1:(end*(ii~=si(1))))  ];
      end
      str = [ str  ']' ];
   end 

 %**************************************************
 % Numeric
 elseif isnumeric(val)

   clear eps
   nk    = floor(abs(log(eps)/log(10)))-1 ;
   kform = ' %%.%.0fg ';
   if prod(si) == 1
     if val == 0
       str = ' 0 ';
     else
        vk = floor(log(abs(val))/log(10))+1;
        vk = vk*(vk>0);
        str = sprintf(sprintf(kform,nk+vk),val);
     end
   else
      str = '[' ; 
      for ii = 1 : si(1)
        for jj = 1 : si(2)
          if val(ii,jj) == 0
            str = [ str  ' 0 ' sep2(1:(end*(jj~=si(2))))  ];
          else
            vk = floor(log(abs(val(ii,jj)))/log(10))+1;
            vk = vk*(vk>0);
            str = [ str sprintf(sprintf(kform,nk+vk),val(ii,jj)) , ...
                      sep2(1:(end*(jj~=si(2))))  ];
          end
        end
        str = [ str  sep1(1:(end*(ii~=si(1))))  ];
      end
      str = [ str  ']' ];
   end 

 %**************************************************
 % CellArray
 elseif iscell(val)

      str = '{' ; 
      for ii = 1 : si(1)
        for jj = 1 : si(2)
          str = [ str  var2mstr(val{ii,jj})           ...
                     sep2(1:(end*(jj~=si(2))))     ];
        end
        str = [ str  sep1(1:(end*(ii~=si(1))))  ];
      end
      str = [ str  '}' ];
 

 %**************************************************
 % StructArray
 elseif isstruct(val)

    fnames = fieldnames(val);
    fnames = fnames(:);
    nf     = size(fnames,1);
   
    fsep = ',';

    str = [ 'struct(' ];
  
    for ff = 1 : size(fnames,1);

      str = [ str  var2mstr(fnames{ff}) fsep '{' ]; 

      for ii = 1 : si(1)
        for jj = 1 : si(2)
          str = [ str  var2mstr(getfield(val,{ii,jj},fnames{ff})) ...
                     sep2(1:(end*(jj~=si(2))))     ];
        end
        str = [ str  sep1(1:(end*(ii~=si(1))))  ];
      end

      str = [ str '}' fsep(1:(end*(ff~=nf)))  ];

    end

    str = [ str ')' ];

 end

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = rmblank(str,dim,cc)

% RMBLANK  Remove Blanks, NewLines at Begin and End of CharacterArrays
%
% String = RMBLANK( CharArray )
%
% CharArray  2-dimensional CharacterArray
%
% further Options:
%
% String = RMBLANK( CharArray , DIM , CHAR )
%
%  
%  DIM  specifies Dimension to work, 
%       default: 2
%
%  CHAR specifies BlankCharacters to remove
%       default:  [ 32  13  10  9 ];  % [ Space CR LF TAB ]
%

  
msg = '';
 nl = char(10);

Nin = nargin;

if Nin < 1
  error('Not enough Input Arguments.')
else
  str0 = str;
  if ischar(str)
    str = double(str);
  end
  ok = isnumeric(str);
  if ok
    ok = all( ( mod(str(:),1) == 0 )  & ...
              ( str(:) >= 0 ) & isfinite(str(:))  );
  end
  if ~ok
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CharArray must be a String or ASCII-Codes.'];
  end
  if size(str,1)*size(str,2) ~= prod(size(str))
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CharArray must be 2-dimensional.'];
  end     
end

if Nin < 2
  dim = 2;
else
  if ~isnumeric(dim)
    msg = [ msg nl(1:(end*(~isempty(msg)))) ...
            'Input DIM must be numeric.' ];
  elseif ~isempty(dim)
    dim = dim(:)';
    if ~all( ( dim == 1 ) |  ( dim == 2 ) )
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
             'Values for Input DIM must Integers larger ZERO.' ];
    end
  end 
end

if Nin < 3
  cc = [ 32  13  10  9 ];  % [ Space CR LF TAB ]
else
  if ischar(cc)
    cc = double(cc);
  end
  ok = isnumeric(cc);
  if ok & ~isempty(cc)
    cc = cc(:)';
    ok = all( ( mod(cc,1) == 0 )  & ...
              ( cc >= 0 ) & isfinite(cc)  );
  end
  if ~ok
      msg = [ msg nl(1:(end*(~isempty(msg)))) ...
              'Input CHAR must be a String or ASCII-Codes.'];
  end
end

if ~isempty(msg)
  error(msg)
end


if isempty(str)
 str = '';
 return
end

if isempty(dim) | isempty(cc)
  str = str0;
  return
end



     jj  = find(str == 0 );
 str(jj) = cc(1);

  blank  = 0*str;

  for ii = cc
    blank = ( blank | ( str == ii ) );
  end

  si = size(str);

  for d = dim

    bad = ( sum(blank,3-d) == si(3-d) );
    jj  = find( bad );
    if ~isempty(jj) 

         p  = [ 3-d  d ];
        str = permute(str,p);

         jj = jj(:)';
         nb = size(jj,2);
        jj1 = find( jj ==   ( 1 : nb ) );       % Blank at Begin
        jj2 = find( jj == ( ( 1 : nb ) + ...    % Blank at End
                            ( si(d) - nb ) ) );
        str(:,jj([jj1 jj2])) = [];

        str = permute(str,p);

    end
    
  end

  str = char(str);



%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [Msg,bb,ref] = char2cell(bb,mnl,ref01);

% CHAR2CELL Converts CharArray to CellStringArray
%
% [Msg,CellString] = CHAR2CELL( String )
%
% Converts the CharacterArray String into a CellStringArry.
%
%---------------------------------------------------------
% CHAR2CELL( String , MarkerNewLine )  
%
% Replace the String's defined by MarkerNewLine with NewLine,
%   default: EMPTY
%
%---------------------------------------------------------
% [Msg,CellString,Reference] = ...
%    CHAR2CELL( String , MarkerNewLine , ReferenceCell)
%  
% Returns References in String, using GET_REF, 
%   ReferenceCell = { S1 S2 C1 C2 }, type: >> help get_ref
%


 Msg = '';
 
 Msg0 = 'CHAR2CELL: ';

 ref = cell(0,2);

 Nin  = nargin;
 Nout = nargout;
 
 nl = char(10);

  if isempty(bb)
     bb = cell(0,1);
     return
  end

  % Marker for NewLine
  if Nin < 2
    mnl = '';
  end

%*****************************************************
% Check Inputs

  if iscellstr(bb)
     bb = cellstr(char(bb));
     bb = strhcat(bb,'',1);
  end
  if isempty(bb)
     bb = cell(0,1);
     return
  end

  ok = ( isnumeric(bb) | ischar(bb) );
  if ok

    if ischar(bb)
       bb = double(bb);
    end

    ok = all( ( mod(bb,1) == 0 )  & ( bb >= 0 ) & isfinite(bb)  );

  end

  if ~ok

      Msg = [ Msg0 ...
              'Input String must be a CharacterArray or ASCII-Codes.'];
 
      bb = cell(0,1);

      return

  end

  %---------------------------------------------------
  % Check MarkerNewLine

  if ~( isempty(mnl) |  ...
        ( ischar(mnl) &  ( prod(size(mnl)) == size(mnl,2) ) ) )

    Msg = [ Msg0  'Input MarkerNewLine must be a String.'];

    return

  end

  %---------------------------------------------------
  % Check Reference

  if ( Nin == 3 )  &  ( Nout == 3 )

     [Msg,ref] = get_ref('',ref01{:});

    if ~isempty(Msg)
       Msg = [ Msg0  'Invalid Input for Reference.' nl Msg ];
       return
    end

  end

%*****************************************************

  if ( size(bb,1) > 1 )  &  ( size(bb,2) > 1 )
     bb = cat( 2 , bb , 10*ones(size(bb,1),1) );
     bb = permute( bb , [ 2 1 ] );
     bb = bb(:);
  end

  if ( size(bb,1) > 1 ) 
     bb = permute( bb , [ 2 1 ] );
  end

  %---------------------------------------------------
  % Check Characters

  ok = all( ( bb ==  9 ) |  ...
            ( bb == 10 ) |  ...
            ( bb == 13 ) |  ...
            (  28 <= bb  &   bb <= 126 ) | ...
            ( 160 <= bb  &   bb <= 255 )        );

  if ~ok
    Msg = [Msg0 'Invalid Characters in String.' ];
    return
  end


%*****************************************************
 

  %---------------------------------------------------
  % Remove CR
  bb( find( bb == 13 ) ) = [];


  bb = char(bb);


  %---------------------------------------------------
  % TAB  --> 8 Blanks 
  bb = strrep( bb , char(9) , char(32*ones(1,8)) ); 

  %---------------------------------------------------
  % mnl --> NewLine   % !!!!!
  if ~isempty(mnl)
    bb = strrep( bb , mnl , char(10) ); 
  end

  %---------------------------------------------------
  % Reference
  if ( Nin == 3 )  &  ( Nout == 3 )

     [MsgR,ref] = get_ref(bb,ref01{:});

  end

  %---------------------------------------------------
  % Form CellString

  % 1. "'"     --> "''"
  bb = strrep( bb , char(39) , char([39  39]) ); 

  % 2. NL --> "';'"
  bb = strrep( bb , char(10) , char([39  59  39]) );
  

  bb = eval([  '{'''  bb  '''}' ]);


%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,c,t] = recpar(h,typ);

% RECPAR  returns ParentHistory of Handle
%
% [Msg,HandleHist,TypeHist] = RECPAR( Handle , StopType )
%
%  recurse trough the Parents unto ( ParentType == StopType )
%
%  StopType starts with UpperCase, returns History excluding
%     Handle with ( ParentType == StopType )
%
%  default: StopType = 'Root'  (recurse unto 'figure')
%
%  HandleHist(end) == Handle
%    TypeHist(end) == HandleType
%

Msg = '';
 c  = zeros(0,1);
 t  =  cell(0,1);

if nargin < 1
   Msg = 'Input Handle is missing.';
   return
end

if nargin < 2
   typ = 'Root';
end

%-----------------------------------------------

if isempty(h)
   return
end

ok = ( isnumeric(h) &  ( prod(size(h)) == 1 ) );
if ok
   ok = ishandle(h);
end

if ~ok
   Msg = 'First Input must be a Single Handle.';
   return
end

if ~( ischar(typ) & ~isempty(typ) & ...
      ( prod(size(typ)) == size(typ,2) ) )
   Msg = 'Type must be a String';
   return
end

%-----------------------------------------------

c = h;
t = { get(h,'type') };

z = 1;

t0 = lower(typ);

while ~( c(1) == 0 )  &  ( ~strcmp(t{1},t0) | ( z == 1 ) )

   z = z + 1;

   c = cat( 1 ,         get(c(1),'parent')  , c );

   t = cat( 1 , { lower(get(c(1),'type')) } , t );

end


if strcmp( t{1} , t0 )

  n = 1 + strcmp( typ(1) , upper(typ(1)) );

  c = c(n:z);
  t = t(n:z);

else

   Msg = [ 'Handle has no Parents with Type '''  typ '''.' ]; 

end

%-----------------------------------------------
   
%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function cnf = layout(Msg)

% LAYOUT GUI-Configuration for using MAKE_GUI

CBFcn = 'startmenu';

wws = 34;   % Select
wwl = 50;   % ListBox

wwc = 12.5;   % Buttons
wwq =  6;     % Quit

cls = [ 1.0   1.0   0.8 ];
cll = [ 1.0   1.0   1.0 ];
clc = [ 0.90  1.0   1.0 ];
clq = [ 1.0   0.9   0.9 ];

%-----------------------------------------------------------------
% SelectPopUp

List = stdgui('Type' , 'popupmenu' , 'Width' , wws , ...
              'Text' , 'Select a Session' , 'Color' , cls , ...
              'CBFcn' , CBFcn, 'Option' , 'Select a Session');

Select = struct( 'List' , { { 1  List } } );

%-----------------------------------------------------------------
% ParameterListbox

List = stdgui('Type' , 'listbox' , 'Width' , wwl+10*i , 'Color' , cll , ...
              'Text' , 'Parameter'  );

 Par = struct( 'List' , { { 1  List } } );


%-----------------------------------------------------------------
% MessageListbox

[m,str] = char2cell(Msg);

List = stdgui('Type' , 'listbox' , 'Width' , wwl+10*i , 'Color' , cll , ...
              'Text' , 'Messages' , 'String' , str  );

Messg = struct( 'List' , { { 1  List } } );


%-----------------------------------------------------------------
% ControlButtons

Ok = stdgui('Type' , 'pushbutton' , 'Width' , wwc , 'Color' , clc , ...
            'UserData' , 0 , 'String' , 'Ok' , 'CBFcn' , CBFcn );

Quit   = stdgui('Type' , 'pushbutton' , 'Width' , wwq , 'Color' , clq , ...
                'UserData' , 0 , 'String' , 'Exit' , 'CBFcn' , CBFcn );

Cancel = stdgui('Type' , 'pushbutton' , 'Width' , wwc , 'Color' , clc , ...
                'UserData' , 0 , 'String' , 'Cancel' , 'CBFcn' , CBFcn );

Control = struct( 'Ok'     , { {  1   Ok     } } , ...
                  'Quit'   , { {  2   Quit   } } , ...
                  'Cancel' , { {  3   Cancel } }        );

%-----------------------------------------------------------------

if isempty(Msg)

  cfg = struct( 'Select'    , { Select  }   , ...
               'Separator1' , { NaN     }   , ... 
                'Par'       , { Par     }   , ...
               'Separator2' , { NaN     }   , ... 
                'Control'   , { Control }         );

else

  cfg = struct( 'Select'    , { Select  }   , ...
               'Separator1' , { NaN     }   , ... 
                'Par'       , { Par     }   , ...
               'Separator2' , { NaN     }   , ... 
                'Messg'     , { Messg   }   , ...
               'Separator3' , { NaN     }   , ... 
                'Control'   , { Control }         );
end

%-----------------------------------------------------------------

cnf = struct( 'Control' , { cfg } );


%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function s = stdgui(varargin);

% STDGUI returns DefaultStructure for single GUI-Item
%
%        Type: ''
%       Width: NaN
%       Color: [1 1 1]
%      String: {''}
%       Value: [1 0 1 0.0100 0.1000]
%       Style: ''
%        Text: ''
%    UserData: []
%       CBFcn: ''
%      Option: []


Nin = nargin;

%       Value  Min  Max  SliderStep
val = [  1     0    1    0.01  0.1  ];

 cl = [ 1  1  1 ]; 


s = struct( 'Type'   , {'text' } , ...
            'Width'  , { NaN } , ...  % ButtoWidth in Character
            'Color'  , { cl  } , ...  % BackGroundColor
            'String' , { {''}} , ...  
            'Value'  , { val } , ...
            'Style'  , { ''  } , ...  % Stlye used for WIdth 
            'Text'   , { ''  } , ...  % TextLabel above
          'UserData' , { []  } , ...
            'CBFcn'  , { ''  } , ...
            'Option' , { []  }       );


Nin = 2*floor(Nin/2);
ind = ( 1 : 2 : Nin );

if ~( ( Nin >= 2 ) & iscellstr( varargin(ind) ) )
   return
end

fs = fieldnames( s );

lv = lower( varargin(ind) );


for ii = 1 : length(fs)

    jj = find( strcmp( lower(fs{ii}) , lv ) );

    if ~isempty(jj)

       jj = max(jj);

       v = varargin{2*jj};

       if strcmp( fs{ii},'Value')
          v0 = v(:)';
          n  = min( length(v) , length(val) );
          v  = val;
          v(1:n) = v0(1:n);
       end

       s = setfield( s , fs{ii} , v );

       switch fs{ii}
        case 'Type'

             switch v
              case 'text'
                   s.Color  =  NaN;
             end

        case 'Width'

             % MultipleLine-Edit
             s.Value(3) = s.Value(3) + 1 * ( imag(v) > 1 ) * ...
                                            strcmp(s.Type,'edit'); 

        case 'String'

         if isnan(real(s.Width))  & ...
            ~any( strcmp( s.Type , {'tag_frame' 'tab_list' 'sel_list'} ) )

           s.Width = size(char(s.String),2) + imag(s.Width);

         end             

       end

    end

end

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [Msg,fig] = make_gui(Config,name)

% MAKE_GUI  Creates a GUI-Figure from Structure
%
% [ Msg , fig ] = MAKE_GUI( Config , FigureName )
%
% requires: TAG_FRAME
%
% see also: STDGUI, TAG_FRAME, TAB_LIST, SEL_LIST, MSG_LIST, MSG_LOGO
%

Msg = '';
fig = '';

nl = char(10);

Msg0 = 'MAKE_GUI: ';


if nargin < 1
   Msg = [ Msg0 'Input Config is missing.' ];
   return
end

if nargin < 2
   name = '';
end


%*******************************************

 form = 'GUI %2.2d/%2.2d/%4.0f %2.0f:%2.2d:%2.2d';
        
   cl = clock;
  ind = [3 2 1 4 5 6];

  tag = sprintf(form,round(cl(ind))) ;

%*******************************************
% PrePositioning of Figure
 
[fs,scr_si] = bestfont;

figpos = [ NaN  NaN  ceil( 2/3 * scr_si([3 4])) ];
figpos(1) = 50;
figpos(2) = scr_si(4)-60-figpos(4);


fig = figure( 'paperunits'   , 'inches'     , ...
              'paperorientation','portrait' , ...
              'units'        , 'pixels'     , ...
              'position'     , figpos       , ...
              'color'        , [1 1 1]      , ...
              'menubar'      , 'none'       , ...
              'toolbar'      , 'none'       , ...
              'numbertitle'  , 'off'        , ...
              'name'         ,  name        , ...
              'colormap'     , [ 1  1  1 ]  , ...
              'createfcn'    , ''           , ...
              'tag'          , tag          , ...
              'resize'       , 'off'        , ...
              'visible'      , 'off'        , ...
              'integerhandle'    , 'on'     , ...
              'handlevisibility' , 'callback'           );

%******************************************************************

 fud = struct( 'Children' , { ...
               struct( 'Menu'     , { [] } , ...
                       'Frame'    , { [] }        ) } );


%******************************************************************

field = fieldnames(Config);

is_menu    = any(strcmp(field,'Menu'));
is_logo    = any(strcmp(field,'Logo'));
is_control = any(strcmp(field,'Control'));


%******************************************************************
% UIMenu's

if is_menu

  try
    [Msg,HM] = mkmenu(fig,'Menu',Config.Menu,fig);
  catch
    Msg = lasterr;
  end

  if ~isempty(Msg)
    Msg = [ Msg0  'Error call MKMENU.' nl Msg ];
    return
  end

  fud.Children.Menu = HM;

end

%******************************************************************

if ~( is_logo | is_control )

  set( fig , 'userdata' , fud );
 
  return

end



%******************************************************************
% Frame arround

try
  [Msg,Frame,hb,pos,bh] = tag_frame(fig,'new','TagNumber' , 0 , ...
                                 'Position'        , [ 0  0 -0 -0 ] , ...
                                 'BackGroundColor' , [ 1  1  1 ]    , ...
                                 'ForeGroundColor' , [ 0  0  0 ]          );
catch
  Msg = lasterr;
end


if ~isempty(Msg)
  Msg = [ Msg0  'Error call TAG_FRAME( New ).' nl Msg ];
  return
end

fud.Children.Frame = Frame;

set( fig , 'userdata' , fud );


%******************************************************************
% FrameChildren

 fch = struct( 'Logo'    , { [] } , ...
               'Control' , { [] }      );

%******************************************************************
% Logo

pos0 = [ 0  -1 ];

if is_logo

  try
    [Msg,hl,hf,ht] = tag_frame(Frame,'add',0,'msg_logo',Config.Logo{:});
  catch
    Msg = lasterr;
  end

  if ~isempty(Msg)
    Msg = [ Msg0  'Error call TAG_FRAME( Add  "msg_logo" ).' nl Msg ];
    return
  end

  pos = get(hl,'position');

  pos0(2) = -1 * ( figpos(4) - pos(2) + 1 );

  fch.Logo = struct( 'Frame'   , { hf } , ...
                     'Text'    , { ht } , ...
                     'Message' , { hl }       );

end

%******************************************************************
% UIControls

if is_control

  cnf = make_cnf(fig);

  try
    [Msg,fch.Control,pos] = mkuic(Frame,0,pos0,Config.Control,cnf);
  catch
    Msg = lasterr;
  end

  if ~isempty(Msg)
    Msg = [ Msg0  'Error call MKUIC.' nl Msg ];
    return
  end

end


%******************************************************************
% Fit FigurePosition

 ud = get(Frame,'userdata');

 pos = pos + 2 * ud.TAG_FRAME.BorderWidth;


 figpos([3 4]) = pos;

 figpos([1 2]) = floor( ( scr_si([3 4]) - figpos([3 4]) ) / 2 );

 figpos(2) = figpos(2) - 20;


 orient = { 'portrait'  'landscape' };
 orient = orient{ 1 + ( figpos(3) >= figpos(4) ) };

 set( fig , 'paperorientation' , orient , ...
            'position'         , figpos       );

 [Msg,pos] = tag_frame(Frame,'resize','absolut');

 if ~isempty(Msg) , return, end

 wygiwys(fig);

%******************************************************************


 ud.Children = fch;


 set( Frame , 'userdata' , ud );



%******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,HM] = mkmenu(par,grp,cfg,fig);


Msg = [];
HM  = []; 

nl = char(10);

if isempty(cfg)
  return
end


%************************************************************


field = fieldnames(cfg);
field = field(:)';

    n = size(field,2);

%------------------------------------------------------------
% HandleStructure

HM      = cell( 2 , n+1 );
HM(1,:) = cat( 2 , field , {'Children'} );
HM(2,:) = { {[]} };

HM      = struct( HM{:} );

HC      = zeros( n , 1 );

%------------------------------------------------------------

HP = epsstr(par);

enable = { 'on'  'off' };

for ii = 1 : n

  cc = getfield( cfg , field{ii} );


  CB    = '';
  usd   = [];

  if ischar( cc{2} )
     if ~isempty(cc{2})
       CB = [ cc{2} '('  HP  ','''  grp  ''',''' field{ii} ''','  ...
              sprintf('%.0f',ii)  ',1);'  ];
     end
  end

  if prod(size(cc)) > 3
     usd = cc{3};
  end

  sets = enable{ 1 + isempty(rmblank(cc{1},2)) };

  HC(ii) = uimenu( ...
       'parent'     , par   , ...
       'callback'   , CB    , ...
       'label'      , cc{1} , ...
       'tag'        , ''    , ...
       'enable'     , sets  , ...
       'checked'    , 'off' , ...
       'Interruptible'   , 'off'    , ...
       'BusyAction'      , 'cancel'  );


  ch = [];

  if isstruct(cc{2})
  % Children

    try
       [Msg,ch] = mkmenu(HC(ii),field{ii},cc{2},fig);
    catch
        Msg = lasterr;
    end

    if ~isempty(Msg) 
       Msg = ['Error call MKMENU( ' field{ii} ' ).' nl Msg ];
       return
    end

  end

  ud = struct( 'Root'     , { get(par,'parent') } , ...
               'Parent'   , { par } , ...
               'Children' , { ch  } , ...
               'UserData' , { usd }       );

  set( HC(ii) , 'userdata' , ud );

  HM = setfield( HM , field{ii} , HC(ii) );

end

HM.Children = HC;



%******************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


function [Msg,HC,pos_out] = mkuic(par,nr,pos0,cfg,cnf);

Msg = [];
HC  = [];

pos_out = [0  0];

if isempty(cfg)
  return
end


%************************************************************
% ParentHandle for CallBack

ud = get( par , 'userdata' );

if nr == 0

   HB = par;

else

   ud = get( par , 'userdata' );

   HB = ud.TAG_FRAME.ButtonHandle(nr);

end


HP = epsstr(HB);



%************************************************************

% Fit for single Characters

is_win = strcmp( upper(computer) , 'PCWIN' );


hd  = cnf.Control.HSpace;   % Horizontal Distance between
vd  = cnf.Control.VSpace;   % Vertical Distance between


fcp = get(par,'ForeGroundColor');  % ForeGroundColor

bcp = get(par,'BackGroundColor');  % TextLabelBackGroundColor

% bcp = 'r';     % Use this to check Extension of Text !!!

ct = cnf.Text;  % Config for TextLabel above

%ct.SingleHeight = 1 * ct.PixelFit(1,2) + ct.PixelFit(2,2);
%ct.FontWeight   = 'bold';

ListBoxExt = '_ListBoxFrame';  % Extension for Frame arround ListBox

%************************************************************


field0 = fieldnames(cfg);
field0 = field0(:)';

    ng = size(field0,2);

%------------------------------------------------------------
% HandleStructure

HC      = cell( 2 , ( ng ) );
HC(1,:) = field0;
HC(2,:) = { {[]} };

HC      = struct( HC{:} );

%************************************************************

pos0   = pos0 + [ hd  -vd ];

wmax_g = 0;   % Max. OffsetWidth


%************************************************************
%------------------------------------------------------------
for iig = 1 : ng
% Group

  cfg1 = getfield( cfg , field0{iig} );

%**********************************************************
if isstruct( cfg1)


  field01 = fieldnames( cfg1 );  % Original

  field1 = field01(:)';          % to expand with fields of TAG_FRAME's

      nc = size(field01,1);

  cc = struct2cell( cfg1 );
  cc = cat(1,cc{:});

  %------------------------------------------------------------
  % Expand HandleStructure
  % 
  % Check for Text only

  text_only = 1;

  for iic = nc : -1 : 1

     field2 = field1(iic);

     typ = lower(cc{iic,2}.Type);

     text_only = ( text_only & strcmp(typ,'text') );

     switch typ 

       case 'tag_frame' 
       % Add TagButton-Fields

         field2 = fieldnames(cc{iic,2}.Option);

         field2 = field2(:)';

       case 'listbox'
       % Add Frame

         field2  = cat( 2 , field1(iic) , ...
                        { cat( 2 , field01{iic} , ListBoxExt ) } ); 

     end


     field1 = cat( 2 , field1(1:iic-1) , field2 , ...
                       field1(iic+1:end) );


  end

  %------------------------------------------------------------
  % Structure of Group

  H1      = cell( 2 , size(field1,2) );
  H1(1,:) = field1;
  H1(2,:) = { {[]} };

  H1      = struct( H1{:} );


  pp = cat(1,cc{:,1});    % Positions
  hp = real(pp);          % HorizNr.

  pos_h = pos0;   % Start for Horizontal

  vmin_h = 0;     % Max. OffsetHight

  nh = max(hp);

  %----------------------------------------------
  % Go Horizontal 
  for iih = 1 : nh

    ih = find( hp == iih );  % Same horizontal Position

    vp = imag( pp(ih) );    % VertNr.

    nv = max(vp);
       
    pos_v = pos_h;  % Start for Vertical

    wmax_v = 0;  % Max. OffsetWidth
 
    %--------------------------------------------
    % Go Vertical
    for iiv = 0 : nv

      iv = find( vp == iiv );

      nn = size(iv,1);

      vmin_n = 0;  % Max. OffsetHight
      
      
      %--------------------------------------------
      % Controls at same Position
  
      for iin = 1 : nn

          iic = ih(iv(iin));  % Index in cc

          ccn = cc{iic,2};

          typ = lower( ccn.Type );

        %--------------------------------------------------------
        % Check Width / Position
 
          if isnan(real(ccn.Width))

             ccn.Width = size( char(ccn.String) , 2 ) + i * imag(ccn.Width);

          end

        %--------------------------------------------------------
        % Check Color

          % BackGroundColor

          if isnan(ccn.Color)
             ccn.Color = bcp;
          end

          % ForeGroundColor

          fc0 = fcp;  

          if isequal( ccn.Color , fc0 )

             fc0 = 1 - fc0;
 
          end

        %--------------------------------------------------------
        % Font's  &  Positioning

          switch typ
           case 'tag_frame'
              c = cnf.Tag;
           case { 'tab_list'  'sel_list'  'listbox' }
              c = cnf.List;
           case 'edit'
              c = cnf.Edit;
           case 'text'
              c = cnf.Text;
           case 'button'
              c = cnf.Button;
              typ = 'pushbutton';
           case 'popupmenu'
              c = cnf.Control;
              c.SliderWidth = cnf.List.SliderWidth;
           otherwise
              c = cnf.Control;
          end

        %--------------------------------------------------------
        % Check for Different Width
     
        c = setstyle(c,ccn.Style,cnf);


        %--------------------------------------------------------

         pos_n = pos_v;

        
        %--------------------------------------------------------
        switch typ

          %*******************************************************************
          case 'tag_frame'

          % Note: TagFrames only 1 in Horizontal !!!!


            % CallBack
            if ~isempty(ccn.CBFcn)

              CB = { ccn.CBFcn HB  field0{iig} };

            else
 
              CB = cell(0,1);

            end


             fieldt = fieldnames(ccn.Option);
                 nt = size(fieldt,1)-1;  
                 nb = prod(size(ccn.String));

               ntag = max(nt,nb);
 
%  'Position'    , [  Left  Bottom  -Right -Top  ] , default: [ 5  5 -5 -5]
%                  [  Left  Bottom   Width  High ] , ... and all Combinations ...
%                  [ -Right -Top     Width  High ]



             % First in Group  or Single in Group ==> LeftOffset == TagFrameOffset
             pos_n(1) = pos_n(1) + ( c.HSpace - hd ) * ( ( iig == 1 ) | ( nc == 1 ) );

             % First Group & Single in Group  ==> TopOffset == TagFrameOffset

             pos_n(2) = pos_n(2) + ( -c.VSpace + vd ) * ...
                                   ( ( iig == 1 )  &  ( nc == 1 ) );

               pos = [ pos_n  -c.HSpace  2*ceil(3*cnf.Text.SingleHeight+2*c.BorderWidth) ];


             VarIn = { 'TagNumber'       , max(nt,nb)    , ...
                       'CBFcn'           , CB            , ...
                       'Position'        , pos           , ...
                       'BorderWidth'     , c.BorderWidth , ...
                       'ForeGroundColor' , fcp           , ...
                       'BackGroundColor' , ccn.Color     , ...
                       'FontName'        , c.FontName    , ...
                       'FontUnits'       , c.FontUnits   , ...
                       'FontSize'        , c.FontSize    , ...
                       'FontWeight'      , c.FontWeight  , ...
                       'Interruptible'   , 'off'         , ...
                       'BusyAction'      , 'cancel'             };
                                            
             [Msg,h,hb,pos,bh] = tag_frame(par,'add',nr,'tag_frame',VarIn{:});

             if ~isempty(Msg) , return, end


             % Set ButtonStrings
             for ib = 1 : nb 
                 set(hb(ib),'string',ccn.String{ib});
             end

             posb = zeros( nt+1 , 2 );  %  [ Width  Hight ]

             hh = cat(1,h,hb);

             % Set ButtonChildren
             for it = 0 : nt 
         
                [Msg,H2,posb(it+1,:)] = mkuic( h , it , [0 0] , ...
                                      getfield(ccn.Option,fieldt{it+1}), cnf );

                 if ~isempty(Msg) , return, end


                 H1 = setfield( H1 , fieldt{it+1} , hh(it+1) );

                 ud = get( hh(it+1) , 'userdata' );

                 ud.Children = H2;


                 CB = '';

                 if ~isempty(ccn.CBFcn)  &  ( it > 0 )
      
                   CB = [ ccn.CBFcn '('  HP ','''  field0{iig}   ''','''  ...
                                               fieldt{it+1} ''',1);' ];

                 end

                 set( hh(it+1) , 'userdata' , ud , ...
                                 'callback' , [ get(hh(it+1),'callback')  CB ] );  

                 
             end  

             %------------------------------------------------------------ 

               ud = get( h , 'userdata' );

             %------------------------------------------------------------ 
             % Set Button to Children

               for it = 1 : nt

                  ud.Children = setfield( ud.Children , fieldt{it+1} , ...
                                          ud.TAG_FRAME.ButtonHandle(it)    );

               end           

             %------------------------------------------------------------ 
             % Adjust Frame in Position in Hight !!!

               pos    = max( posb , [] , 1 ) + 2*ud.TAG_FRAME.BorderWidth;
               pos(2) = pos(2) + bh;

               ud.TAG_FRAME.Position(4) = pos(2);

             %------------------------------------------------------------ 

               set( h , 'userdata' , ud );

             %------------------------------------------------------------ 
             % Adjust OriginalPosition in UserData of Parent-Frame !!!

               pud = get(par,'userdata');
               
               if nr == 0

                  jj = find( pud.TAG_FRAME.FrameChildren.Handle == h );

                  pud.TAG_FRAME.FrameChildren.Position{jj}(4) = pos(2);
 
                  set(par,'userdata',pud);

               else

                   hud = get( pud.TAG_FRAME.HideHandle(nr) , 'userdata' );

                    jj = find( hud.Handle == h );

                   hud.Position{jj}(4) = pos(2);

                   set( pud.TAG_FRAME.HideHandle(nr) , 'userdata' , hud )

               end
 

             %------------------------------------------------------------ 

               % Last Group & Single in Group  ==> Offset == TagFrameOffset
               
               last = ( ( iig == ng )  &  ( nc == 1 ) );

               % wmax_v will added to pos_h, move offs same like pos_n above !!!
               % Last Group & Single in Group  ==> LeftOffset == TagFrameOffset

               c.HSpace = (~last) * hd + ...
                            last  * ( 2*c.HSpace - hd );


               c.VSpace = (~last) * vd + ...
                            last  * c.VSpace ;

     
               pos_n(2) = pos_n(2) - pos(2);


               pos = pos([1 2 1 2]);  % Need pos(3) below !!!


          %*******************************************************************
          otherwise
          % uicontrol


            CB = '';

            if ~isempty(ccn.CBFcn)
      
               CB = [ ccn.CBFcn '('  HP ','''  field0{iig}  ''',''' ...
                                           field01{iic} ''',1);' ];

            end


            %------------------------------------------------------
            % Position

            cw = real(ccn.Width);  % CharacterWidth
            ch = imag(ccn.Width);  % CharacterHight


            ww0 = ceil( cw * c.PixelFit(1,1) + c.PixelFit(2,1) + ...
                        c.LineOffset(1) + c.SliderWidth + 2 * c.BorderWidth );

            hh0 =  ch * ( c.PixelFit(1,2) + c.LineOffset(2) ) + ...
                   c.PixelFit(2,2);

            hh0 = ceil( ( ch >  0 ) * hh0 + ...
                        ( ch == 0 ) * c.SingleHeight + 2 * c.BorderWidth );



            %------------------------------------------------------
            % SeparatorPosition
            if strcmp(typ,'separator');

               typ = 'frame';

               ww0 = ~( cw == 0 ) * ww0  + ...
                      ( cw == 0 ) * cnf.Separator.SingleHeight;

               hh0 = ~( ch == 0 ) * hh0  + ...
                      ( ch == 0 ) * cnf.Separator.SingleHeight;

               if ~isempty(ccn.Style)

                  if strcmp( lower(ccn.Style) , 'list' )

                    ww0 = ww0 + ~( cw == 0 ) * 2 * cnf.Control.BorderWidth;
                    hh0 = hh0 + ~( ch == 0 ) * 2 * cnf.Control.BorderWidth;

                  end

               end


               fc0       = fcp;   % ForeGroundColor
               ccn.Color = bcp;   % BackGroundColor


            %------------------------------------------------------
            elseif ~strcmp(typ,'listbox')

               if ~isempty(ccn.Style)

                  if strcmp( lower(ccn.Style) , 'list' )

                    ww0 = ww0 + 2 * cnf.Control.BorderWidth;

                  end

               end

            end
            %------------------------------------------------------

            % HorizontalAlignment
            ht = 'left';
            if   strcmp(typ,'pushbutton')  | ...
               ( strcmp(typ,'edit') & ( ch == 0 ) )
              ht = 'center';
            end


            %------------------------------------------------------
            % Frame around if ListBox !!!

            bw = cnf.Control.BorderWidth * strcmp( typ , 'listbox' );


            hht  = ceil( ct.SingleHeight  * (~isempty(ccn.Text)) );

            %------------------------------------------------------
            % TextLabel above
            if ~isempty(ccn.Text)

                pos = [ pos_n+[bw 0]  ww0 hht ];

               [Msg,htext] = tag_frame(par,'add',nr,'uicontrol' , ...
                           'position'    ,  pos          , ...
                           'style'       , 'text'        , ...
                           'string'      , ccn.Text      , ...
                           'fontname'    , ct.FontName   , ...
                           'fontunits'   , ct.FontUnits  , ...
                           'fontsize'    , ct.FontSize   , ...
                           'fontweight'  , ct.FontWeight , ... 
                           'cdata'       , []            , ...
                           'foregroundcolor'     , fcp   , ...
                           'backgroundcolor'     , bcp   , ...
                           'horizontalalignment' , 'center'         , ...
                           'callback'            , ''   , ...
                           'tag' , [ field0{iig} '.' field01{iic}  '_Text' ]   );

               if ~isempty(Msg) , return, end

               pos_n(2) = pos_n(2) - pos(4);

            end
 

            %------------------------------------------------------

            if ~ischar(ccn.Option)
               ccn.Option = '';
            end


            %------------------------------------------------------

            pos = [ pos_n  ww0  hh0 ];

            %------------------------------------------------------
            % Move vertical if Text single in Vertical

            if strcmp( typ , 'text' )

              dh = ( 2*cnf.Control.BorderWidth + cnf.Control.SingleHeight - ...
                     c.SingleHeight ) / 2 * strcmp( typ , 'text' ) ;
 
              dh = ( 1 - 1/2*text_only ) * dh;

              pos(2) = pos(2) - dh * ( 1 + ( nv == 0 ) * ( nn == 1 ) );

              pos_n(2) = pos_n(2) + dh;

            end

            %------------------------------------------------------
            % Frame around if ListBox !!!

            if strcmp( typ , 'listbox' )

               %------------------------------------------------------
               % Search for following Buttons, cc(:,1) == NaN

               dh = 0;

               ind = ( iic+1 : nc );
 
                jj = find( cumprod( double( isnan( hp(ind) ) ) ) );
 
               if ~isempty(jj)

                 ind = ind(jj);

                 for jj = ind

                     ccb = cc{jj,2};

                    typ1 = lower(ccb.Type);

                   %--------------------------------------------
                   switch typ1
                      case 'listbox'
                        c1 = cnf.List;
                      case 'edit'
                        c1 = cnf.Edit;
                      case 'text'
                        c1 = cnf.Text;
                      otherwise
                        c1 = cnf.Control;
                   end

                   %--------------------------------------------
                   % Check for use Height
                   if ~isempty(ccb.Style)

                      n = size(ccb.Style,2);

                      if strcmp( ccb.Style(n) , upper(ccb.Style(n)) );

                          ccb.Style = cat( 2 , lower(ccb.Style(1:n-1)) , ccb.Style(n) );

                          c1 = setstyle( c1 , ccb.Style , cnf );

                      end

                   end

                    ch1 = imag(ccb.Width);  % CharacterHight

                    hb = ch1 * ( c1.PixelFit(1,2) + c1.LineOffset(2) ) + c1.PixelFit(2,2);

                    hb = ceil( ( ch1 >  0 ) * hb + ...
                               ( ch1 == 0 ) * c1.SingleHeight + 2 * c1.BorderWidth );

                    hb = hb + ( c1.SliderHeight - hb ) * strcmp(typ1,'slider');

                    dh = dh + hb;

                  end
                  % jj
                end
                % ~isempty(jj)

                posf        = pos;
                posf([3 4]) = posf([3 4]) + [ 0  dh ] + 2*bw;
 
                pos_n       = pos_n + [ 1  -1 ] * 2*bw;

                pos([1 2])  = pos([1 2]) + [ 1  -1 ] * bw;

                [Msg,hf] = tag_frame(par,'add',nr,'uicontrol' , ...
                            'position'    , posf         , ...
                            'style'       , 'frame'      , ...
                            'string'      ,  ''          , ...
                            'cdata'       , []           , ...
                            'userdata'    , []           , ...
                            'tooltipstring'       , ''   , ...
                            'foregroundcolor'     , fc0       , ...
                            'backgroundcolor'     , ccn.Color , ...
                            'callback'            , ''        , ...
                            'Interruptible'       , 'off'     , ...
                            'BusyAction'          , 'cancel'  , ...
                            'tag'  , [ field0{iig} '.' field01{iic}  ListBoxExt ]  );

               if ~isempty(Msg) , return, end

               H1 = setfield( H1 , [ field01{iic}  ListBoxExt ] , hf );
                
            end
            % listbox-frame

            %------------------------------------------------------

            [Msg,h] = tag_frame(par,'add',nr,'uicontrol' , ...
                            'position'    , pos          , ...
                            'style'       , typ          , ...
                            'string'      , ccn.String   , ...
                            'fontname'    , c.FontName   , ...
                            'fontunits'   , c.FontUnits  , ...
                            'fontsize'    , c.FontSize   , ...
                            'fontweight'  , c.FontWeight , ...
                            'cdata'       , []           , ...
                            'value'       , ccn.Value(1) , ...
                            'min'         , ccn.Value(2) , ...
                            'max'         , ccn.Value(3) , ...
                            'sliderstep'  , ccn.Value([4 5])/diff(ccn.Value([2 3])) , ...
                            'userdata'    , ccn.UserData      , ...
                            'tooltipstring'       , ccn.Option , ...
                            'foregroundcolor'     , fc0       , ...
                            'backgroundcolor'     , ccn.Color , ...
                            'horizontalalignment' , ht        , ...
                            'callback'            , CB        , ...
                            'Interruptible'       , 'off'     , ...
                            'BusyAction'          , 'cancel'  , ...
                            'tag'                 , [ field0{iig} '.' field01{iic} ]  );

             if ~isempty(Msg) , return, end


             H1 = setfield( H1 , field01{iic} , h );


             pos_n(2) = pos_n(2) - pos(4);


             %------------------------------------------------------
             % Search for following Buttons, cc(:,1) == NaN

             ind = ( iic+1 : nc );
 
             jj = find( cumprod( double( isnan( hp(ind) ) ) ) );
 
             if ~isempty(jj)
      
               ind = ind(jj);

               for jj = ind

                   ccb = cc{jj,2};

                 typ = lower(ccb.Type);

                 if isnan(ccb.Color)
                    ccb.Color = bcp;
                 end

                 fc1 = fcp;
                 if isequal( fc1 , ccb.Color )
                    fc1 = 1 - fc1;
                 end

                 %-------------------------------------
                 switch typ
                    case 'listbox'
                      c1 = cnf.List;
                    case 'edit'
                      c1 = cnf.Edit;
                    case 'text'
                      c1 = cnf.Text;
                    otherwise
                      c1 = cnf.Control;
                 end

                 %-------------------------------------
                 % Check for use Height
                 if ~isempty(ccb.Style)

                    n = size(ccb.Style,2);

                    if strcmp( ccb.Style(n) , upper(ccb.Style(n)) );

                        ccb.Style = cat( 2 , lower(ccb.Style(1:n-1)) , ccb.Style(n) );

                        c1 = setstyle( c1 , ccb.Style , cnf );

                    end

                 end

                  pos(2) = pos_n(2);


                  cw = real(ccb.Width);  % CharacterWidth
                  ch = imag(ccb.Width);  % CharacterHight

                  pos(4) = ch * ( c1.PixelFit(1,2) + c1.LineOffset(2) ) + c1.PixelFit(2,2);

                  pos(4) = ceil( ( ch >  0 ) * pos(4) + ...
                                 ( ch == 0 ) * c1.SingleHeight + 2 * c1.BorderWidth );

                  pos(4) = pos(4) + ( c1.SliderHeight - pos(4) ) * ...
                                       strcmp(typ,'slider');

                  pos_n(2) = pos_n(2) - pos(4);
  

                 CB = '';

                 if ~isempty(ccb.CBFcn)

                   CB = [ ccb.CBFcn '('  HP ','''  field0{iig} ''',''' ...
                                               field01{jj} ''',1);' ];

                 end


                 if ~ischar(ccb.Option)
                    ccb.Option = '';
                 end


                 [Msg,hb] = tag_frame(par,'add',nr,'uicontrol' , ...
                            'position'    , pos           , ...
                            'style'       , typ           , ...
                            'string'      , ccb.String    , ...
                            'fontname'    , c1.FontName   , ...
                            'fontunits'   , c1.FontUnits  , ...
                            'fontsize'    , c1.FontSize   , ...
                            'fontweight'  , c1.FontWeight , ...
                            'cdata'       , []            , ...
                            'value'       , ccb.Value(1)  , ...
                            'min'         , ccb.Value(2)  , ...
                            'max'         , ccb.Value(3)  , ...
                            'sliderstep'  , ccb.Value([4 5])/diff(ccb.Value([2 3])) , ...
                            'userdata'    , ccb.UserData      , ...
                            'tooltipstring'       , ccb.Option , ...
                            'foregroundcolor'     , fc1       , ...
                            'backgroundcolor'     , ccb.Color , ...
                            'horizontalalignment' , ht        , ...
                            'callback'            , CB        , ...
                            'Interruptible'       , 'off'     , ...
                            'BusyAction'          , 'cancel'  , ...
                            'tag'                 , [ field0{iig} '.' field01{jj} ] );

                 if ~isempty(Msg) , return, end

                 H1 = setfield( H1 , field01{jj} , hb );

               end                 
               % jj

             end
             % ~isempty(jj)

             pos(3) = pos(3) + 2*bw;  % Used for wmax           
      
        end
        % typ


        drawnow

        %*******************************************************************

        wm =  pos(3) + c.HSpace + ( hd - c.HSpace ) * ...
                                  ( hd < c.HSpace ) * ...
                         ~strcmp(typ,'tag_frame') * ( iih == nh );

        vm = pos_n(2) - c.VSpace;
  
        wmax_v = max( wmax_v , wm );
        vmin_n = min( vmin_n , vm );

        %*******************************************************************

      end
      % iin

      pos_v(2) = vmin_n;            % Move vertical before next V-Element
    
    end    
    % iiv, Go Vertical

    pos_h(1) = pos_h(1) + wmax_v;   % Move horizontal before next H-Element

    vmin_h = min( vmin_h , pos_v(2) );
      
  end
  % iih, Go Horizontal 

  wmax_g = max( wmax_g , pos_h(1) );

  pos0(2) = vmin_h;    % Move Vertical before next Group


%**********************************************************
elseif isnumeric(cfg1)

  %---------------------------------------------------------
  % Separator

  ok = ( prod(size(cfg1)) == 1 );
  if ok
     ok = isnan(cfg1);
  end

  if ok

     c = cnf.Separator;

     pos  = [ c.HSpace  pos0(2)  -c.HSpace  c.SingleHeight ];

     [Msg,H1] = tag_frame(par,'add',nr,'uicontrol', ...
               'position'        ,  pos          , ...
               'style'           , 'frame'       , ...
               'backgroundcolor' ,  bcp          , ...
               'foregroundcolor' ,  fcp          , ...
               'tag'             ,  'Separator'        );

      if ~isempty(Msg) , return, end


      pos0(2) = pos0(2) - pos(4) - vd;

  end
        
end
%**********************************************************

  HC = setfield( HC , field0{iig} , H1 );

end
% iig, Group

%------------------------------------------------------------ 
% Return Position

pos_out = [ wmax_g  -pos0(2) ];

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function c = setstyle(c,style,cnf);

% style  use Width
% Style  use Width and Hight
% stylE  use Hight
% STYLE  use all  

if isempty(style)
   return
end

n = size(style,2);

UseWidth = strcmp( style(n) , lower(style(n)) );

UseHight = ( strcmp( style(1) , upper(style(1)) ) | ...
             strcmp( style(n) , upper(style(n)) )    );
        
UseAll   = strcmp( style , upper(style) );
   
try

  cf = fieldnames(cnf);
  jj = find( strcmp( lower(style) , lower(cf) ) );

  if ~isempty(jj)

      c1 = getfield( cnf , cf{jj(1)} );

      if UseAll

         c = c1;

      else

         if UseWidth

             c.SliderWidth   = c1.SliderWidth;
             c.PixelFit(:,1) = c1.PixelFit(:,1);
             c.LineOffset(1) = c1.LineOffset(1);
             c.HSpace        = c1.HSpace;

         end

         if UseHight

             c.SliderHeight  = c1.SliderHeight;
             c.PixelFit(:,2) = c1.PixelFit(:,2);
             c.LineOffset(2) = c1.LineOffset(2);
             c.SingleHeight  = c1.SingleHeight;
             c.VSpace        = c1.VSpace;


         end

      end

  end

end 

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function cnf = make_cnf(fig)

% MAKE_CNF   Font-, PixelConfiguration for UIControls


%************************************************

[fs,scr_si,ppi,is_win,fn] = bestfont;


FontName   = 'helvetica';
FontUnits  = 'points';
FontSize   =  fs;
FontWeight = 'normal';

        fs = fs - is_win;

SliderWidth  = 0;
SliderHeight = 10 + is_win;
BorderWidth  = 3;

 TagFontSize    = fs + 2;
 TagFontWeight  = 'bold';
 TagFrameOffset = 3;

ListFontSize    = FontSize;
ListFontName    = fn;

ListSliderWidth = 15 + ( 3 + 3*(ppi>96) ) * is_win;
ListLineOffset  = 2 - [ 0  2 ] * is_win;

TextFontSize    = fs + 2;
TextFontWeight  = 'bold';
TextBorderWidth = 0;

SeparatorHeight  = 3;
SeparatorOffset  = 5;


z2  = zeros(1,2);
z22 = zeros(2,2);

%************************************************

%-----------------------------------------------
Tag = struct(  ...
               ...
'FontName'     , {    FontName     } , ...
'FontUnits'    , {    FontUnits    } , ...
'FontSize'     , { TagFontSize     } , ...
'FontWeight'   , { TagFontWeight   } , ...
'BorderWidth'  , {  BorderWidth    } , ...
'LineOffset'   , { z2              } , ...
'SliderWidth'  , {  SliderWidth    } , ...
'SliderHeight' , {  SliderHeight   } , ...
'PixelFit'     , { z22             } , ...
'SingleHeight' , { 0               } , ...
'HSpace'       , { TagFrameOffset  } , ...
'VSpace'       , { TagFrameOffset  }        );

%-----------------------------------------------
List = struct( ...
               ...
'FontName'     , { ListFontName    } , ...
'FontUnits'    , {     FontUnits   } , ...
'FontSize'     , { ListFontSize    } , ...
'FontWeight'   , {     FontWeight  } , ...
'BorderWidth'  , {     BorderWidth } , ...
'LineOffset'   , { ListLineOffset  } , ...
'SliderWidth'  , { ListSliderWidth } , ...
'SliderHeight' , { ListSliderWidth } , ...
'PixelFit'     , { z22             } , ...
'SingleHeight' , { 0               } , ...
'HSpace'       , { 0               } , ...
'VSpace'       , { 0               }        );

%-----------------------------------------------
Text = struct( ...
               ...
'FontName'     , {     FontName    } , ...
'FontUnits'    , {     FontUnits   } , ...
'FontSize'     , { TextFontSize    } , ...
'FontWeight'   , { TextFontWeight  } , ...
'BorderWidth'  , { TextBorderWidth } , ...
'LineOffset'   , { z2              } , ...
'SliderWidth'  , {     SliderWidth } , ...
'SliderHeight' , { 0               } , ...
'PixelFit'     , { z22             } , ...
'SingleHeight' , { 0               } , ...
'HSpace'       , { 0               } , ...
'VSpace'       , { 0               }        );

%-----------------------------------------------
Edit = struct( ...
               ...
'FontName'     , { FontName    } , ...
'FontUnits'    , { FontUnits   } , ...
'FontSize'     , { FontSize    } , ...
'FontWeight'   , { FontWeight  } , ...
'BorderWidth'  , { BorderWidth } , ...
'LineOffset'   , { z2          } , ...
'SliderWidth'  , { SliderWidth } , ...
'SliderHeight' , { 0           } , ...
'PixelFit'     , { z22         } , ...
'SingleHeight' , { 0           } , ...
'HSpace'       , { 0           } , ...
'VSpace'       , { 0           }        );

%-----------------------------------------------
Control = struct( ...
                  ...
'FontName'     , { FontName     } , ...
'FontUnits'    , { FontUnits    } , ...
'FontSize'     , { FontSize+2   } , ...
'FontWeight'   , { FontWeight   } , ...
'BorderWidth'  , { BorderWidth  } , ...
'LineOffset'   , { z2           } , ...
'SliderWidth'  , { SliderWidth  } , ...
'SliderHeight' , { SliderHeight } , ...
'PixelFit'     , { z22          } , ...
'SingleHeight' , { 0            } , ...
'HSpace'       , { 0            } , ...
'VSpace'       , { 0            }        );

%-----------------------------------------------
Button = struct( ...
                 ...
'FontName'     , {     FontName    } , ...
'FontUnits'    , {     FontUnits   } , ...
'FontSize'     , { TextFontSize    } , ...
'FontWeight'   , { TextFontWeight  } , ...
'BorderWidth'  , { BorderWidth  } , ...
'LineOffset'   , { z2           } , ...
'SliderWidth'  , { SliderWidth  } , ...
'SliderHeight' , { 0            } , ...
'PixelFit'     , { z22          } , ...
'SingleHeight' , { 0            } , ...
'HSpace'       , { 0            } , ...
'VSpace'       , { 0            }        );

%-----------------------------------------------
Separator = struct( ...
                    ...
'FontName'     , { FontName     } , ...
'FontUnits'    , { FontUnits    } , ...
'FontSize'     , { FontSize     } , ...
'FontWeight'   , { FontWeight   } , ...
'BorderWidth'  , { 0            } , ...
'LineOffset'   , { z2           } , ...
'SliderWidth'  , { SliderWidth  } , ...
'SliderHeight' , { 0            } , ...
'PixelFit'     , { z22          } , ...
'SingleHeight' , { SeparatorHeight } , ...
'HSpace'       , { SeparatorOffset } , ...
'VSpace'       , { 0               }        );


%************************************************

cnf = struct( 'Tag'       , { Tag       } , ...
              'List'      , { List      } , ...
              'Text'      , { Text      } , ...
              'Edit'      , { Edit      } , ...
              'Control'   , { Control   } , ...
              'Button'    , { Button    } , ...
              'Separator' , { Separator }       );
              

%***************************************************
% Right Positioning of UIControls

  % DummyAxes

   axe = axes( 'parent'   , fig          , ...
               'units'    , 'normalized' , ...
               'position' , [ 0 0 1 1 ]  , ...
               'visible'  , 'off'              );
  % DummyText

  ht = text('parent'     , axe        , ...
            'units'      ,'pixels'    , ...
            'position'   , [ 1 1 0 ]  , ...
            'string'     , ''         , ...
            'interpreter', 'none'     , ...
            'visible'    , 'off'               );


%------------------------------------------------

for ff = { 'Tag'  'List'  'Text'  'Edit'  'Control'  'Button' }

    c = getfield( cnf , ff{1} );

    set( ht , 'fontname'   , c.FontName   , ...
              'fontunits'  , c.FontUnits  , ...
              'fontsize'   , c.FontSize   , ...
              'fontweight' , c.FontWeight        );

    ext = zeros(2,4);

    for ii = [ 1  2 ]

      set( ht , 'string' , char( double('H') * ones(ii,ii) ) );

      ext(ii,:) = get( ht , 'extent' );

    end

    m =   ext(2,[3 4]) - ext(1,[3 4]);
    n = ( ext(1,[3 4]) - m );

    c.PixelFit = cat( 1 , m , n );

    cnf = setfield( cnf , ff{1} , c );

end

delete(ht);
delete(axe);

%***************************************************

c = cnf.Control.PixelFit;

hd   = ceil( 1.2 * ( 1 * c(1,1) + c(2,1) ) );  % Horizontal Distance between
vd   = ceil( 2/3 * ( 1 * c(1,2) + c(2,2) ) );  % Vertical Distance between



cnf.List.PixelFit(2,2) = 0;  % List: no vertical Offset


for ff = { 'Separator'  'List'  'Text'  'Edit'  'Control'   'Button' }

    c = getfield( cnf , ff{1} );

    switch ff{1}

     %------------------------------------------
     case 'List'

       c.HSpace = 2*hd;
       c.VSpace =   vd;

       c.SingleHeight = c.PixelFit(1,2) + c.LineOffset(2);

     %------------------------------------------
     case { 'Control'  'Text'  'Edit'  'Button' }

       c.HSpace = hd + hd * strcmp( ff{1} , 'Button' );
       c.VSpace = vd;

       f = 1 + 0.2 * ( ( strcmp( ff{1} , 'Control' ) & ~is_win ) + ...
                       strcmp( ff{1} , 'Button'  )                 ) ;

       c.SingleHeight = ceil( f * ( 1 * c.PixelFit(1,2) + c.PixelFit(2,2) ) );

       c.SingleHeight = c.SingleHeight - c.PixelFit(2,2) * ...
                                ( strcmp( ff{1} , 'Edit' ) & is_win );

     %------------------------------------------
     case 'Tag'

       f = 1 + 0.2;

       c.SingleHeight = ceil( f * ( 1 * c.PixelFit(1,2) + c.PixelFit(2,2) ) );


     %------------------------------------------
     case 'Separator'

       c.VSpace =   vd;


   end

    cnf = setfield( cnf , ff{1} , c );

end

%***************************************************

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [fs,ssi,ppi,is_win,fn] = bestfont

% BESTFONT  Returns an optimal FontSize for this computer
%
% [ FontSize , ScreenSize , ScreenPixelsPerInch ] = BESTFONT
%
% [ ... , IsWin , FixedWidthFontName ] = BESTFONT
%
% returns true for PCWIN-System and the FixedWidthFontname
%

is_win = strcmp( upper(computer) , 'PCWIN' );

uni = get(0,'units');       
      set(0,'units','pixels')
ssi = get(0,'ScreenSize');  
      set(0,'units',uni);
          
ppi = get(0,'ScreenPixelsPerInch');

is_tall = -1 + ( ssi(4) >=  480 ) + ...
               ( ssi(4) >=  600 ) + ...
             1*( ssi(4) >= 1050 );

fs =  8 + 2 * is_tall - 1 * is_win;


if is_win
   fn = 'courier';
else
   fn = { 'arrial'  get(0,'fixedwidthfontname') };
   fn = fn{ 1 + ( ssi(4) >= 1050 ) } ;
end


%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [Msg,pfade,mpath] = recpath(pfad,varargin)

% RECPATH   Add's Directory recursivley to Matlab's-SearchPath
%
% [Msg,Directories,MatlabPath] = RECPATH(Directory,options)
%
% Msg           contains a ErrorMessage, empty if all ok.
% Directories   is the CellstrinArray of the added Directories
%
% valid options:
%
%  '-begin'   appends  the directories
%  '-end'     prepends the directories
%
%  '-remove'  removes  the directories
%
%  '-echo'    display's added Directories in MatlabCommand
%
%



Msg0    = 'RECPATH: ';
Msg     = '';
pfade   = {};
nl      = char(10);
tb      = char(32*ones(1,size(Msg0,2)));  % TabSpace
fsp     = filesep;
psp     = pathsep;
comp    = computer;
mpath   = matlabpath;

Nin = nargin;

if Nin < 1
 pfad = cd;
end

if ~( ischar(pfad)  &  ...
      ( prod(size(pfad)) == size(pfad,2) ) ) 
 Msg = 'Pfad-Input must be a nonempty String';
 return
end

VarArgIn = varargin;
VarArgIn = VarArgIn(:);
if ~iscellstr(VarArgIn)
 Msg = [ 'Option-Inputs must be strings.' ];
 return
end

%-------------------------------------

is_win = strcmp( upper(comp) , 'PCWIN' );

%-------------------------------------
% Check Pfad

if isempty(pfad)
 if is_win
   pfad = [ 'C:'  fsp ];
 else
   pfad = fsp;
 end
elseif is_win
 pfad = strrep(pfad,'/',fsp);  % Take care before accepted UNIX-FileSeperator
end
 


if ~( ( exist(pfad,'dir') == 7 )  |  ...
      ( ~is_win  & strcmp(pfad,fsp) ) )
 Msg = [ 'Invalid Directory: '  pfad ];
 return
end

if strcmp(pfad,'.');
  pfad = cd;
end

if strcmp(pfad,'..');
 pfad = cd;
 if isempty(pfad)
   pfad = fsp;
 else
   jj = findstr(pfad,fsp);
   if ~isempty(jj)
     jj(find(jj==size(pfad,2))) = [];
   end
   if ~isempty(jj)
      pfad = pfad(1:max(jj));
   end
 end
end


if ~isempty(findstr(pfad,'..'))
  orgpfad = cd;
  ok = 1;
  try
     cd(pfad);
  catch
     ok = 0;
  end
  if ~ok
    Msg = ['Invalid Directory:' pfad ];
    return
  end
  pfad = cd;
  cd(orgpfad)
  if isempty(pfad)
    if is_win
      pfad = [ 'C:'  fsp ];
    else
      pfad = fsp;
    end
  end
end

pfad = cat( 2 , pfad , fsp(1:(~strcmp(pfad(end),fsp))) );

remsep = (~strcmp(pfad,fsp)) * (~strcmp(comp([1 2]),'MA'));

  pfad = pfad( 1 : (end-remsep) );


%-------------------------------------
% Check Options

Option = '-end';
  mode = 0;

if ~isempty(VarArgIn)
  jj = find(strcmp(VarArgIn,'-begin') | ...
            strcmp(VarArgIn,'-end')   | ...
            strcmp(VarArgIn,'-remove')      );
  if ~isempty(jj)
   Option = VarArgIn{max(jj)};
  end
  mode = any(strcmp(VarArgIn,'-echo'));
end
 

%-------------------------------------
%  GetPath

if mode

  fprintf([ nl Msg0  'Get DirectoryStructure of:' ...
            nl nl tb  strrep(pfad,'\','\\') nl ]);

end


[Msg,pfade] = getpath0([pfad fsp(1:remsep)],remsep,fsp);


%-------------------------------------
%  SetPath

if mode & ~strcmp(Option,'-remove')

    str = { 'Append'  'Prepend' };

    str = str{ 1 + strcmp(Option,'-begin') };

    fprintf([ nl Msg0  str ' Directories to Matlab''s SearchPath:' ...
              nl nl tb ...
              strrep( strhcat(pfade,'',1,[nl tb]) , '\' , '\\' ) nl ])

end

mpath = setpath(pfade,psp,Option,is_win);

if mode & ~strcmp(Option,'-remove')

   fprintf(nl)

end


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@qq

function  mpath = setpath(pfade,psp,Option,is_win);

% SETPATH  Set MatlabPath
%
%  MatlabPath = SETPATH( pfade , PathSep , Option , IsWin );
%

mpath = matlabpath;  % Current MatlabPath

%-------------------------------------------------
% Append and Prepend PathSep

mpath = cat(  2   , psp(1:(end*(~strcmp(mpath( 1 ),psp)))) , ...
            mpath , psp(1:(end*(~strcmp(mpath(end),psp)))) );

%-------------------------------------------------
% Check for existing Path's in current MatlabPath

ind = zeros(0,2);  % [ Start  Lenght ] of existing Path's

mp = mpath;
pf = pfade;

if is_win
   mp = lower(mp);  % !!!
   pf = lower(pf);  % !!!
end

for p = pf(:)'

    jj = findstr( mp , cat( 2 , psp , lower(p{1}) , psp ) );

    if ~isempty(jj)
       for kk = jj(:)'
           ind = cat( 1 , ind , [ kk  size(p{1},2)+1 ] );
       end
    end

end

%-------------------------------------------------
% Remove Path's from current MatlabPath

if ~isempty(ind)

   ii = grp2ind(ind(:,1),ind(:,2));

   mpath(ii) = [];

end

%-------------------------------------------------
% Append or Prepend Pfade

pfade = strhcat( pfade , psp , size(pfade,1)+1 );

switch Option

  case '-begin'

     mpath = cat( 2 , psp , pfade , mpath );
    
  case '-end'

     mpath = cat( 2 , mpath , pfade , psp );

end

%-------------------------------------------------
% Remove PathSep from begin and End

mpath = mpath( 2 : ( end-1) );

%-------------------------------------------------
% Set new MatlabPath

matlabpath( mpath );


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@qq

function [ Msg , pfade ] = getpath0(pfad,remsep,fsp);

% function [ Msg , pfade ] = getpath0(pfad,RemoveSeparator,FileSeperator);


Msg     = '';
pfade   = {};
nl      = char(10);


%-------------------------------------
%  Recurse

d = dir(pfad);

if isempty(pfad)
  return
end

is_dir = cat(1,d.isdir);

if ~any(is_dir)
  return
end

is_dir = find(is_dir);

 pfade    = cell(size(is_dir(:),1),1);
 pfade(:) = { {} };

for ii = 1 : size(pfade,1)

  jj = is_dir(ii);

  if ~any(strcmp(d(jj).name,{ '.'  '..'  'private' }))  & ...
     ~strcmp(d(jj).name(1),'@')

    p = cat( 2 , pfad , d(jj).name , fsp(1:(~strcmp(d(jj).name(end),fsp))) );

    [msg,pfade{ii}] = getpath0( p , remsep , fsp );

    Msg = [ Msg  nl(1:(end*(~isempty(Msg))))  msg ];

  end

end


pfade = cat( 1 , {pfad(1:(end-remsep))} , cat(1,pfade{:}) );


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ii = grp2ind(i0,l);

% GRP2IND  Built IndexVector from StartIndex and Length
%
% Index = GRP2IND( StartIndex , GroupLength )
%

if isempty(i0);
   ii = [];
   return
end

si = size(i0);

if ( sum( si > 1 ) > 1 )
   error('StartIndex must be a Vector.');
end

i0 = i0(:);
l  =  l(:);

if ~isequal(size(i0,1),size(l,1))
   error('Size of StartIndex and GroupLenght must be the same.');
end

n = size(l,1);

ii = ones(sum(l),1);
jj = cumsum( cat(1,1,l) , 1 );

ii(jj(1:n)) = i0;

if n > 1
   ii(jj(2:n)) = ii(jj(2:n))-(i0(1:n-1,1)+l(1:n-1)-1);
end

ii = cumsum(ii,1);

if n == 1
   return
end

jj = find( si > 1 );

if jj == 1
   return
end
  
perm = cat( 2 , (1:jj-1)+1 , 1 , (jj+1:size(si,2)) );

ii = permute(ii,perm);

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

str(    size(str,1),2) = { '' };


str = permute( str , [ 2  1 ] );

str = cat(2,str{:});


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function wygiwys(fig)

% WYGIWYS  WhatYouGetIsWhatYouSee
%
% Switch the PaperPosition to the same View like on the Screen
% The Figure will positioned in the Center of the Paper
%
% WYGIWYS( FigureHandle )
%
%


if nargin < 1
  fig = get(0,'currentfigure');
end

if isempty(fig)
 return
end

ok = ( isnumeric(fig)  &  ( prod(size(fig)) == 1 ) );
if ok
   ok = ishandle(fig);
   if ok
      ok = strcmp( get(fig,'type') , 'figure' );
   end
end

if ~ok
   error('Input must be a FigureHandle.');
end


figuni = get(fig,'units');
papuni = get(fig,'paperunits');

set(fig,     'units' , 'pixels' , ...
        'paperunits' , 'inches'       );

figpos = get(fig,'position');
pappos = get(fig,'paperposition');
pap_si = get(fig,'papersize');

ppi    = get(0,'screenpixelsperinch');

pappos([3 4]) = figpos([3 4]) / ppi;
pappos([1 2]) = (pap_si-pappos([3 4])) / 2;

set(fig,'paperposition',pappos);

set(fig,     'units' , figuni   , ...
        'paperunits' , papuni         );


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [ h , form ] = epsstr( h );

% EPSSTR  Transform Number exact into String,
%          
% using Matlab's floating point relative accuracy.
%
%  Form = EPSSTR;   
%    returns Format for using with SPRINTF
%
%  [ String , Form ] = EPSSTR( Number ); 
%    returns exact String for Number
%
%  Form = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 )
%


form = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 );

if nargin < 1

  h = form;

  return

end


if ~isnumeric(h)  |  ( prod(size(h)) > 1 )
 error('Handle must be a single Numeric.');
end


  h = sprintf(form,h);

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


function  [Msg,varargout] = tag_frame(h,action,varargin);

% TAG_FRAME   Creates a TagFrame
%
% [ Msg, ... ] = TAG_FRAME( Handle , Action , ... );
%
%  Note:  The first Output contains ErrorMessages for invalid Inputs.
%         If all is ok, Msg is empty.
%
%-------------------------------------------------------------------------
% 1. Create a new TagFrame in a Figure, Action 'NEW'
%
%   
% [ Msg, FrameHandle, ButtonHandle, PixelPosition, ButtonHight ] = ...
%    TAG_FRAME( FigureHandle, 'New', Property1, PropertyValue1, ... );
%
% Creates a new TagFrame in the Figure, defined by 
%   FigureHandle, and returns the Handle of the FrameBox, TabButtons, and 
%   the Position in Pixels. The Units for the Position of the ListBox 
%   in the Figure are set to 'pixels'.
% 
% Additional Inputs are Properties and their Values:
%
%  'TagNumber'   , Number of Tab's                 , default: 3
%  'Position'    , [  Left  Bottom  -Right -Top  ] , default: [ 3  3 -3 -3]
%                  [  Left  Bottom   Width  High ] , ... and all Combinations ...
%                  [ -Right -Top     Width  High ]
%  'BorderWidth' , PixelWidth of UIControl-Border  , default: 3
%  'HighPadding' , Padding of TabButton in Hight   , default: 1.5
%  'CBFcn'       , CallBackFcn, called if a Tag is activated:
%                  feval( CBFcn{:} , 'Activate' , 1 , Number )
%
% More Inputs are FontProperties and their Values, by default 
%  the DefaultUIControlFontProperties of the Root will used.
% 
%
% Please:  don't change the Property 'Tag' of the TagFrameHandle,
%           this Property will used to identify the Handle as valid 
%           TagFrame.
%
%
% The UserData of the TagFrameHandle contains all Properties in a StructArray.
%
%-------------------------------------------------------------------------
% 2. Add an Object, Action 'ADD'
%
%  [ Msg , Handles ] = TAG_FRAME( FrameHandle , 'Add' , nr , type , varargin )
%  
%   nr = [ 0 .. ud.TagNumber ]  |  FrameHandle  |  ButtonHandle
%
% type = 'msg_logo' |  'tag_frame'  |  'tab_list'  |  'sel_list'  |  'uicontrol'
%
% Following OutPuts: 
%
%   'msg_logo' :  MessageListHandle, LogoFrameHandle, LogoTextHandle
%   'tag_frame':  TagFrameHandle,  TagButtonHandle, PixelPosition, ButtonHight
%   'sel_list' :  SelFrameHandle,  TextHandle, PixelPosition
%
% UIControls can be positioned equal to the Positioning of the TagFrame, see under NEW
%
%-------------------------------------------------------------------------
% 3. Activate a Tag, Action 'ACTIVATE'
%
%  [ Msg , ActiveNumber ] = TAG_FRAME( FrameHandle , 'Activate' , nr )
%
%  nr = [ 0 .. ud.TagNumber ]  |  FrameHandle  |  ButtonHandle
%
%  nr = NaN sets TagFrame visible off
%
%-------------------------------------------------------------------------
% 4. Set Enability of Tag's, Action 'ENABLE'
%
%  Msg  = TAG_FRAME( FrameHandle , 'Enable' , sets , nr )
%
%  sets = 'off'  |  'on'
%
%   nr = [ 0 .. ud.TagNumber ]  |  FrameHandle  |  ButtonHandle
%  
%-------------------------------------------------------------------------
% 5. Set Visibility of TagFrame, Action 'VISIBLE'
%
%  Msg  = TAG_FRAME( FrameHandle , 'Visible' , sets )
%
%  sets = 'off'  |  'on'
%
%
%-----------------------------------------------------------------------------
% 6. Resize the TagFrame after the Figure was resized, Action 'RESIZE' 
%
%  [ Msg, PixelPosition, ButtonHight ] = TAG_FRAME( FrameHandle, 'Resize', Mode )
%
%  Resize the TagFrame in the Figure,  depending on Mode
%   call this from the ResizeFcnProperty of the Figure.  
%
%  Mode:  'normalized'  ( short:  'n' )
%
%   The TagFramePosition will set normalized to the Figure, with the normalized
%   Position when the TagFrame was created. 
%
%  Mode:  'absolut'     ( short:  'a' )
%
%   The TagFramePosition will set absolut, 
%    defined by the original Position.
%
%-----------------------------------------------------------------------------
% 7. Delete the TagFrame, Action 'DELETE'
%
%  MSG = TAG_FRAME( FrameHandle , 'Delete' )
%
%-----------------------------------------------------------------------------
%
% see also: MAKE_GUI, TAB_LIST, SEL_LIST, MSG_LIST, MSG_LOGO
%
%

Nout = nargout - 1;

Nout = Nout * ( Nout > 0 );

varargout = cell(1,Nout);



Msg = '';

nl = char(10);


Msg0 = 'TAG_FRAME: ';

Nin = nargin;

if Nin < 2
  Msg = [ Msg0  'Inputs H and Action are undefined.' ];
  return
end


ok = ( isnumeric(h) & ( prod(size(h)) == 1 ) );
if ok
 ok = ishandle(h);
end

if ~ok
 Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
         ' First Input must be a Handle.' ];
end


if ~ischar(action)
 Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
          'Action must be a String,'  ];
end

if ~isempty(Msg)
  Msg = [ Msg0   Msg ];
  return
end


action = upper(action);

Msg0 = [ 'TAG_FRAME( ' action ' ): ' ];

%------------------------------------------------------------
% Check Handle

if strcmp(action,'NEW')

 if ~strcmp(get(h,'type'),'figure')
   Msg = [ Msg0  'Handle must be a Figure.' ];
   return
 end

 fig = h;
   h = [];

%------------------------------------------------------------
else

 ok = strcmp(get(h,'type'),'uicontrol');
 if ok
   ok = ( strcmp(get(h,'style'),'frame')      &  ...
          strcmp(get(h,'tag'),'TAG_FRAME')          );
 end

 if ~ok
   Msg = [ Msg0  'Handle must be a ' ...
           'valid TagFrame-UIControl.' ];
   return
 end

 fig = get(h,'parent');

 ud0 = get(h,'userdata');

 ud  = ud0.TAG_FRAME;
  
end
%------------------------------------------------------------


  VarIn = varargin;
  VarIn = VarIn(:);
 

switch action

%*********************************************************
case 'NEW'

  
  if ( mod(size(VarIn,1),2) ~= 0 )
    Msg = [ 'Additional Inputs must contain UIControl-Property-Value-Pairs.' ];
  else 
    VarIn = reshape(VarIn,2,size(VarIn,1)/2)';
    if ~iscellstr(VarIn(:,1))
     Msg = [ 'Additional Inputs must contain UIControl-Property-Value-Pairs.' nl ...
             'Properties must be Strings.' ];
    end
  end

  if ~isempty(Msg)
     Msg = [ Msg0  Msg ];
     return
  end


  % FigurePosition in Pixels

  figuni = get(fig,'units');
           set(fig,'units','pixels');
  figpos = get(fig,'position');
           set(fig,'units',figuni);


  % DefaultUserData

  e1 = zeros(0,1);

  ch = struct( 'Handle'   , { e1        }  , ...
               'Position' , { cell(0,1) }        );

  ud = struct( 'TagNumber'     , {  3  }         , ... % Nr of Tab's
               'Position'      , {[ 3  3 -3 -3 ]} , ... % Position 
               'BorderWidth'   , {  3 }          , ... % BorderPixelWidth of UIControl
               'HighPadding'   , { 1.5 }         , ... % Position per ListboxLine
               'CBFcn'         , { {} }          , ... % CallBackFcn
               'PixelHight'    , {  0 }          , ... % PixelHight of Font
               'FrameHandle'   , { NaN }         , ... % Handle of Frame arround
               'ButtonHandle'  , { e1  }         , ... % Handle of Buttons
               'HideHandle'    , { e1  }         , ... % Handle of TextObject, using for width
               'ActiveNr'      , {  0  }         , ... % Nr of Active Button
               'FrameChildren' , { ch  }         , ... % Children of Frame
               'Visibility'    , { 'on' }        , ... % Visibility of TagFrame
               'PixelPosition' , { zeros(1,4) }  , ... % Actual PixelPosition [ Left Bottom Width High ]
               'NormPosition'  , { zeros(1,4) }  , ... % Default normalized Position
               'FigPos0'       , { figpos     }  );   % Default FigurePosition


  %-------------------------------------------------------------------
  % Look, if the first 5 Properties are given in Inputs

  fields = fieldnames(ud);
  fields = fields([1 2 3 4 5]);

  nv = size(VarIn,1);

  is_ud = zeros(nv,1);

  for ii = 1 : nv

     jj = find( strcmp( lower(VarIn{ii,1}) , lower(fields) ) );

     if ~isempty(jj)

       val = VarIn{ii,2};
       msg = '';

       if strcmp( fields{jj} , 'CBFcn' )

           if ischar(val)
              val = cellstr(val);
           end

           ok = iscell(val);
           if ok & ~isempty(val)
              ok = ( ischar(val{1}) & ~isempty(val{1}) & ...
                    ( prod(size(val{1})) == size(val{1},2) )  );
           end
 
           if ~ok
               msg = [  msg  nl(1:(end*(~isempty(msg)))) ...
                        'Value for '  fields{jj} ' must be a CharacterArray' ...
                        ' or CellArray with a String in the 1. Element.' ];

           end

       elseif ~isnumeric(val) | isempty(val)

            msg = [ 'Value for '  fields{jj} ' must be numeric and not empty.' ];

       elseif strcmp(fields{jj},'Position')

            ok =  all(isfinite(val));
            if ~ok
                msg = [ 'Value for '  fields{jj} ' must contain finite numerics.' ];
            end
            si0 = size(getfield(ud,fields{jj}));
            if ~isequal(size(val),si0)
               str = sprintf(' %.0f by',si0);
               str = [ '[' str(1:end-2) ']' ];
               msg = [  msg  nl(1:(end*(~isempty(msg)))) ...
                        'Value for '  fields{jj} ' must be of Size ' str '.' ];
            elseif ok
               if any( ( val([1 2]) < 0 )  &  ( val([3 4]) < 0 ) )

                  msg = [ 'Invalid Values for '  fields{jj} ...
                          ', type "help tag_frame" for more Informations.' ];

               end
            end

       else
            val_min = 1 - strcmp(fields{jj},'TagNumber');

            ok1 =  ( all( val >= val_min )  &  all(isfinite(val)) );
            ok2 =  ( all( mod(val,1) == 0 )  |  strcmp(fields{jj},'HighPadding') );
            if ~( ok1 & ok2 )
                msg = [ 'Value for '  fields{jj} ' must contain Integers >= ' ...
                         sprintf('%.0f',val_min) '.' ];
            end
            si0 = size(getfield(ud,fields{jj}));
            if ~isequal(size(val),si0)
               str = sprintf(' %.0f by',si0);
               str = [ '[' str(1:end-2) ']' ];
               msg = [  msg  nl(1:(end*(~isempty(msg)))) ...
                        'Value for '  fields{jj} ' must be of Size ' str '.' ];
            end
       end

        if isempty(msg)
          is_ud(ii) = 1;
             ud     = setfield(ud,fields{jj},val);
        else
             Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) msg ];
        end

     end
     % ~isempty(jj)
  end
  % ii

  if ~isempty(Msg)
    Msg = [ Msg0  Msg ];
    return
  end

  
  VarIn(find(is_ud),:) = [];

  %-------------------------------------------------------------------
  % Check UIControlProperties in the other Inputs

  VarIn(:,3) = { [] };  % DefaultValues


  % Set DefaultProperties 

  for ii = 1 : size(VarIn,1)

    msg = '';

    try

      VarIn{ii,3} = get(0,['DefaultUIControl' VarIn{ii,1}]);

    catch

       msg = lasterr;

    end

    if isempty(msg)
       set(0,['DefaultUIControl' VarIn{ii,1}],VarIn{ii,2});
    else
       Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
               'Invalid UIControlProperty: ' VarIn{ii,1} ];
       break
    end

  end

  if ~isempty(Msg)
     % Set DefaultProperties back
    for jj = 1 : ii-1
      set(0,['DefaultUIControl' VarIn{jj,1}],VarIn{jj,3});
    end

    Msg = [ Msg0  Msg ];
    return

  end

  %--------------------------------------------------------
  % Define DummyAxes, containing TextObject to get PixelHigh

  axe = axes('parent'  , fig , ...
                      'units'   , 'pixels' , ... 
                      'position',[ 0 0 1 1 ] , ...
                      'xlim'    , [ 0  1 ] , ...
                      'ylim'    , [ 0  1 ] , ...
                      'xtick'   , [] , ...
                      'ytick'   , [] , ...
                      'visible' , 'off' , ...
                      'nextplot' , 'add' , ... 
                      'handlevisibility' , 'callback'   );

  ht = text('parent' , axe , ...
          'units' , 'pixels' , ...
          'position',[ 1 1 0 ] , ...
          'string'  , 'H' , ... 
          'interpreter' , 'none' , ...
          'fontunits'  , get(0,'DefaultUIControlFontUnits') , ... 
          'fontsize'   , get(0,'DefaultUIControlFontSize')  , ... 
          'fontname'   , get(0,'DefaultUIControlFontName')  , ... 
          'fontangle'  , get(0,'DefaultUIControlFontAngle') , ... 
          'fontweight' , get(0,'DefaultUIControlFontWeight') , ...
          'visible'    , 'off'       );

   % Determine PixelHight

    e = get( ht , 'extent' );

   delete(ht);
   delete(axe);

   ud.PixelHight = e(4);

   %--------------------------------------------------------
 

   ud.ButtonHandle = NaN*ones(ud.TagNumber,1);
   ud.HideHandle   = NaN*ones(ud.TagNumber,1);


   for ii = 1 : ud.TagNumber

     ud.ButtonHandle(ii) = uicontrol('parent'   , fig     , ...
                                     'visible'  , 'off'   , ...
                                     'selected' , 'off'   , ...
                                     'units'    , 'pixels', ...
                                     'position' ,  [ 0  0  1  1 ] , ...
                                     'style'    , 'pushbutton' , ...
                                     'string'   , ''           , ...
                                     'enable'   , 'on'         , ...
                                     'visible'  , 'on'         , ...
                                     'tag'      , 'TAG_BUTTON'  );
 
   end
   
   ud.FrameHandle = uicontrol('parent'   , fig     , ...
                              'visible'  , 'off'   , ...
                              'selected' , 'off'   , ...
                              'units'    , 'pixels', ...
                              'position' ,  [ 0  0  1  1 ] , ...
                              'style'    , 'frame' , ...
                              'enable'   , 'on'         , ...
                              'visible'  , 'on'         , ...
                              'tag'      , 'TAG_FRAME'  );

                      

   form = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 );

     HF = sprintf(form,ud.FrameHandle);


   DelCB = [ 'tag_frame('  HF  ',''Delete'',-1);' ];

   set(ud.FrameHandle,'DeleteFcn',DelCB);


   % ButtonUserData

   udb = struct( 'Root'     , { fig           }  , ...
                 'Parent'   , { ud.FrameHandle } , ...
                 'Children' , { [] }                     );

   for ii = 1 : ud.TagNumber

     CB = [ 'tag_frame('  HF  ',''Activate'','  sprintf('%.0f',ii)  ',1);'  ];

     set( ud.ButtonHandle(ii) , 'callback' , CB , ...
                                'userdata' , udb          );
    
     ud.HideHandle(ii) = uicontrol(  'parent'   , fig     , ...
                                     'visible'  , 'off'   , ...
                                     'selected' , 'off'   , ...
                                     'units'    , 'pixels', ...
                                     'position' ,  [ 0  0  1  1 ] , ...
                                     'style'    , 'text'     , ...
                                     'string'   , { '' }     , ... 
                                     'enable'   , 'on'         , ...
                                     'userdata' ,   ch       , ...
                                     'tag'      , 'TAG_HIDE'  );

    
   end


   % Set DefaultProperties back
   for ii = 1 : size(VarIn,1)
     set(0,['DefaultUIControl' VarIn{ii,1}],VarIn{ii,3});
   end


   ud0 = struct( 'TAG_FRAME' , { ud  } , ...
                 'Root'      , {  0  } , ...
                 'Parent'    , { fig } , ...
                 'Children'  , { []  }         );

   set( ud.FrameHandle , 'userdata' , ud0 );

   %-------------------------------------------
   % Set correct Size of UIControls

   try
     [Msg,pos,ButtonHigh] = tag_frame(ud.FrameHandle,'Resize','absolut');
   catch
      Msg = lasterr;
   end

   if ~isempty(Msg)

     Msg = [ Msg0 'Error call TAG_FRAME( Resize ).' nl Msg ];

     set(ud.FrameHandle,'DeleteFcn','');

     delete(ud.FrameHandle);
     delete(ud.ButtonHandle);
     delete(ud.HideHandle);

     return

   end
 
   ud0.TAG_FRAME.PixelPosition = pos;
   ud0.TAG_FRAME.NormPosition  = pos ./ figpos([ 3  4  3  4 ]);

   set( ud.FrameHandle , 'userdata' , ud0 );

   %-------------------------------------------
   
   try
      Msg = tag_frame(ud.FrameHandle,'Activate');
   catch
      Msg = lasterr;
   end

   if ~isempty(Msg)
     Msg = [ Msg0 'Error call TAG_FRAME( Activate ).' nl Msg ];
     delete(ud.FrameHandle);
     delete(ud.ButtonHandle);
     delete(ud.HideHandle);
     return
   end

   %-------------------------------------------

   out = { ud.FrameHandle  ud.ButtonHandle  pos  ButtonHigh }; 

   n = min(Nout,size(out,2));

   varargout(1:n) = out(1:n);


%*********************************************************
case 'RESIZE'

  if isempty(VarIn)
    mode = 'a';   % Absolut
  else
    mode = VarIn{1};

    ok = ( ischar(mode)  &  ~isempty(mode) );
    if ok
        mode = lower(mode(1,1));
          ok = any(strcmp(mode,{'a'  'n'}));
    end
    if ~ok
      Msg = [  Msg0 'Input must be a String: ''normalized'' or ''absolut''.' ];
      return
    end
  end


  % FigurePosition in Pixels

  figuni = get(fig,'units');
           set(fig,'units','pixels');
  figpos = get(fig,'position');
           set(fig,'units',figuni);


  ButtonHigh = ceil(ud.PixelHight*ud.HighPadding) + 2*ud.BorderWidth;

  ButtonHigh = ButtonHigh * ( ud.TagNumber > 0 );

  %--------------------------------------------------------------------------
  if strcmp(mode,'n')
  % Normalized

     pos = ud.NormPosition .* figpos([3 4 3 4]);

     pos([1 3]) =  ceil( pos([1 3]) );
     pos([2 4]) = floor( pos([2 4]) );
 

  %--------------------------------------------------------------------------
  else
  % Absolut

    %  'Position'  , [  Left  Bottom  -Right -Top  ] , default: [ 5  5 -5 -5]
    %                [  Left  Bottom   Width  High ]
    %                [ -Right -Top     Width  High ]


     po = ud.Position;

     pos = ( [ 1  1  0  0 ] + po ) .* ( po >= 0 );

     pos([3 4]) = pos([3 4]) + ...
                   ( figpos([3 4]) - pos([1 2]) + 1 + po([3 4]) ) .* ...
                   ( po([1 2]) >= 0 ) .* ( po([3 4]) <= 0 );

     pos([1 2]) = pos([1 2]) + ...
                   ( figpos([3 4]) - pos([3 4]) + 1 + po([1 2]) ) .* ...
                   ( po([1 2]) < 0 ) .* ( po([3 4]) > 0 );

     pos([1 3]) =  ceil( pos([1 3]) );
     pos([2 4]) = floor( pos([2 4]) );
 
             
  end
  %--------------------------------------------------------------------------

     pos(4) = pos(4);

     ud.PixelPosition = pos;

  % Frame

     posf = pos;
     posf(4) = posf(4) - ButtonHigh + ud.BorderWidth * ( ButtonHigh > 0 );  % !!!

     
  % Buttons & Hide

  n = ud.TagNumber;

  if n > 0

     % Buttons

     pos(1) = pos(1)-0;
     pos(3) = pos(3)+0;

     posb = [ pos(1)          pos(2)+pos(4)-ButtonHigh  ...
             floor(pos(3)/n)  ButtonHigh                    ];

     posb = posb(ones(1,n),:);

     dw =  pos(3) - n * posb(1,3);

     ii = ( n-dw+1 : ud.TagNumber );

     posb(ii,3) = posb(ii,3) + 1;

     posb(2:n,1) = posb(2:n,1)+cumsum(posb(1:n-1,3));

  
     % HideHandle

     posh = posb;
     posh(:,1) = posh(:,1) + ud.BorderWidth; 
     posh(:,2) = posf(2)+posf(4)-ud.BorderWidth;
     posh(:,3) = posh(:,3) - 2 * ud.BorderWidth; 
     posh(:,4) = 1 * ud.BorderWidth;

  end
  % n > 0

  ud0.TAG_FRAME = ud;

     posf([3 4]) = max( posf([3 4]) , 1 );

     set(ud.FrameHandle , 'units'     , 'pixels' , ...
                          'position'  , posf     , ...
                          'userdata'  , ud0            )

   for ii = 1 : n

     posb(ii,[3 4]) = max( posb(ii,[3 4]) , 1 );

     set(ud.ButtonHandle(ii) , 'units'    , 'pixels' , ...
                               'position' , posb(ii,:)    )

     posh(ii,[3 4]) = max( posh(ii,[3 4]) , 1 );

     set(ud.HideHandle(ii) , 'units'    , 'pixels' , ...
                             'position' , posh(ii,:)       )

   end


  % Resize Children

    for ii = 1 : n

      ch = get(ud.HideHandle(ii),'userdata');

      for jj = 1 : size(ch.Handle,1)

         tag_frame(ud.FrameHandle,'ResizeChildren',ch.Handle(jj),ch.Position{jj},mode);

      end

    end

    for jj = 1 : size(ud.FrameChildren.Handle,1)

       tag_frame(ud.FrameHandle,'ResizeChildren', ...
            ud.FrameChildren.Handle(jj),ud.FrameChildren.Position{jj},mode);

    end

   %-------------------------------------------

   out = { ud.PixelPosition  ButtonHigh-ud.BorderWidth };  % !!!

   n = min(Nout,size(out,2));

   varargout(1:n) = out(1:n);


%*********************************************************
case 'ACTIVATE'

  clb = 0;   % CallBack

  if ~isempty(VarIn)

    nr = VarIn{1};

    ok = ( isnumeric(nr) &  ( prod(size(nr)) == 1 ) );
    if ok
       ok = ( any( nr == ( 0 : ud.TagNumber ) )  |  isnan(nr) );
       if ok
          if nr == ud.ActiveNr
             return
          end
       else
          hh = [ ud.FrameHandle ; ud.ButtonHandle ];
          ok = any( nr == hh );
          if ok
             nr = find( nr == hh );
             nr = nr(1)-1;
          end
       end
    end

    if ~ok

      Msg = [ Msg0  'Input must be a ButtonNumber or ZERO ' nl ...
                    ' or a ButtonHandle or TagFrameHandle.' ];
      return

    end

    if prod(size(VarIn)) > 1
       clb = VarIn{2};
    end

  else
     
    nr = ud.ActiveNr;

  end


    if isnan(nr)  

       sets = { 'off'  'off' };

    else

       ud.ActiveNr = nr;
 
       sets = { 'off'  ud.Visibility };

    end

    ud0.TAG_FRAME = ud;

    set(ud.FrameHandle  , 'visible'  , sets{2} , ...
                          'userdata' ,  ud0       );
    set(ud.ButtonHandle , 'visible'  , sets{2}       );


    %-----------------------------------------

    for ii = 1 : ud.TagNumber

      set( ud.HideHandle(ii) , 'visible' , sets{1+(ii==nr)} );

      ch = get(ud.HideHandle(ii),'userdata');

      ch = ch.Handle;

      set( ch , 'visible' , sets{1+(ii==nr)} );


      cht = get(ch,'tag');

        %---------------------------------------
        % Check for SelList

        tt = find( strcmp(cht,'SEL_FRAME') );

        for jj = tt(:)'

          sel_list(ch(jj),'visible',sets{1+(ii==nr)});
  
        end

        %---------------------------------------
        % Check for Tabular

        tt = find( strcmp(cht,'TAB_FRAME') );

        for jj = tt(:)'

          tab_list(ch(jj),'visible',sets{1+(ii==nr)});
  
        end

        %---------------------------------------
        % Check for TagFrame

        tt = find( strcmp(cht,'TAG_FRAME') );

        for jj = tt(:)'
    
          tag_frame(ch(jj),'visible',sets{1+(ii==nr)});
  
        end

        %---------------------------------------
        % Check for MsgLogo

        tt = find( strcmp(cht,'LOGO_FRAME') );

        for jj = tt(:)'
    
          msg_logo(ch(jj),'visible',sets{1+(ii==nr)});
  
        end

        %---------------------------------------
        % Check for MsgList

        tt = find( strcmp(cht,'MESSAGE_FRAME') );

        for jj = tt(:)'
    
          msg_list(ch(jj),'visible',sets{1+(ii==nr)});
  
        end


    end


    %-----------------------------------------
    % Check for activated Frame

    ch = ud.FrameChildren.Handle;

    set( ch , 'visible' , sets{1+(0==nr)} );


    cht = get(ch,'tag');


      %---------------------------------------
      % Check for Tabular

      tt = find( strcmp(cht,'TabFrame') );

      for jj = tt(:)'

         tab_list(ch(jj),'visible',sets{1+(0==nr)});

      end

      %---------------------------------------
      % Check for TagFrame

      tt = find( strcmp(cht,'TAG_FRAME') );

      for jj = tt(:)'
   
         tag_frame(ch(jj),'visible',sets{1+(0==nr)});
  
      end



    %-----------------------------------------
    % CallBack

     if isequal(clb,1)  &  ~isempty(ud.CBFcn)
        try
           feval( ud.CBFcn{:} , 'Activate' , 1 , nr );
        catch
           fprintf([nl Msg0 'Error call CBFcn.' nl lasterr nl ]);
        end
     end

    %-----------------------------------------

    out = { nr };

    n = min(Nout,size(out,2));

    varargout(1:n) = out(1:n);


%*********************************************************
case 'ADD'

  if ( prod(size(VarIn)) < 2 )

      Msg = [ Msg0  'Not enough Input Arguments.' ];

      return

  end

  %----------------------------------------------------------
  % Check 1. Input

    nr = VarIn{1};

    ok = ( isnumeric(nr) &  ( prod(size(nr)) == 1 ) );
    if ok
       ok = any( nr == ( 0 : ud.TagNumber ) );
       if ~ok
          hh = [ ud.FrameHandle ; ud.ButtonHandle ];
          ok = any( nr == hh );
          if ok
             nr = find( nr == hh );
             nr = nr(1)-1;
          end
       end
    end

    if ~ok

      Msg = [ Msg0  '1. Input must be a ButtonNumber or ZERO ' nl ...
                    ' or a ButtonHandle or TagFrameHandle.' ];
      return

    end

  %----------------------------------------------------------
  % Check 2. Input

    type = VarIn{2};

    ok = ( ischar(type) & ~isempty(type) & ...
           ( prod(size(type)) == size(type,2) ) );
    if ok
       type = lower(type);
       typs = {'uicontrol' 'sel_list'  'tab_list' 'tag_frame' 'msg_list' 'msg_logo'};
       ok = any(strcmp(type,typs));
    end

 
    if ~ok

      Msg = [ Msg0  '2. Input must be ' ...
                    '''uicontrol'', ''tab_list'', ''tag_frame'', ' ...
                    '''sel_list'' , msg_list'' or ''msg_logo'' .'];
      return

    end

  %----------------------------------------------------------
  % Check other Input

   VarIn = VarIn(3:end);

   if ~isempty(VarIn)

      VarIn = VarIn(:);

      is_logo = strcmp(type,'msg_logo');

      ok = is_logo;

      if ~ok

        ok = ( mod(size(VarIn,1),2) == 0 );

        if ok
  
           VarIn = reshape(VarIn,2,size(VarIn,1)/2);

              ok = iscellstr(VarIn(1,:));

        end

      
        if ~ok
  
          Msg = [ Msg0  'Following Inputs must be ' ...
                        'Property-Value-Pairs, Properties must be Strings.'];
          return

        end
 
      end

   end

  %----------------------------------------------------------
  % Build Control

   switch type

     %-------------------------------------------------------------
     case 'uicontrol'

       % Search for "position"
       if ~isempty(VarIn)

         VarIn0 = VarIn;

         ii = find( strcmp( lower(VarIn(1,:)) , 'position' ) );

         if ~isempty(ii)
           for jj = ii(:)'
             pos         = VarIn{2,jj};
             pos([3 4])  = abs(pos([3 4]));
             pos([3 4])  = pos([3 4]) + ( pos([3 4]) == 0 );
             VarIn{2,jj} = pos;
           end
         end
         
       end

       try
          h = uicontrol('parent',fig,VarIn{:});
       catch
          Msg = [ Msg0 'Error call UICONTROL.' nl lasterr ];
          return
       end

       if ~isempty(VarIn)

         VarIn = VarIn0;

       end

     %-------------------------------------------------------------
     case 'sel_list'

       try
         [Msg,h,ht] = sel_list(fig,'New',VarIn{:});
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)

          Msg = [ Msg0 'Error call SEL_LIST.' nl Msg ];
          return

       end

     %-------------------------------------------------------------
     case 'tab_list'

       try
         [Msg,h] = tab_list(fig,VarIn{:});
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)

          Msg = [ Msg0 'Error call TAB_LIST.' nl Msg ];
          return

       end


     %-------------------------------------------------------------
     case 'tag_frame'

       try
         [Msg,h,hb,p,bh] = tag_frame(fig,'New',VarIn{:});
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)
          Msg = [ Msg0 'Error call TAG_FRAME( New ).' nl Msg ];
          return

       end

       % Add Root to Buttons

         if nr == 0
            hr = ud.FrameHandle;
         else
            hr = ud.ButtonHandle(nr);
         end

         for hh = hb(:)'
           bud = get( hh , 'userdata' );
           bud.Root = hr;
           set( hh , 'userdata' , bud );
         end
  
       % Parent to FrameHandle

         hud = get( h , 'userdata' );
         hud.Parent = hr;
         hud.Root   = get( hr , 'Parent' );
       
         set( h , 'userdata' , hud );

     %-------------------------------------------------------------
     case 'msg_list'

       try
         [Msg,h] = msg_list(fig,'New',VarIn{:});
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)

          Msg = [ Msg0 'Error call MSG_LIST( New ).' nl Msg ];
          return

       end

     %-------------------------------------------------------------
     case 'msg_logo'

       try
         [Msg,hl,h,ht] = msg_logo(fig,'New',VarIn{:});
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)

          Msg = [ Msg0 'Error call MSG_LOGO( New ).' nl Msg ];
          return

       end

   end


  %----------------------------------------------------------
  % Get Position

   switch type

     %-------------------------------------------------------------
     case 'uicontrol'

       if isempty(VarIn)

         pos = get(h,'position');

       else

         ii = find( strcmp( lower(VarIn(1,:)) , 'position' ) );

         if isempty(ii)
           pos = get(h,'position');
         else
           pos = VarIn{2,ii(end)};
         end
         
       end

     %-------------------------------------------------------------
     case 'sel_list'
   
       tud = get(h,'userdata');
    
       pos = tud.Position;

     %-------------------------------------------------------------
     case 'tab_list'
   
       tud = get(h,'userdata');
    
       pos = tud.position;

     %-------------------------------------------------------------
     case 'tag_frame'
   
       tud = get(h,'userdata');
    
       pos = tud.TAG_FRAME.Position;
   
     %-------------------------------------------------------------
     case 'msg_list'
   
       tud = get(h,'userdata');
    
       pos = tud.PixelOffset;
   
     %-------------------------------------------------------------
     case 'msg_logo'
   
       tud = get(h,'userdata');
    
       pos = tud.Offset;

   end


  %-----------------------------------------
  % Resize

    Msg = tag_frame( ud.FrameHandle,'ResizeChildren',h,pos,'absolut',nr);

    if ~isempty(Msg)

        delete(h);

        return

    end

  %-----------------------------------------
  % Add new Children
 
    if nr == 0

      ud.FrameChildren.Handle   = cat( 1 , ud.FrameChildren.Handle , h );
      ud.FrameChildren.Position = cat( 1 , ud.FrameChildren.Position , {pos} );

      ud0.TAG_FRAME = ud;

      set( ud.FrameHandle , 'userdata' , ud0 );

    else

      ch = get(ud.HideHandle(nr),'userdata');
 
      ch.Handle = cat(1,ch.Handle,h);

      ch.Position = cat( 1 , ch.Position , {pos} );

      set( ud.HideHandle(nr) , 'userdata' , ch );

    end


  %-----------------------------------------

    if strcmp(type,'msg_logo');
       out = { hl h ht };        % MessageListHandle, LogoFrameHandle, LogoTextHandle
    elseif strcmp(type,'msg_logo');
       out = { h ht };           % FrameHandle TextHandle
    elseif strcmp(type,'tag_frame')
       tud = get(h,'userdata');
       out = { h  hb tud.TAG_FRAME.PixelPosition bh }; % TagFrameHandle TagButtonHandle
    else
       out = { h };
    end

    n = min(Nout,size(out,2));

    varargout(1:n) = out(1:n);


%*********************************************************
case 'RESIZECHILDREN'

  if prod(size(VarIn)) < 2

    Msg = [ Msg0 'ChildrenHandle and Position required.' ];
    return
 
  end


  h = VarIn{1};

%  ok = ( isnumeric(h)  &  ( prod(size(h)) == 1 ) );
%  if ok 
%     ok = ishandle(h);
%     if ok
%        ok = strcmp( get(h,'type') , 'uicontrol' );
%     end
%  end
%  if ~ok
%      Msg = [ Msg0  'Input must be a ChildrenHandle.' ];
%      return
%  end


  pos = VarIn{2};

  mode = lower(VarIn{3}(1));

  if prod(size(VarIn)) < 4

     nr = ud.ActiveNr;

  else

    nr = VarIn{4};

%    ok = ( isnumeric(nr) &  ( prod(size(nr)) == 1 ) );
%    if ok
%       ok = any( nr == ( 0 : ud.TagNumber ) );
%       if ~ok
%          hh = [ ud.FrameHandle ; ud.ButtonHandle ];
%          ok = any( nr == hh );
%          if ok
%             nr = find( nr == hh );
%             nr = nr(1)-1;
%          end
%       end
%    end
%    if ~ok
%      Msg = [ Msg0  '2. Input must be a ButtonNumber or ZERO ' nl ...
%                    ' or a ButtonHandle or TagFrameHandle.' ];
%      return
%    end

  end

 
  %-----------------------------------------------------
  % FramePosition

    % FigurePosition in Pixels

     figuni = get(fig,'units');
              set(fig,'units','pixels');
     figpos = get(fig,'position');
              set(fig,'units',figuni);

    fpos        = get(ud.FrameHandle,'position');

    opos        = zeros(1,4);  % [  Left  Bottom  -Right -Top  ]
    opos([1 2]) = fpos([1 2]) - 1 + ud.BorderWidth;
    opos([3 4]) = - ( figpos([3 4]) - ( fpos([1 2]) + fpos([3 4]) - 1 ) + ...
                      ud.BorderWidth );
 
             ff = figpos([3 4])./ud.FigPos0([3 4]);
 

  switch(get(h,'tag'));

     %-------------------------------------------------------------
     case 'TAG_FRAME'

    %  'Position'  , [  Left  Bottom  -Right -Top  ] , default: [ 5  5 -5 -5]
    %                [  Left  Bottom   Width  High ]
    %                [ -Right -Top     Width  High ]


       if strcmp(mode,'n');

          pos([1 3]) = pos([1 3]) * ff(1);
          pos([2 4]) = pos([2 4]) * ff(2);

       end
  
       pos([1 2]) = pos([1 2]) + opos([1 2]) .* ( pos([1 2]) >= 0 ) + ...
                                 opos([3 4]) .* ( pos([1 2]) <  0 );
       pos([3 4]) = pos([3 4]) + opos([3 4]) .* ( pos([3 4]) <= 0 );

       tud = get(h,'userdata');

       tud.TAG_FRAME.Position = pos;

       set(h,'userdata',tud);
 
       try
          Msg = tag_frame(h,'Resize','absolut');
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)

          Msg = [ Msg0 'Error call TAG_FRAME( Resize ).' nl Msg ];
        
       end


       if nr ~= ud.ActiveNr
 
          tag_frame( h , 'visible' , 'off' );

       end

     %-------------------------------------------------------------
     case 'SEL_FRAME'

       if strcmp(mode,'n');

          pos = pos .* ff;
 
       end

       pos([1 2]) = pos([1 2]) + opos([1 2]) .* ( pos([1 2]) >= 0 ) + ...
                                 opos([3 4]) .* ( pos([1 2]) <  0 );

       tud = get(h,'userdata');

       tud.Position = pos;

       set(h,'userdata',tud);

       try
          Msg = sel_list(h,'Resize','absolut');
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)

          Msg = [ Msg0 'Error call SEL_LIST( RESIZE ).' nl Msg ];

       end

       if nr ~= ud.ActiveNr
 
          sel_list( h , 'visible' , 'off' );

       end

     %-------------------------------------------------------------
     case 'TAB_FRAME'

       if strcmp(mode,'n');

          pos = pos .* ff;
 
       end

       pos([1 2]) = pos([1 2]) + opos([1 2]) .* ( pos([1 2]) >= 0 ) + ...
                                 opos([3 4]) .* ( pos([1 2]) <  0 );

       try
          Msg = tab_list(h,'position',pos);
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)

          Msg = [ Msg0 'Error call TAB_LIST( Position ).' nl Msg ];

       end

       if nr ~= ud.ActiveNr
 
          tab_list( h , 'visible' , 'off' );

       end

     %-------------------------------------------------------------
     case 'MESSAGE_FRAME'

        % pos = [  Left  Right  Top     ] 
        %       [  Left  Right -Bottom  ]
 
       if strcmp(mode,'n');

          pos([1 2]) = pos([1 2]) * ff(1);
          pos([ 3 ]) = pos([ 3 ]) * ff(2);

       end

       pos([1 2]) = pos([1 2]) + opos([1 3]).*[1 -1];
       pos([ 3 ]) = ( pos(3) - opos(4) ) * ( pos(3) >= 0  ) + ...
                    ( pos(3) - opos(2) ) * ( pos(3) <  0  );

       tud = get(h,'userdata');

       tud.PixelOffset = pos;

       set(h,'userdata',tud);

       try
          Msg = msg_list(h,'resize','absolut');
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)

          Msg = [ Msg0 'Error call MSG_LIST( Resize ).' nl Msg ];

       end

       if nr ~= ud.ActiveNr
 
          msg_list( h , 'visible' , 'off' );

       end


     %-------------------------------------------------------------
     case 'LOGO_FRAME'

       %  pos.MsgOffset  =  [ Left Right Top ]
       %  pos.LogoOffset =  [ Left  Right ]

       if strcmp(mode,'n');

          pos.Msg([1 2]) = pos.Msg([1 2]) * ff(1);
          pos.Msg([ 3 ]) = pos.Msg([ 3 ]) * ff(2);

          pos.Logo = pos.Logo * ff(1);

       end

       pos.Logo(1) = pos.Logo(1) + opos(1);
       pos.Msg(2)  = pos.Msg(2) - opos(3);
       pos.Msg(3)  = ( pos.Msg(3) - opos(4) ) * ( pos.Msg(3) >= 0  ) + ...
                        ( pos.Msg(3) - opos(2) ) * ( pos.Msg(3) <  0  );

       tud = get(h,'userdata');

       tud.Offset = pos;

       set(h,'userdata',tud);

       lud = get(tud.Handle.List,'userdata');
       lud.PixelOffset = pos.Msg;
       set(tud.Handle.List,'userdata',lud);

       try
          Msg = msg_logo(h,'resize','absolut');
       catch
          Msg = lasterr;
       end
       
       if ~isempty(Msg)

          Msg = [ Msg0 'Error call MSG_LOGO( Resize ).' nl Msg ];

       end

       if nr ~= ud.ActiveNr
 
          msg_logo( h , 'visible' , 'off' );

       end



     %-------------------------------------------------------------
     otherwise
     % uicontrol

    %  'Position'  , [  Left  Bottom  -Right -Top  ] , default: [ 5  5 -5 -5]
    %                [  Left  Bottom   Width  High ]
    %                [ -Right -Top     Width  High ]

       if strcmp(mode,'n');

          pos([1 3]) = pos([1 3]) * ff(1);
          pos([2 4]) = pos([2 4]) * ff(2);

       end

       poo  = pos;

       pos([1 2]) = pos([1 2]) + opos([1 2]) .* ( pos([1 2]) >= 0 ) + ...
                                 opos([3 4]) .* ( pos([1 2]) <  0 );
       pos([3 4]) = pos([3 4]) + opos([3 4]) .* ( pos([3 4]) <= 0 );


       po = pos;

       pos = ( [ 1  1  0  0 ] + po ) .* ( po >= 0 );

       pos([3 4]) = pos([3 4]) + ...
                     ( figpos([3 4]) - pos([1 2]) + 1 + po([3 4]) ) .* ...
                     ( po([1 2]) >= 0 ) .* ( po([3 4]) <= 0 );

       pos([1 2]) = pos([1 2]) + ...
                     ( figpos([3 4]) - pos([3 4]) + 1 + po([1 2]) ) .* ...
                     ( po([1 2]) < 0 ) .* ( po([3 4]) > 0 );

       uni = get(h,'units');

       pos([3 4]) = max( pos([3 4]) , 1 );

       set( h , 'units'    , 'pixels' , ...
                'position' , pos            );

       set(h,'units',uni);

       if nr ~= ud.ActiveNr
 
         set(h,'visible','off');

       end

 end


%***********************************************************************
case { 'VISIBLE'  'ENABLE' }

   
  ok = ~isempty(VarIn);

  if ok

    sets = VarIn{1};

    ok = ( ischar(sets) & ~isempty(sets) & ...
           ( prod(size(sets)) == size(sets,2) ) );
    if ok
       ok = any(strcmp(sets,{'on'  'off'}));
    end

  end

  if ~ok

    Msg = [ Msg0  'Input for VISIBLE must be ''on'' or ''off''.'];
    return

  end

  %-------------------------------------------------- 
  if strcmp(action,'VISIBLE')

    ud.Visibility = sets;

    ud0.TAG_FRAME = ud;

    set ( ud.FrameHandle , 'userdata' , ud0 );

    if strcmp(sets,'off')

       tag_frame( ud.FrameHandle , 'Activate' , NaN );
  
    else

       tag_frame( ud.FrameHandle , 'Activate' );

    end
     
    return

  end

  %--------------------------------------------------
  if size(VarIn,1) > 1

    nr = VarIn{2};

    ok = ( isnumeric(nr) &  ( prod(size(nr)) == 1 ) );
    if ok
       ok = ( any( nr == ( 0 : ud.TagNumber ) )  |  isnan(nr) );
       if ~ok
          hh = [ ud.FrameHandle ; ud.ButtonHandle ];
          ok = any( nr == hh );
          if ok
             nr = find( nr == hh );
             nr = nr(1)-1;
          end
       end
    end

    if ~ok

      Msg = [ Msg0  'Input must be a ButtonNumber or ZERO ' nl ...
                    ' or a ButtonHandle or TagFrameHandle.' ];
      return

    end

  else
     
    nr = ud.ActiveNr;

  end


  %---------------------------------------------------------
  if strcmp(sets,'off')  &  any( nr == [ ud.ActiveNr  0 ] )

         tag_frame( ud.FrameHandle , 'Activate' , 0 );

  end


  hh = ud.ButtonHandle;

  if nr > 0

     hh = hh(nr);

  end

  set( hh , 'enable' , sets );


%***********************************************************************
case 'DELETE'

   
  % Delete Children

    for ii = 1 : ud.TagNumber

      if ishandle(ud.HideHandle(ii))

        ch = get(ud.HideHandle(ii),'userdata');

        ch = ch.Handle;

        for h = ch(:)'

           try
             delete(h);
           end

        end

      end

    end


    hh = [ ud.ButtonHandle ; ud.HideHandle ; ud.FrameChildren.Handle ];

    for h = hh'

         try
           delete(h);
         end

    end
    
    ok = isempty(VarIn);
    if ~ok
       ok = ~isequal(VarIn{1},-1);
    end

    if ok
       set(ud.FrameHandle,'deletefcn','');
       delete(ud.FrameHandle)
    end

end

