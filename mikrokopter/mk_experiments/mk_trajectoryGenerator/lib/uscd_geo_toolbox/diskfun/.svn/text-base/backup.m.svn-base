function [msg,ARC_File,OPT_File,TMP_File,ini] = backup(varargin)

% BACKUP   creates a Backup-Archive of an Directory from selected Extensions
%
% BACKUP(SourcePath,ArchiveFileName,[options],'.ext', ... ,'|dir', ... )
% BACKUP(SourcePath,ArchiveFilePath,[options],'.ext', ... ,'|dir', ... )
%   
% <ext> select Extensions of Files to archive
% <dir> select SubDirectoryNames to exclude from Archive
%      
% If a ArchiveFilePath is given, it has to end with FILESEP: "/" !!!
%                                   
% Notes for using Extensions
%
%    give '.-' to archive Files without any Extensions
%    give '.*' to archive Files with all Extensions
%    use  '.r' to archive README-Files too
%         
%    if Extension contains only UpperCase or only LowerCase Letters
%       both cases will archived.
%
% The Options:
%
%   '-m method','-r','-nr','-l','-nl','-u','-nu','f','-nf', ...
%   '-d#','-K#','-M#','-T DateTime','-o options', ...
%   '+ Suffix' ,':s',':i',':a',':d',':c','*'
%
% '-m method'  Defines Type of Archive
%
%     method == ArchiveMethod:
%                'zip'
%                'tar'
%                'tgz'                 no Update, no Recurse
%                'gz'   --> 'tar.gz'   no Update
%                'bz'   --> 'tar.bz2'  no Update
%
%  '-r#' ==  recurse trough SourcePath, # = Number of RecursionDepth
%  '-nr' ==  no recurse
%  '-l'  ==  follow DirectoryLinks
%  '-nl' ==  no follow 
%  '-u'  ==  update ArchiveFile
%  '-nu' ==  no update
%  '-f'  ==  forcing overwrite existing ArchiveFile
%  '-nf' ==  no forcing
%
%  '-a'  ==  adds TempFile with FileList
%  '-na' ==  don't add TempFile with FileList
%
%  '-s'  ==  simulation mode, no archiv created!!!
%  '-ns' ==  no simulation mode
%
%  '-K#' ==  Max. Size of Files to archive in KBytes 
%  '-M#' ==  Max. Size of Files to archive in MBytes 
%
%  '-T DateTime' == Archive only Files newer then DateTime
%
%      DateTime  == 'YYYY MM DD  hh mm ss'    Numbers only
%                   'DD-Month-YYYY hh:mm:ss'  DateString
%
%  '-o options'  == additional Options for selected ArchiveMethod
%
%  '+ Suffix'    == Suffix for ArchiveFileName, appended with '_'
%                   a nonempty Suffix overwrites the automaticly Suffixes!
%
%  ':#'  Create Automaticly Archive, possible Modes (Options see below):
%
%        ':s' ==  ScriptArchive, Suffix "scr"
%        ':i' ==   ImageArchive, Suffix "img"
%        ':a' ==    ArchArchive, Suffix "arc"
%        ':d' ==    DataArchive, Suffix "dat"
%        ':c' ==   ComplArchive, Suffix "cpl"
%
%  '*'  == Creates a History of ArchiveUpdates, the History starts with 
%            the Number 000: [ArchiveFileName]_[Suffix]_000.***
%            each following Archive contains only newer Files.
%            To check the Date and Options, the OPT_File of the older Archives 
%            is required, to check the Date to renew the Files and the Options
%             and Extensions, which must be equal in the History.
%            The OPT_File is created by BACKUP, see below.
%             It looks like: [ArchiveFileName]_[Suffix]_###.***.bcp
%
%
%------------------------------------------------------------------------
%  defaults:  no recurse
%             no follow DirectoryLinks
%             no update
%             no force
%             no simulation
%        '-m method' ==  '-m gz'         % TAR.GZ-Archive
%             '-d#'  == '-dinf'          % Unlimited
%             '-K#'  == '-K512'          % 512-KBytes MaxSize
%             '-T*'  == '-T 00-Jan-0000'
%             '.ext' == '.m'             % M-Files only
%
% The Defaults are defined at the Begin of the BACKUP.
%
%------------------------------------------------------------------------
%
% If the Input SourcePath is missing, the actual Directory will archived.
%
% If the Input ArchiveFileName is missing, the FileName will build from the
%    - last SubDirectoryName of SourcePath and 
%    - the Suffix of the Option '+ Suffix' or 
%                 of actual DateTime: 'YYYYMMDD_hhmmss'.
% 
%   If the SourcePath locates in the HOME-Directory, the UserName will added.
%
%
%  example: >> backup('/work/project/matlab/')
%           creates an Archive:  matlab_YYYYMMDD_hhmmss.zip
%  
%  example: >> backup('/home/USER/matlab/toolbox','-mgz')
%           creates an Archive:  USER_toolbox_YYYYMMDD_hhmmss.tar.gz
%  
%
%------------------------------------------------------------------------
% Options for Automaticly Archives:
%
%  ':s' ==  ScriptArchive: '+scr','-M 1','-r',exdir,scr_ext
%  ':i' ==   ImageArchive: '+img','-M10','-r',exdir,img_ext
%  ':a' ==    ArchArchive: '+arc','-M10','-r',exdir,arc_ext
%  ':d' ==    DataArchive: '+dat','-M50','-r',exdir,dat_ext
%  ':c' ==   ComplArchive: '+cpl','-M10','-r',exdir,scr_ext,img_ext,...
%                                                   arc_ext,dat_ext
%
%            The InputOptions overwrite the automaticly Options!
%
%  exdir = '|tmp','|temp','|cache','|.xvpics','|backup'
%
%  scr_ext:
%   '.m','.c','.cpp','.f','.h','.py','.pas','.bas','.par','.dat','.asc', ...
%   '.ini','.cnf','.cfg','.cdl','.res','.opt','.log','.obj','.def','.mod', ...
%   '.txt','.text','.tex','.toc','.bib','.sty','.tab','.docu','.pdf, ...
%   '.inf','.info','.inp','.lst','.list','.edt','.rtf','.sdw','sxw','.doc', ...
%   '.csv','.xls','.hlp','.dll','.lnk','.kdelnk','.spk','.xpm','.ray', ...
%   '.shtml','.html','.htm','.rdf','.css','.cgi','.java','.js','.pl','.php', ...
%   '.mexlx','mexglx','.mexaxp','.mexsol','.mexsg64','.lx','.axp','.sol','.sg64', ...
%   '.cal','.cap','.grp','.sh','.bat','.sys','.conf','.config','.rc','.r','.-'
%
%  img_ext:
%   '.xpm','.gif','.jpg','.jpeg','.bmp','.png','.tif','.tiff','.xcf','.pcx', ...
%   '.ps','.ps1','.ps2','.eps','.psc','.epsc','.epsi','.fig','.xfig','.pdf', ...
%   '.ppm','.pbm','.fli','.mpg','.mpeg','.mov','.wmv','.avi'
%
%  arc_ext: '.zip','.tar','.tgz','.gz','.bz2', ...
%           '.arj','.jar','.rar','.rpm','.bcp'
%
%  dat_ext: '.mat','.dat','.asc','.cdf','.nc','.mnc','.hdf','.bin','.ray', ...
%           '.ctd','.mc','.rcm','.mtd','.adcp','.arg','.sami','.raw','.gps', ...
%           '.nav','.cnv','.btl','.cap','.db','.cal','.csv','.xls',
%           '.dxf','.stl','.0*'
%
%
% The Modes for Automaticly Archives are defined at the Begin of the BACKUP.
%
%
%------------------------------------------------------------------------
%
% BACKUP creates/updates an OptionFile: [ArchiveFile].bcp 
%        The OptionFile contains following Informations:
%
%         >F ArchiveFileName
%         >P SourcePath
%         >D DateTimeString
%         >C USER@HOST (ARCH)
%         >I Info about ArchiveContents
%         >B Size of ArchiveFile [bytes]
%         >H Date of previous ArchiveFile in History
%         >O VARARGIN of BACKUP
%
%
% [Msg,ArchiveFile,OptionFile,TemporaryFile,Info] = BACKUP( ... )
%
%  returns the FileNames and the Information: 
%
%    Info = [ NDirectory NFiles Bytes ArchiveBytes ]
%
% The Temporary File will deleted after creating the archive,
%  so normally the String for TemporaryFile should be empty!
%
%
% In SimulationMode the List of the Files to archive are returned:
%
% [Msg,Files,Bytes,Dates,Info] = BACKUP( ... , '-s' , ... )
%
%
%------------------------------------------------------------------------
%
% To REPEAT the BACKUP with the Options, saved in the OptionFile
%  in the Line ">O ... ", give the OptionFileName as first Input, 
%  surrounded with "<>":
%
%   BACKUP('<OptionFile>',[options],...)
%
%  example: 
%
%    pfad = fullfile(getenv('HOME'),'matlab');
%    name = 'my_matlab';
%
%    backup(pfad,name,'-mgz','-r','.m');    
%      % Create Archive:  my_matlab.tar.gz
%      %     OptionFile:  my_matlab.tar.gz.bcp
%
%    backup('<my_matlab.tar.gz.bcp>');  % Rewrites the Archive
%
%    backup('<my_matlab.tar.gz.bcp>','-mzip','.cnf','|mexcdf');
%      % Create a ZIP-Archive with same Options, 
%      % additional added CNF-Files, Exludes Directory "mexcdf"
%    
%------------------------------------------------------------------------
% NOTES for using BACKUP under WINDOWS
%
%  To create Archives under Windows, the Executable WTAR.EXE and GZIP.EXE
%   are provided in the "private"-Directory of BACKUP.
%
%   Additional the DLL-File CXGWIN.DLL is required.
%
%  With this tools only the Methods "gz" and "tar" are available.
%
%------------------------------------------------------------------------
% Examples:
%  
%  Source = fullfile( getenv('HOME') , 'matlab' );   % SourceDirectory
%  Destin = fullfile( getenv('HOME') , 'backup' );   % DestinationDirectory
%
%  backup( Source , fullfile(Destin,'backup') , '-mgz' , '-r' , '.m' , '.dat' );
%  % Creates in Directory $Destin an Archive:  backup.tar.gz
%                              an OptionFile:  backup.tar.gz.bcp   
%
%  backup( Source , fullfile(Destin,'backup') , '-mgz' , '-r' , ':scr'  );
%  % Creates in Directory $Destin an Archive:  backup_src.tar.gz
%                              an OptionFile:  backup_src.tar.gz.bcp   
%
%  backup( Source , '-mzip' , '-r' , 


%****************************************************************************

if 1 %%% isunix

   [msg,ARC_File,OPT_File,TMP_File,ini] = main(varargin{:});

else

   msg = 'UNIX-Enviroment required.';

   ARC_File = '';
   OPT_File = '';
   TMP_File = '';

   ini = zeros(0,4);  %  [ NDir NFiles Bytes ArchBytes ] 

end

%****************************************************************************
% Return

if ~isempty(msg) & ( nargout == 0 )

    fprintf(1,'\n%s\n\n',msg);

    clear msg

    return

end

%****************************************************************************
%
% Definition of defaults
%
%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [methods,no_recurse_method,no_update_method,isdos] = getmethods

% Returns valid Methods and Extensions

%        Method  MethodExt  ArchiveExt  Commands
methods = { 'gz'    '.tar'  '.tar.gz'   {'tar' 'gzip' }
            'tar'   '.tar'  '.tar'      {'tar'}
            'zip'   '.zip'  '.zip'      {'zip'}
            'tgz'   '.tgz'  '.tgz'      {'tar'}
            'bz'    '.tar'  '.tar.bz2'  {'tar' 'bzip2'}       };

no_recurse_method = { 'tgz' };

no_update_method = { 'tgz'  'gz'  'bz' };

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% WINDOWS: TAR.GZ only, using WTAR.EXE + GZIP.EXE
%          Programs in PRIVATE-Directory of BACKUP!

isdos = strcmp(upper(computer),'PCWIN');

if isdos

   pfd = fileparts(which(mfilename));
   prg = {'wtar.exe' 'gzip.exe' 'cygwin.dll'};
   for ii = 1 : size(prg,2)
       prg{ii} = fullfile(pfd,'private',prg{ii}); 
   end

   methods = { 'gz'  '.tar'  '.tar.gz'  prg 
               'tar' '.tar'  '.tar'     prg([1 3]) };

end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ini = getwiz

% Options for wizard, Automaticly Archives

% { Flag Suffix Options Extensions }

dir_ext = { 'tmp' ; 'temp' ; 'cache' ; '.xvpics' ; 'backup' };

ini = { 's'  'scr' { 'M 1' 'r' }  scr_ext  dir_ext
        'i'  'img' { 'M10' 'r' }  img_ext  dir_ext
        'a'  'arc' { 'M10' 'r' }  arc_ext  dir_ext
        'd'  'dat' { 'M50' 'r' }  dat_ext  dir_ext
        'c'  'cpl' { 'M10' 'r' }  {}       dir_ext  };

ini{end,4} = cat( 1 , ini{1:end-1,4} );     % Complete

 %----------------------------------------------
 function ext = scr_ext

    ext = { '.m'
            '.ini'
            '.cnf'
            '.cfg'
            '.cdl'
            '.res'
            '.opt'
            '.obj'
            '.def'
            '.log'
            '.par'
            '.txt'
            '.text' 
            '.tex'
            '.toc'
            '.bib'
            '.sty'
            '.tab'
            '.rtf'
            '.pdf'
            '.sdw'
            '.sxw'
            '.csv'
            '.xls'
            '.hlp'
            '.doc'
            '.docu'
            '.inf'
            '.info'
            '.lst'
            '.list'
            '.edt'
            '.mod'
            '.dat'
            '.asc'
            '.cap'
            '.grp'
            '.cal'
            '.ray'
            '.r'  
            '.c'
            '.cpp'
            '.f'
            '.h'
            '.spk'
            '.xpm' 
            '.shtml'
            '.html'
            '.htm' 
            '.rdf'
            '.css' 
            '.cgi' 
            '.java'
            '.js'
            '.py'
            '.pas'
            '.bas'
            '.pl'
            '.php'
            '.inp'
            '.mexlx'
            '.mexglx'
            '.mexaxp'
            '.mexsol'
            '.mexsg64'
            '.lx'
            '.axp'
            '.sol'
            '.sg64'
            '.dll'
            '.lnk'
            '.kdelnk'
            '.sh'
            '.bat'
            '.sys'
            '.conf'
            '.config'
            '.rc'
            '.-'       };

 %----------------------------------------------
 function ext = img_ext

    ext = { '.xpm'  
            '.gif'  
            '.jpg'
            '.jpeg'
            '.bmp'
            '.png'
            '.tif'
            '.tiff'
            '.xcf'
            '.pcx'
            '.ps'
            '.ps1'
            '.ps2'
            '.eps'
            '.psc'
            '.epsc'
            '.epsi'
            '.pdf'
            '.fig'
            '.xfig'
            '.ppm'
            '.pbm'
            '.fli'
            '.mpg'
            '.mpeg'
            '.mov'
            '.wmv'
            '.avi'  };


 %----------------------------------------------
 function ext = arc_ext

    ext = { '.tar'
            '.gz'
            '.tgz'
            '.zip'
            '.bz2' 
            '.arj' 
            '.rpm' 
            '.jar'
            '.rar'
            '.bcp' };


 %----------------------------------------------
 function ext = dat_ext

    ext = { '.mat'
            '.dat'
            '.asc'
            '.cdf'
            '.nc'
            '.mnc'
            '.hdf'
            '.bin'
            '.ray'
            '.ctd'
            '.mc'
            '.rcm'
            '.mtd'
            '.adcp'
            '.arg'
            '.sami'
            '.raw'  
            '.gps'
            '.nav'
            '.cnv'
            '.btl'
            '.cap'
            '.cal'
            '.db'
            '.csv'
            '.xls'
            '.dxf'
            '.stl'
            '.0*'  };


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ext = optext

% Returns OPT_File-Extension

ext = '.bcp';


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ext = defext

% Returns Default-Extension

ext = { '.m' };

%****************************************************************************
%
% End of Definition of defaults
%
%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,ARC_File,OPT_File,TMP_File,ini] = main(varargin)

% MAIN   Main-Function for BACKUP
  
%--------------------------------------------------------

Nin = nargin;

msg = '';
nl  = char(10);

ARC_File = '';
OPT_File = '';
TMP_File = '';

ini = zeros(0,4);  %  [ NDir NFiles Bytes ArchBytes ] 

%*********************************************************

isdos = strcmp(upper(computer),'PCWIN');

qt = '';
if isdos
   qt = '"';
end

fs  = filesep;

p0  = cd;

if isempty(p0)
   p0 = fs;
end

%****************************************************************************
% Basic Check Inputs

msg = checkbase(varargin);

%****************************************************************************
% Check for Repetition

if isempty(msg)  &  ~isempty(varargin)
   bcp = varargin{1};
   if size(bcp,2) >= 3
      if strcmp( bcp([1 end]) , '<>' )
         bcp = bcp( 2 : (end-1) );
         [msg,ARC_File,OPT_File,TMP_File,ini] = repeat(bcp,varargin(2:end));
         if ~isempty(msg)
             msg = cat(2,'Error call BACKUP:REPEAT',nl,msg);
          end
          return
       end
   end
end

%****************************************************************************
% Get and Check Inputs

simulate = 0;

if isempty(msg)

  try
     [msg,Source,File,ARC_File,TMP_File,cl,mode,prev,method,ext,exl,...
      comm,addtmp,simulate,recurse,follink,update,force,depth,bytes,day,option] = ...
      checkin(varargin,qt,isdos);
  catch
     msg = cat(2,'Error call BACKUP:CHECKIN',nl,lasterr);
  end

end

%****************************************************************************
% Special for WTAR.EXE

if isempty(msg)

   wtar = 'wtar.exe';
   nc   = size(comm{1},2);
   nw   = size(wtar,2);

   chk  = ( isdos & ( File(2) == ':' ) & ...
            strcmp( comm{1}(max(1,nc-nw+1):nc) , wtar ) );

   if chk
      if ~strcmp( Source([1 2]) , File([1 2]) )
          msg = sprintf('For using of %s the Drives of SourcePath (%s) and ArchiveFile (%s) must be the same.', ...
                         upper(wtar),upper(Source([1 2])),upper(File([1 2])));
      end
   end

   % Create relative Directory for File below !!!

end

%****************************************************************************
% Create FileList to Archive

file =  cell(0,1);
byte = zeros(0,1);
date =  cell(0,1);


if isempty(msg)

  try
     [msg,StartPfad,ini,file,byte,date] = makelist(Source,TMP_File,...
                                            ext,exl,depth,bytes,day, ...
                                      addtmp,simulate,recurse,follink,isdos);
  catch
      msg = cat(2,'Error call BACKUP:MAKELIST',nl,lasterr);
  end
 
end

%----------------------------------------
% Check TMP_File

if isempty(msg) & ( simulate == 0 )

  d = dir(TMP_File);

  ok = ~isempty(d);
  if ok
     ok = ~( d(1).bytes == 0 );
  end

  if ~ok
      msg = 'Temporary File is empty or does not exist.';
  end

end

%****************************************************************************
% Create Archive, using Method

if isempty(msg) 

   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   % Special for WTAR.EXE
   if chk
      File = relpath(File,StartPfad);
   end
   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   try
      command = makecomm(File,TMP_File,method,comm,update,option,qt,isdos);
   catch
      msg = cat(2,'Error call BACKUP:MAKECOMM',nl,lasterr);
   end

end

arc_ok = 0;
d0     = dir(ARC_File);

if isempty(msg)

   try
      msg = makearch(StartPfad,command,simulate);
   catch
      msg = cat(2,'Error call BACKUP:MAKEARCH',nl,lasterr);
   end

   arc_ok = simulate;

   if ~arc_ok

       d1 = dir(ARC_File);
       arc_ok = ( prod(size(d1)) == 1 );

       if ~arc_ok
          msg = [ 'Error creating Archive ' ARC_File nl 'File doesn''t exist.' ];
       end

       if arc_ok & ( prod(size(d0)) == 1 )
          arc_ok = ~isequal(d0,d1);
       end

       if arc_ok
          ini(4) = d1.bytes; 
       end

   end

end

cd(p0);


%****************************************************************************
% Delete TMP_File

if ~isempty(TMP_File) 

   if ( exist(TMP_File,'file') == 2 )
      try 
         delete(TMP_File)
      end
   end

   if ( exist(TMP_File,'file') == 2 )
      fprintf(1,'\nWarning: Can''t remove TMP-File:%s\n\n',TMP_File);
   else
      TMP_File = '';
   end

end

%****************************************************************************
% Info & OptionFile

if arc_ok

   fprintf(1,'\n%s successfull written\n\n',ARC_File);

   if ini(1) == 1
      dstr = 'Directory';
   else
      dstr = 'Directories';
   end

   if ini(2) == 1
      fstr = 'File';
   else
      fstr = 'Files';
   end

   form = cat( 2 , '%.0f ' , dstr , ', %.0f ' , fstr , ': ' , ...
                   '%.0f kB / %.1f MB  --->  %.0f kB / %.1f MB' );

   if update

     form = cat( 2 , 'update ' , form , ' total' );

     str = sprintf( form , ini(1) , ini(2) , ini(3)/2^10 , ini(3)/2^20 , ...
                           ini(4)/2^10 , ini(4)/2^20 );

   else

     form = cat( 2 , 'added ' , form , ' (%.0f%%)' );

     str = sprintf( form , ini(1) , ini(2) , ini(3)/2^10 , ini(3)/2^20 , ...
                        ini(4)/2^10 , ini(4)/2^20 , 100*ini(4)/(ini(3)+(ini(3)==0)) );

   end

   fprintf(1,'%s\n\n',str);

   %-------------------------------------------------------
   % Write OPT_File

      OPT_File = cat( 2 , ARC_File , optext );

      fprintf(1,'Write OptionFile: %s ...',OPT_File);

      if simulate
         fid = 1;   % STDOUT
         fprintf(fid,'\n\n');
      else
         fid = fopen(OPT_File,mode);
      end

      if fid == -1

         fprintf(1,'error, can''t open File.');

         OPT_File = '';

      else

         dat = datestr(datenum(cl(1),cl(2),cl(3),cl(4),cl(5),cl(6)),0);
         byt = sprintf('%11.0f',ini(4));

         usr = getenv('USER');
         hst = getenv('HOSTNAME');
         arc = getenv('ARCH');

         if ~isempty(hst)
             hst = cat(2,'@',hst);
         end
         if isempty(arc)
            arc = computer;
         end

         cmp = cat(2,usr,hst);
         if ~isempty(cmp)
             cmp = sprintf('%s (%s)',cmp,arc);
         else
             cmp = arc;
         end

         if isempty(varargin)
            opt = '';
         else
            opt = cat( 2 , '''' , strhcat(varargin,''',''') , '''' );
         end

         fprintf(fid,'>F %s\n',ARC_File);
         fprintf(fid,'>P %s\n',Source);
         fprintf(fid,'>D %s\n',dat);
         fprintf(fid,'>C %s\n',cmp);
         fprintf(fid,'>I %s\n',str);
         fprintf(fid,'>B %s\n',byt);
         fprintf(fid,'>H %s\n',prev);
         fprintf(fid,'>O %s\n',opt);

         if ~( fid == 1 )
             fclose(fid);
             fprintf(1,'ok.');
         end

      end 

      fprintf(1,'\n');

      if simulate
         fprintf(1,'Simulation finished\n');
      end

      fprintf(1,'\n');

else

   OPT_File = '';

   if simulate
      fprintf(1,'\nSimulation finished\n\n');
   end

end

if simulate
   ARC_File = file;
   OPT_File = byte;
   TMP_File = date;
end

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  command = makecomm(File,TMP_File,method,comm,update,option,qt,isdos);


switch method

  %------------------------------------------------------------------------
  case 'zip'

    opt = '';

    if update
       opt = '-u';
    end

    command = sprintf( '! %s%s%s "%s" %s %s -@ < "%s"' , qt , comm{1} , qt , ...
                       File , opt , option , TMP_File );  

  %------------------------------------------------------------------------
  case { 'tar'  'tgz'  'gz'  'bz' }

    opt = '';

    if update
       opt = '-uvf';
    elseif strcmp(method,'tgz')
       opt = '-czvf';
    else
       opt = '-cvf';
    end

    incl = { '-T'  '-I' };

    incl = incl{ 1 + strcmp(mexext,'mexsol') };

    command = sprintf( '! %s%s%s %s %s "%s" %s "%s"' , qt , comm{1} , qt , ...
                        opt , option , File , incl , TMP_File );

end

if prod(size(comm)) > 1
   % gz, bz
   
   if isdos
      command = cat( 2 , {command} , {sprintf('! %s%s%s "%s"',qt,comm{2},qt,File)} );
   else
      command = sprintf( '%s ; %s%s%s "%s"' , command , qt , comm{2} , qt , File );
   end

end


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  msg = makearch(StartPfad,command,simulate);

msg = '';
nl  = char(10);

try
  cd(StartPfad)
catch
  msg = ['Can''t change to Directory: ' StartPfad ];
  return
end

if ischar(command)
   command = {command};
else
   command = command(:)';
end

for cc = command

    fprintf(1,'\n%s\n\n',cc{1});

    if ~simulate
        eval(cc{1});
    end

end

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,StartPfad,ini,file,byte,date] = makelist( Pfad , TMP_File , ...
                                                 ext,exl,depth,bytes,day,...
                                         addtmp,simulate,recurse,follink,isdos);

msg = '';
nl  = char(10);

fs  = filesep;

ini    = zeros(1,4);  % [ NDir NFiles Bytes  ArchBytes ]
ini(4) = NaN;

file =  cell(0,1);
byte = zeros(0,1);
date =  cell(0,1);

%******************************************************************

[pre,StartPfad] = subname(Pfad);

Pfad = cat(2,Pfad,fs(1:(end*(~strcmp(Pfad(end),fs)))));

if ~isempty(pre)
   pre  = cat(2,pre,fs(1:(end*(~strcmp(pre(end),fs)))));
end

%******************************************************************
% Get DirectoryContents, write to TMP-File

%------------------------------------------------------------------
% Create TMP_File

if simulate

   fid = 1;  % stdout

else

   fid = fopen(TMP_File,'wt+');

   if fid == -1
      TMP_File = '';
      msg = [ 'Can''t open TMP_File ' TMP_File ' for FileList.' ];
      return
   end

end

%******************************************************************
% Get List

fprintf(1,'\nWrite Files of %s to Archive into:\n%s\n\n',Pfad,TMP_File);

try
  [ini,file,byte,date] = writelist(fid,Pfad,pre,fs,ext,exl,depth,bytes,day,...
                                    recurse,follink,ini,isdos);
catch
   msg = ['Error call BACKUP:MAKELIST:WRITELIST' nl lasterr];
end

if simulate
   return
end

if ( ini(2) > 0 ) & addtmp
   if isdos
      fprintf(fid,'%s\n',strrep(TMP_File,fs,'/'));
   else
      fprintf(fid,'%s\n',TMP_File);
   end
end

fclose(fid);

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ini,file,byte,date] = writelist(fid,Pfad,pre,fs,ext,exl,depth,bytes,day,...
                                         recurse,follink,ini,isdos)

nl = char(10);

file =  cell(0,1);
byte = zeros(0,1);
date =  cell(0,1);

% Pfad = cat(2,Pfad,fs(1:(end*(~strcmp(Pfad(end),fs)))));
% pre  = cat(2,pre,fs(1:(end*(~strcmp(Pfad(end),fs)))));

% Special for README-Files

rdm = { 'readme*'
        'Readme*'  };

is_rdm = any(strcmp(ext,'.r'));

if is_rdm
   ext = cat( 1 , ext(:) , rdm );
end

n = cell(0,1);  % Names

for e = ext(:)'

    d = dircont(Pfad,cat(2,'*',e{1}));

    if ~isempty(d)

       ok = ( ( cat(1,d{:,1}) == 0     )  & ...
              ( cat(1,d{:,3}) <= bytes )  & ...
              ( cat(1,d{:,5}) >= day   )        );

       d = d( find(ok) , : );

       nd = size(d,1);

       chk_rdm = ( is_rdm & ~isempty(n) & ~isempty(d) );
       if chk_rdm
          chk_rdm = any(strcmp(e{1},rdm));
       end

       if chk_rdm
       % Check for Duplicate README-Files

             ok = zeros(nd,1);
             for ii = 1 : nd
                 ok(ii) = ~any(strcmp(n,d{ii,4}));
             end

             d = d( find(ok) , : );
                           
       end

       if ~isempty(d)

          ini(1:3) = ini(1:3) + [ 1  size(d,1)  sum(cat(1,d{:,3})) ];

          n = cat( 1 , n , d(:,4) );  % Names

          ff = d(:,4);
          bb = cat(1,d{:,3});
          dd = d(:,2);

          d = d(:,[4 4 4]);  % FileNames

          d(:,1) = { pre };

          d(:,3) = { nl };

          d = permute(d,[2 1]);

          d = cat(2,d{:});

          if isdos
             d = strrep(d,fs,'/');
          end

          if fid == 1
             ff   = cat( 2 , pre(ones(size(ff,1),1),:) , char(ff) );
             ff(find(ff==fs)) = '/';
             file = cat(1,file,cellstr(ff));
             byte = cat(1,byte,bb);
             date = cat(1,date,dd);
          end

          fprintf(fid,'%s',d);

       end 

    end

end

if ~recurse | ( imag(depth) == real(depth) )
    return
end

depth = real(depth) + i * ( imag(depth) + 1 );

%----------------------------------------------------
% Recurse

d = dircont(Pfad,fs);

if ~isempty(d)

   bad = ( strcmp( d(:,4) , '..' )  |  strcmp( d(:,4) , '..' ) );

   d = d( find( ( cat(1,d{:,1}) == 1 ) & ( ~bad ) ) , : );

   if ~isempty(d)

      for ii = 1 : size(d,1)

          pf = cat( 2 , Pfad , d{ii,4} , fs );
          pr = cat( 2 , pre  , d{ii,4} , fs );

          ok = ( follink | isdos );

          if ~ok
              ok = ~islink(pf);
          end

          if ~isequal( exl , {char(0)} )

             chk2 = ( strcmp(exl,lower(exl)) | strcmp(exl,upper(exl)) );

             ok = ( ok & ~any(strcmp(exl,d{ii,4})) );

             if any(chk2)

                jj = find(chk2);

                ok = ( ok & ~any(strcmp(lower(exl(jj)),lower(d{ii,4}))) );
           
             end

          end

          if ok
             [ini,ff,bb,dd] = writelist(fid,pf,pr,fs,ext,exl,depth,bytes,day,...
                                      recurse,follink,ini,isdos);
              file = cat(1,file,ff);
              byte = cat(1,byte,bb);
              date = cat(1,date,dd);
          end

       end

    end
   
end


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function msg = checkbase(Vin);

% Basic Check Inputs

msg = '';

if isempty(Vin)
   Vin = cell(0,1);
end

ok = ( iscellstr(Vin) | isempty(Vin) );
if ok & ~isempty(Vin)
   try
     cat(2,Vin{:});
   catch
     ok = 0;
   end
end

if ~ok
    msg = 'Inputs must be Strings.';
end

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,Source,File,ARC_File,TMP_File,cl,mode,prev, ...
          method,ext,exl,comm,addtmp,simulate,recurse,follink,update, ...
          force,depth,bytes,day,option] = checkin(Vin,qt,isdos);


% CHECKIN  Check of Inputs 
%

msg = '';
nl  = char(10);

%*********************************************************
% Defaults

Source = '';
File   = '';

ARC_File = '';
TMP_File = '';

mode = 'w+';

cl   = zeros(1,6);  % [ YY MM DD hh mm sss ]

prev = '';

%---------------------------------------------------------
% String for DateTime

form  = '_%4.0f%2.2d%2.2d_%2.2d%2.2d%2.2d';

cl  = round(clock);

def = sprintf(form,round(clock));

%*********************************************************
% Get Options

[ok,v0,v1,ext,exl,use_prev,suf] = chkopt(Vin);

%*********************************************************
% Get Options

[msg,method,mext,aext,comm,addtmp,simulate,recurse,follink,update,force, ...
 depth,bytes,day,option] = getopt(v1,isdos);

if ~isempty(msg)
   return
end

%*********************************************************

if ~isempty(suf)
   ok = checkname(suf,'file');
   if ~ok
      msg = 'Invalid Characters in Suffix.';
      return
   end
end

%*********************************************************
% Get Source and File

[msg,Source,Target,File,auto] = getfile(v0);

if ~isempty(msg)
    return
end

%*********************************************************
% Check File for Extension and Suffix

File = File( 1 : ( end - 1*strcmp(File(end),'.') ) );

ee = { mext aext };
ee = ee( 1 : 1+(~strcmp(ee{1},ee{2})) );

for e = [ ee {suf} ]

   i1 = size(File,2);
   i0 = max( 1 , i1-size(e{1},2)+1 );

   i1 = i1 + ( i0-1 - i1 ) * strcmp(File(i0:i1),e{1});

   File = File(1:i1);

end

%------------------------------------------------------------------
% Add Suffix 

File = cat( 2 , File , suf );

%*********************************************************
% Check for Previous

if use_prev

    nr = 0;

   [nd,opt,form] = getprev(Target,File,aext);

   % nd = [ Nr Day Ok ]

   if ~isempty(nd)

      [nd,opt] = checkprev(nd,opt,cl,ext,recurse,follink,depth,bytes,option);

      ii = find(nd(:,3));

      ok = ~isempty(ii);
      if ok
         ii = min(ii);
         %     Highest Number  |  Newest Date
         ok = ( ( ii == 1 )  |  ( nd(ii,2) == max(nd(1:ii,2)) ) );
      end

      %------------------------------------------------------------
      if ~( ok | force )

         txt = [ nl 'Inconsistent previous ArchiveHistory.' ];

         int = [ txt nl '     <C>ontinue | <R>enew History | {<A>bort} : ' ];

         s = input(int,'s');

         s = lower( cat( 2 , s , 'a' ) );

         if ~any(strcmp(s(1),{'c' 'r'}))
            msg = [ txt ' Aborted.' ];
            return
         end

         if strcmp( s(1) , 'r' )
         % Renew
           nd(:,3) = 0;
         end

         force = 1;

      end 
      %------------------------------------------------------------

      ii = find(nd(:,3));

      if ~isempty(ii)
          ii = min(ii);
          nr = nd(ii,1) + 1;
         day = nd(ii,2);
        prev = datestr(day);
      end

   end

   File = cat( 2 , File , sprintf(form,nr) );

   update = 0;   % !!!!!!


elseif auto & isempty(suf)

   File = cat( 2 , File , def );

   def = '';  % To not use for TMP_File !!!

end

%------------------------------------------------------------------
% Add TargetDirectory & Extension

File = cat( 2 , Target , File , mext );

%------------------------------------------------------------------
% Check File for Directory

for e = ee
   f = cat( 2 , File , e{1} );
   if exist(f,'dir') == 7
      msg = ([ 'A Directory with Name of File: ' f ' allready exist.' ])
      return
   end
end

%-------------------------------------------------------------------

i1 = size(File,2) - size(mext,2);

ARC_File = cat( 2 , File(1:i1) , aext );

TMP_File = cat( 2 , File(1:i1), def , '.tmp' );

%-------------------------------------------------------------------
% Check for Existing File

exst = ( exist(ARC_File,'file') == 2 );

if exst 

  %-----------------------------------------------------
  % Check Archive to UpDate
  if update 

      switch method
        case 'zip'
          command = sprintf('%s%s%s -T "%s"',qt,comm{1},qt,ARC_File);
        case 'tar'
          command = sprintf('%s%s%s -tvf "%s"',qt,comm{1},qt,ARC_File);
      end

      fcn = 'unix';
      if isdos
         fcn = 'dos';
      end

      [s,w] = feval(fcn,command);

      update = ( s == 0 );

      if ~( update | force )

            txt = [  'File ' ARC_File ' is not a '  upper(method) ...
                     '-Archive.' ];

            int = [ nl upper(fcn) ': ' command nl w nl txt nl ...
                    '     <O>verwrite | {<A>bort} : ' ];

            s = input(int,'s');

            s = cat( 2 , s , 'a' );

            if ~strcmp(lower(s(1)),'o')
               msg = [ txt ' Aborted.' ];
               TMP_File = '';
               ARC_File = '';
               return
            end

      end
      % force

  %-----------------------------------------------------
  elseif ~force

      txt = [ nl 'File ' ARC_File ' allready exist.' ];

      int = [ txt nl  '     <O>verwrite | {<A>bort} : ' ];

      s = input(int,'s');

      s = cat( 2 , s , 'a' );

      if ~strcmp(lower(s(1)),'o')
         msg = [ txt ' Aborted.' ];
         TMP_File = '';
         ARC_File = '';
         return
      end

  end
  % update

end

%-------------------------------------------------------------------
% Try to open File

update = ( update  &  exst );

mode = { 'w+'  'a' };

mode = mode{ ( 1 + update ) };
 
fid = fopen(ARC_File,mode);
 
if fid == -1
   if update
      msg = (['Can''t open File: '  ARC_File ]);
      return
   end
else
   fclose(fid);
end

if ~update
   if exist(ARC_File,'file') == 2
      try,delete(ARC_File);end
   end
end


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,Source,Target,File,auto] = getfile(v);

% Returns SourceDirectory and FileName
%

msg = '';

nl = char(10);

Source = '';
Target = '';
File   = '';

auto = 0;

fs = filesep;

p0 = cd;

if isempty(p0)
   p0 = fs;
end

%*********************************************************
% Get Source and File

   nv = prod(size(v));

   if nv >= 1
      Source = v{1};
   end

   if nv >= 2
      File = v{2};
   end

if ~isempty(File)
   ok = checkname(File,'dir');
   if ~ok
      msg = 'Invalid Characters in ArchiveFileName.';
      return
   end
end
 
%*********************************************************
% Check Source

if isempty(Source)
   Source = p0;
end

try
  cd(Source)
catch
  msg = (['Can''t change to SourceDirectory: '  Source  ]);
  return
end
  
if ~strcmp(Source,fs)
    Source = cd;
end

cd(p0);


Source = cat( 2 , Source , fs(1:(end*(~strcmp(Source(end),fs)))) );

d = dir(Source);

if isempty(d)
   msg = (['Can''t read SourceDirectory: ' Source ]);
   return
end

%*********************************************************
% DefaultFileName from Source

%------------------------------------------------------
% Check for HOME-Directory

File0 = subname(Source,'backup');

File0 = File0( 1+strcmp(File0(1),'.') : end );

usr = getenv('USER');
hm  = getenv('HOME');

%------------------------------------------------------
% Check for HOME-Directory

ok = ~isempty(hm);

if ok

   n = min( size(hm,2) , size(Source,2) );

  ok = strcmp( Source(1:n) , hm );

end

%------------------------------------------------------
% Check for USER-Name

ok = ( ok & ~isempty(usr) );

if ok

   n = min( size(usr,2) , size(File0,2) );

  ok = ( ok & ~strcmp( File0(1:n) , usr ) );

 usr = cat( 2 , usr , '_' );

end

%------------------------------------------------------

File0 = cat( 2 , usr(1:end*ok) , File0 );

%*********************************************************
% Get Directory for FileName

% Check for File == "Directory/"

ok = ~isempty(File);
if ok
   ok = ( ( exist(File,'dir') == 7 ) & strcmp(File(end),fs) );
   if ok
      try
        ok = ~isempty(dir(File));
     catch
        ok = 0;
     end
   end
end

if ok
   Target = File;
   File   = File0;
   auto   = 1;
else
   [File,Target,auto] = subname(File,File0);
    File = File( 1+strcmp(File(1),'.') : end );
end

if isempty(Target)
   Target = p0;
end

try
   cd(Target)
catch
   msg = (['Can''t change to ArchiveDirectory: '  Target  ]);
   return
end

if ~strcmp(Target,fs)
   Target = cd;
end

cd(p0);

Target = cat( 2 , Target , fs(1:(end*(~strcmp(Target(end),fs)))) );

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,method,mext,aext,comm,addtmp,simulate,recurse,follink, ...
          update,force,depth,bytes,day,option] = getopt(Vin,isdos)

% Returns BACKUP-Options
%

msg = '';
nl  = char(10);

%**********************************************************
% Defaults

recurse = 0;
follink = 0;
update  = 0;
force   = 0;
simulate= 0;

addtmp  = 0;

depth   = NaN;
bytes   = 512 * 2^10;
day     =   0;

[methods,no_recurse_method,no_update_method,isdos] = getmethods;


mext   = methods(:,2);  % MethodExtension
aext   = methods(:,3);  % ArchiveExtension
comm   = methods(:,4);  % Commands

methods = methods(:,1);

method = methods{1};    % Default

option = '';


%**********************************************************

  msgT = '';
  msgD = '';
  msgK = '';
  msgM = '';

for ii = 1 : prod(size(Vin))

    v = Vin{ii};

    s2 = size(v,2);

    i2 = min( 2 , s2 );

     recurse = ( ( recurse | strcmp(v(1),'r') ) & ~strcmp(v(1:i2),'nr') );
     follink = ( ( follink | strcmp(v(1),'l') ) & ~strcmp(v(1:i2),'nl') );
     update  = ( ( update  | strcmp(v(1),'u') ) & ~strcmp(v(1:i2),'nu') );
     force   = ( ( force   | strcmp(v(1),'f') ) & ~strcmp(v(1:i2),'nf') );
     simulate= ( ( simulate| strcmp(v(1),'s') ) & ~strcmp(v(1:i2),'ns') );
     addtmp  = ( ( addtmp  | strcmp(v(1),'a') ) & ~strcmp(v(1:i2),'na') );

 
    %--------------------------------------------------------
    if ( s2 > 1 )

       v1 = v(1);
       v  = v(2:end);

       switch v1

          %---------------------------------------
          case 'T'

             try

               try
               % Try DateVec

                 d = zeros(1,6);

                 v = eval([ '['   v   ']' ]);

                 nv = size(v,2);
                 n  = min( nv , size(d,2) );

                 d(1:n) = v(1:n);

                 day    = datenum(d(1),d(2),d(3),d(4),d(5),d(6));

               catch
               % Try DateString

                 day = datenum( ger2eng(v) );

               end

             catch

               msgT = 'Invalid Format';

             end
  
          %---------------------------------------
          case 'r'

            eval(['depth = '  v  ';'],'msgD = lasterr;'); 

            if isempty(msgD)
               if isempty(depth)
                  depth = NaN;
               else
                  ok = ( isnumeric(depth) & ( prod(size(depth)) == 1 ) );
                  if ok
                     if ~isfinite(depth)
                         depth = NaN;
                     else
                         ok = ( ( mod(depth,1) == 0 ) & ( depth >= 0 ) );
                     end
                  end
                  if ~ok
                      msgD = 'Single positive integer or infinite required.';
                  end
               end
            end

          %---------------------------------------
          case 'K'

            eval(['bytes = '  v  '*2^10;'],'msgK = lasterr;'); 

            if isempty(msgK)
                  ok = ( isnumeric(bytes) & ( prod(size(bytes)) == 1 ) );
                  if ok
                     ok = isfinite(bytes);
                  end
                  if ~ok
                      msgK = 'Single finite numeric required.';
                  end
            end

          %---------------------------------------
          case 'M'

            eval(['bytes = '  v  '*2^20;'],'msgM = lasterr;'); 

            if isempty(msgM)
                  ok = ( isnumeric(bytes) & ( prod(size(bytes)) == 1 ) );
                  if ok
                     ok = isfinite(bytes);
                  end
                  if ~ok
                      msgM = 'Single finite numeric required.';
                  end
            end

          %---------------------------------------
          case 'o'

            option = rmblank(v,2);

          %---------------------------------------
          case 'm'

            method = rmblank(v,2);

       end
       % switch

    end
    % s2 > 1

end
% ii

  %---------------------------------------------------------------------------
  if ~any(strcmp(method,methods))

      msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                     'Method must be any of ' , strhcat(methods,', ') , '.' );    

  else

    if recurse & any(strcmp(method,no_recurse_method))

        msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                'Option ''-r'' (Recurse) not allowed for method ''' , method ,'''.' );    
 
    end

    if update & any(strcmp(method,no_update_method))

        msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                'Option ''-u'' (UpDate) not allowed for method ''' , method ,'''.' );    
 
    end

    mm = find( strcmp( method , methods ) );

    mext = mext{mm};
    aext = aext{mm};
    comm = comm{mm};

    for cc = comm(:)'

        if isdos
           s = ~any( exist(cc{1},'file') == [ 2  3 ] );
        else
           [s,w] = unix(['which ' cc{1} ]);
        end

        if ~( s == 0 )

            msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
             'Can''t find requested command ' , upper(cc{1}) , ' on machine.' );

        end

    end

  end   

  %---------------------------------------------------------------------------
  if ~isempty(msgD)
      msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                      'Invalid Value for Option ''-r''.' , nl ,  msgD );
  end

  %---------------------------------------------------------------------------
  if ~isempty(msgK)
      msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                      'Invalid Value for Option ''-K''.' , nl ,  msgK );
  end

  %---------------------------------------------------------------------------
  if ~isempty(msgM)
      msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                      'Invalid Value for Option ''-M''.' , nl ,  msgM );
  end

  %---------------------------------------------------------------------------
  if ~isempty(msgT)
      msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                      'Invalid Value for Option ''-T''.' , nl ,  msgT );
  end


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,v0,v1,ext,exl,prv,suf] = chkopt(Vin);

% Search and Separate  Options in Vin
%

v0 = {};
v1 = {};

ext = defext;
exl = {char(0)};

prv = 0;

suf = '';

n = prod(size(Vin));

ok = zeros(n,1);

if n == 0 
   return
end

wizini = getwiz;  % { Flag Suffix Options Extensions }

ext = '.';
exl = '|';
opt = '-';
suf = '+';
wiz = ':'; 
prv = '*';

fs = filesep;

for ii = 1 : n

  v = Vin{ii};

  ok(ii) = 3 * isempty(v);

  if ~ok(ii)
     
    s2  = size(v,2);
    s22 = ( s2 > 1 );

    is_prv = strcmp( v     , prv );
    is_wiz  = strcmp( v(1) , wiz );
    is_suf  = strcmp( v(1) , suf );
    is_opt  = strcmp( v(1) , opt );
    is_ext  = strcmp( v(1) , ext );
    is_exl  = strcmp( v(1) , exl );

    is_file = ~( is_prv | is_wiz | is_suf | is_opt );

    v2 = v(1+s22);

    is_wiz = ( is_wiz & s22 & any(strcmp(wizini(:,1),v2)) );
    is_suf = ( is_suf & s22 );
    is_opt = ( is_opt & s22 );
    is_exl = ( is_exl & s22 & ~strcmp(v2,fs) );
    is_ext = ( is_ext & s22 & ~any(strcmp(v2,{ '.' fs})) );

    is_file = ( is_file & ~( is_exl | is_ext ) );

    ok(ii) =   1 * is_opt + 2 * is_suf + 3 * is_file ...
             - 1 * is_ext - 2 * is_wiz - 3 * is_prv - 4 * is_exl;

    i0 = 1 + ( is_opt | is_suf | is_wiz | is_exl );

    Vin{ii} = rmblank( Vin{ii}( i0 : end ) , 2 );

  end

end

v0 = Vin( find( ok == 3 ) );
v1 = Vin( find( ok == 1 ) );

ext = Vin( find( ok == -1 ) );
exl = Vin( find( ok == -4 ) );

prv = any( ok == -3 );

%----------------------------------------------------
% Get Suffix

if any( ok == 2 )
   is_suf = max( find( ok == 2 ) );
   suf = Vin{is_suf};
   if ~isempty(suf)
      suf = cat( 2 , '_' , suf );
   end
else
   suf = '';
end

%----------------------------------------------------
% Get Wizard

if any( ok == -2 )
   is_wiz = find( ok == -2 );

   w = char(Vin(is_wiz));
   w = w(:,1);

  % Remove Duplicate Wizards
   [ww,si] = sort(w);
       bad = find( diff(double(ww),1,1) == 0 ) + 1;

   w(si(bad)) = [];

   nw = size(w,1);

   for ii = nw : -1 : 1

     is_wiz = find( strcmp( wizini(:,1) , w(ii) ) );
     v1  = cat( 1 , wizini{is_wiz,3}(:) , v1(:) );
     ext = cat( 1 , wizini{is_wiz,4}(:) , ext(:) );
     exl = cat( 1 , wizini{is_wiz,5}(:) , exl(:) );
     if isempty(suf) & ( nw == 1 )
        suf = cat( 2 , '_' , wizini{is_wiz,2} );
     end

   end

end


%----------------------------------------------------
% Check Extensions

if isempty(ext)

   ext = defext;

elseif any( strcmp( ext , '.*' ) )

   ext = { char(ones(1,0)) };

else

  % Remove Duplicate Extensions

  [ee,si] = sort(ext);

  ee = char(ee);

  bad = find( sum( ( diff(double(ee),1,1) == 0 ) , 2 ) == size(ee,2) ) + 1;

  ext(si(bad)) = [];

end

%----------------------------------------------------

if isempty(exl)
   exl = { char(0) };
end


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [nd,opt,form] = getprev(Target,File,aext);

% Search for Archives in UsePrevius-Mode: [Target]/[File]_###.ext

pre = '_';

form = cat(2,pre,'%3.3d');

i00  = size(File,2) + size(pre,2) + 1;  % StartIndex for ArchiveNumber, without '_' !!!

oext = optext;  % Extension for OPT_File

nd  = zeros(0,3);    % [ Nr Day Ok ]
opt =  cell(0,1);    % Options

%***************************************************************
% Look for Archives

d = dir(cat(2,Target,File,pre,'*',aext));

if ~isempty(d)
   d = d( find(~cat(1,d.isdir)) );
end

if isempty(d)
   return
end

%***************************************************************
% Check Archives

n = prod(size(d));

ref_in = { char(10) char(13) [char(32) char(9)] };

nd  = NaN*zeros(n,3);  % [ Nr  Date Ok ]
opt = cell(n,1);       % Options

opt(:) = { cell(0,1) };


%---------------------------------------------------------------

for ii = 1 : n

    %*************************************************************
    % Check ArchiveNumber

    nr = d(ii).name(i00:end);
     
    if ~isempty(nr)

       i1 = find( double(nr) == double('.') );

       if ~isempty(i1)

          i1 = min(i1) - 1;

          try, nr = eval(nr(1:i1)); catch, nr = []; end

          if ( isnumeric(nr) & ( prod(size(nr)) == 1 ) );
             if isfinite(nr) & ( nr >= 0 ) & ( mod(nr,1) == 0 )
                nd(ii,1) = nr;
             end 
          end

       end

    end

    %*************************************************************
    % Check for OPT-File

    ok = ~isnan(nd(ii,1));

    if ok

       o = dir( cat(2,Target,d(ii).name,oext) );

       ok = ( prod(size(o)) == 1 );
       if ok
          ok = ~o.isdir;
       end

    end

    %---------------------------------------------------------
    % Load OPT-File

    if ok

       if isempty(fileparts(o.name))
          o.name = cat(2,Target,o.name);
       end

       [msg,bb] = loadfile(o.name,1024,'char');

       ok = ( isempty(msg) & ~isempty(bb) );

       if ok
          bb = cat( 2 , bb , ref_in{1} );
       end

    end

    %---------------------------------------------------------
    % Get and Check Size

    if ok

       [msg,ref] = get_ref(bb,'>B',ref_in{:});

       ok = ( isempty(msg) & ~isempty(ref) );

       if ok
          try
             ref = eval(ref{end,2});
          catch
             ok = 0;
          end
       end
 
       if ok
          ok = ( ref == d(ii).bytes );
       end

    end

    %---------------------------------------------------------
    % Get Date

    if ok

       [msg,ref] = get_ref(bb,'>D',ref_in{:});

       ok = ( isempty(msg) & ~isempty(ref) );
 
       if ok
          try
             nd(ii,2) = datenum(ger2eng(ref{end,2})); 
          end
       end
 
       ok = ~isnan(nd(ii,2));

    end

    %---------------------------------------------------------
    % Get Options

    if ok

       [msg,ref] = get_ref(bb,'>O',ref_in{:});

       ok = ( isempty(msg) & ~isempty(ref) );

       if ok
          try
             ref = eval([ '{' ref{end,2} '}' ]); 
          catch
             ok = 0;
          end
       end

       if ok                
          msg = checkbase(ref);
          ok = isempty(msg);
          if ok
             opt{ii} = ref;
          end
       end

    end

    %*************************************************************
    % OK ???

    nd(ii,3) = ok;

end

%***************************************************************

ok = find( ~isnan(nd(:,1)) );  % Nr. Ok

nd  =  nd(ok,:);
opt = opt(ok);


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [nd,opt] = checkprev(nd,opt,cl,ext,recurse,follink,depth,bytes,option)

% Compare Extensions and Options of previous Files with actual

[h,si] = sort(-nd(:,1));   % Sort reverse by Nr !!!

nd  =  nd(si,:);
opt = opt(si);

%--------------------------------------------------------
% DateLimit: 01.01.2000 .. CurrentDate

d0 = datenum(2000,01,00);  
d1 = datenum(cl(1),cl(2),cl(3),cl(4),cl(5),cl(6));

ok = ( ( d0 < nd(:,2) ) & ( nd(:,2) < d1 ) );

nd(:,3) = ( nd(:,3) & ok );
 
%--------------------------------------------------------
% Analyse Options

good = find(nd(:,3));

for ii = good(:)'

    [ok,v0,v1,e] = chkopt(opt{ii});

    ok = isequal(e,ext);

    if ok
 
       [msg1,mtd,m,x,c,a,s,r,f,upd,frc,d,b,t,o] = getopt(v1);

       ok = isempty(msg1);

       if ok
          dok = ( ( d == depth   ) | ( isnan(d) & isnan(depth) ) );
          ok  = ( ( r == recurse ) & ...
                  ( f == follink ) & ...
                    dok            & ...
                  ( b == bytes   ) & ...
                  isequal(o,option)          );
       end

    end

    nd(ii,3) = ok;

    if ok
       break
    end

end

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,ARC_File,OPT_File,TMP_File,ini] = repeat(bcp,vin);

% Repeat the Backup using the Options from existing BCP-File

msg = '';
nl  = char(10);

ARC_File = '';
OPT_File = '';
TMP_File = '';

ini = zeros(0,4);  %  [ NDir NFiles Bytes ArchBytes ] 

ref_in = { char(10) char(13) [char(32) char(9)] };

%---------------------------------------------------------
% Load OptionFile, check with BCP-Extension

ext = { ''  optext };

msg    = cell(1,2);
msg(:) = {''};

for ee = [ 1  2 ]
    file = cat(2,bcp,ext{ee});
    if ~( exist(file,'file') == 2 )
        msg{ee} = sprintf('OptionFile "%s" doesn''t exist.',bcp);
    else
        [msg{ee},bb] = loadfile(file,1024,'char');
        if ~isempty(msg{ee})
            msg{ee} = sprintf('Error call LOADFILE(''%s'').\n%s',bcp,msg{ee});
        elseif isempty(bb)
            msg{ee} = sprintf('OptionFile "%s" is empty.',bcp);
        else
            [msg{ee},opt] = get_ref(bb,'>O',ref_in{:});
            if ~( isempty(msg{ee}) & ~isempty(opt) );
                 msg{ee} = sprintf('Can''t get Options from OptionFile "%s".',bcp);
            end
        end
    end
    if isempty(msg{ee})
       break
    end
end

if ~isempty(msg{ee})
    msg = msg{1};
    return
end

msg = '';

%---------------------------------------------------------
% Get Options

try
   opt = eval([ '{' opt{end,2} '}' ]); 
catch
   msg = sprintf('Can''t evaluate Options from OptionFile "%s".',bcp);
   return
end

msg = checkbase(opt);
if ~isempty(msg)
    msg = sprintf('Invalid Options from OptionFile "%s".\n%s',bcp,msg);
    return
end

if isempty(opt)
   opt = {''};
end

%-------------------------------------------------
% Check for SourceDirectory in Options

v = opt{1};

if isempty(v)
   [m,v] = get_ref(bb,'>P',ref_in{:});
   ok = isempty(m) & ~isempty(v);
   if ok
      v  = v{end,2};
      ok = ( exist(v,'dir') == 7 );
      if ok
         opt{1} = v;
      end
   end
   if ~ok
       opt = opt{2:end};
   end
elseif size(v,2) >= 3
   if strcmp( v([1 end]) , '<>' )
      msg = sprintf('Recursive Repetition in Options of OptionFile "%s".',bcp);
      return
   end
end

%-------------------------------------------------
% Call BACKUP

try
   [msg,ARC_File,OPT_File,TMP_File,ini] = backup(opt{:},vin{:});
catch
    msg = lasterr;
end

if ~isempty(msg)
    msg = sprintf('Error call BACKUP.\n%s',msg);
end

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


function f = relpath(f,d)

% Makes relative Path

isdos = strcmp(upper(computer),'PCWIN');

fs = filesep;

[p,n,e] = fileparts(f);

p = rmblank(p,2,fs);
d = rmblank(d,2,fs);

is = ( p == fs );
ns = sum(is)+1;
is = cat(2,0,find(is),size(p,2)+1);
pp = cell(1,ns);
for ii = 1 : ns
    pp{ii} = p((is(ii)+1):(is(ii+1)-1));
end
 
is = ( d == fs );
ns = sum(is)+1;
is = cat(2,0,find(is),size(d,2)+1);
dd = cell(1,ns);
for ii = 1 : ns
    dd{ii} = d((is(ii)+1):(is(ii+1)-1));
end

np = size(pp,2);
nd = size(dd,2);

p = cell(1,0);

for ii = nd : -1 : 1
    ok = ( ii <= np );
    if ok
       if isdos
          ok = strcmp(lower(pp{ii}),lower(dd{ii}));
       else
          ok = strcmp(pp{ii},dd{ii});
       end
    end

    if ~ok
        p = cat(2,p,{'..'});
    else
        break
    end
end

p = cat(2,p,pp(ii+1:np));

f = cat(2,fullfile(p{:},n),e);

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = checkname(name,chk)

% Checks for valid Characters

if nargin < 2
   chk = 'file';
end

  %---------------------------------------------------
  % Check for valid Characters
  % 0 .. 9  |  A .. Z   |  a .. z  |  .  |   _  |  FileSeparator

  name = double(name); 
    
    fs = double(filesep);

  ok = all(  ( (  48 <= name  &  name <= 57  )  | ...
               (  65 <= name  &  name <= 90  )  | ...
               (  97 <= name  &  name <= 122 )  | ...
                 name == 46   |  name == 95     |   name == 126   | ...
               ( strcmp(chk,'dir') & ( name == fs  | ...
                  ( name == 58  & strcmp(upper(computer),'PCWIN') ) ) )   ) );


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [n,p,ok] = subname(p,def);

% SUBNAME   Split Directory- or FileName in Root and Sub  
%
% [ SubName , RootName , UseDefault ] = SUBNAME( DirName , DefaultSubName )
%
%
% see also: FILEPARTS
%

  Nin = nargin;

  n = '';

  ok = 1;

  if Nin < 1
     p = '';
  end

  if Nin < 2
     def = '';
  end

  if ~( ischar(p)  &  ( prod(size(p)) == size(p,2) ) &  ...
        ischar(def)  &  ( prod(size(def)) == size(def,2) )         );
     error('Inputs must be Strings.');
  end

  fs = filesep;

  n  = def;

  if isempty(p)
     return
  end

  %-------------------------------------
  % Remove Duplicate FileSep

  ii = find( double(fs) == double(p) );

  if ~isempty(ii)

       jj = find( diff(ii) == 1 );

       p(ii(jj+1)) = [];

  end

  if strcmp(p,fs)
     return
  end

  p = p( 1 : ( end - strcmp(p(end),fs) ) );

  %--------------------------------------------
  % Split
 
  ok = 0;

  ii = find( double(fs) == double(p) );

  if isempty(ii)
     ii = 0;
  else
     ii = max(ii);
  end

  n = p(ii+1:end);
  p = p(1:ii);


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [c,l,t] = dircont(pfad,varargin)

% DIRCONT  Returns Contents of Directory
%
% [ C , LISTE , TEXT ] = DIRCONT( Pfad , Mode , WildCard )
%
% gives Contents of Directory "Pfad"
%
%  C    = { IsDir  Date  Bytes  Name  DateNum };
%
%  LISTE = CellString(C)
%
%  TEXT  = STRHCAT( LISTE,'',1);
%
%  Mode == {'n'} , sort by Name            of c(:,4), default
%          {'a'}           Alphabetical
%          {'c'}           Character
%          {'l'}           Letter
%          {'s'}           String
%          {'b'}           Buchstaben
%          {'w'}           Woerter
%
%  Mode == {'t'} , sort by Time == DateNum of c{:,5}
%          {'d'}           DateTime
%          {'z'}           Zeit
%
% If WildCard ends with FileSep, only Directory will returned.
%
%
% !!! SPECIAL EDITION FOR BACKUP, CHECK LOWERCASE and UPPERCASE of WILDCARD !!!
% !!!    CHECK for Files without Extension with '*.-' !!!
%

Nin  = nargin;
Nout = nargout;


Nout = nargout;

fs = filesep;

c = cell(0,5);  % { IsDir Date  Bytes  Name DateNum}
l = cell(0,1);
t = '';


if nargin < 1
  pfad = cd;
else
  if ~( ischar(pfad)  &  ( prod(size(pfad)) == size(pfad,2) )  )
    error('Input must be a String.');
  end
  if isempty(pfad)
     pfad = fs;
  end
end
 
mode = 'a';
wc   = '';


for ii = 1 : Nin-1

   val = varargin{ii};

   if iscellstr(val)

      if ~isempty(val{1})
        mode = lower(val{1}(1));
      end

   elseif ( ischar(val)  &  ~isempty(val) & ...
            ( prod(size(val)) == size(val,2) )     )

      wc = val;

   end

end

get_dir = 0;
if ~isempty(wc)
   get_dir = strcmp(wc(end),fs);
   if get_dir
      wc = wc(1:end-1);
   end
end

% Format for Bytes
form = '%8.0f';  

% Seperator: [ sep0 Date sep1 Bytes sep2 Name ];

sep0 = ' ';
sep1 = '  ';
sep2 = '  ';

%------------------------------------------------
% Start

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Special

chk0 = 0;   % Check for Files without Extension
chk2 = 0;   % Ccheck for Lower- and UpperCase Extension

if isempty(wc)

   d = dir(pfad);

else

   pfad = cat( 2 , pfad , fs(1:(end*(~strcmp(pfad(end),fs)))) );

   chk0 = strcmp(wc,'*.-');

   if chk0
      wc = '*';
   end

   lc = lower(wc);
   uc = upper(wc);

   chk2 = ( isequal(wc,lc) | isequal(wc,uc) );

   chk2 = ( chk2 & (~chk0) & ~isequal(lc,uc) );

   chk2 = ( chk2 & ~strcmp(upper(computer),'PCWIN') );

   if chk2
      d = cat( 1 , dir(cat(2,pfad,lc)) , dir(cat(2,pfad,uc)) );
   else
      d = dir(cat(2,pfad,wc));
   end

end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


N = size(d,1);


if N == 0
  return
end


sep0 = ' ';
sep1 = '   ';
sep2 = '   ';

l = cell(N,1);
c = cell(N,5);

for ii = 1 : N

    if chk0 & ~d(ii).isdir
       n = d(ii).name(2:end);  % '.' at begin is allowed
       if ~isempty(n)
          if any( double(n) == double('.') )
             d(ii).name = '.';   % Will removed later from List !!!
          end
       end
    end
 
    c{ii,1} = d(ii).isdir * ( 1 - 2 * strcmp(d(ii).name,'..') );
    c{ii,2} = d(ii).date;
    c{ii,3} = d(ii).bytes;
    c{ii,4} = d(ii).name;
 
    if isempty(d(ii).date)
       c{ii,5} = NaN;
    else
       c{ii,5} = datenum(ger2eng(d(ii).date));  
    end

end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if 0 %%% strcmp(upper(computer),'PCWIN')

   c(:,4) = lower( c(:,4) );   % Make LowerCase

   % Remove Duplicate Names

   [cc,si]= sort(c(:,4));
    cc    = char(cc);
    jj    = ( sum( ( diff(double(cc),1,1) == 0 ) , 2 ) == size(cc,2) );
    if any(jj)
       jj = find(jj) + 1;
       c(si(jj),:) = [];
    end

end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


nd = size(char(c{:,2}),2);

for ii = 1 : N*( Nout >= 2 )

     fs_end = size(fs,2) * ( c{ii,1}                      &  ...
                            ~strcmp(c{ii,4}(end) ,  fs  ) &  ...
                            ~strcmp(c{ii,4}      , '..' )        );

    bytes  = sprintf( form , c{ii,3});
    if d(ii).bytes == 0
       bytes = char( 32 * ones( 1 , size(bytes,2) ) );
    end

    if isempty(c{ii,2})
       date = char( 32 * ones(1,nd) );
    else
       date = cat( 2 , c{ii,2} , char( 32 * ones(1,nd-size(c{ii,2},2)) ) );
    end

    l{ii}  = cat( 2 , sep0 , date                ,  ...
                      sep1 , bytes               , ...
                      sep2 , c{ii,4} , fs(1:fs_end)        );

 
end


c0 = cell(0,5);
l0 = cell(0,1);


if get_dir

   bad = ( strcmp( c(:,4) , '.' )  |  ~cat(1,c{:,1}) );

else

   jj = find( strcmp( c(:,4) , '..' ) );

   if ~isempty(jj);

     c0 = c(jj,:);
     l0 = l(jj);

     c(jj,:) = [];
     l(jj)   = [];

   end

   bad = strcmp( c(:,4) , '.' );
   
end


bad = find(bad);

c(bad,:) = [];
l(bad,:) = [];

% Sort alphabetical

is_dir = cat(1,c{:,1});

ii = find( is_dir);
jj = find(~is_dir);

switch mode
  case { 't' 'd' 'z' }

    % via DateNum

    [h,s_ii] = sort( (-1)*cat(1,c{ii,5}) );
    [h,s_jj] = sort( (-1)*cat(1,c{jj,5}) );

  case { 'n' 'a' 'c' 'l' 's' 'b' 'w' }

    % via Name

    [h,s_ii] = sort( c(ii,4) );
    [h,s_jj] = sort( c(jj,4) );

  otherwise

       s_ii  = ( 1 : prod(size(ii)) );
       s_jj  = ( 1 : prod(size(jj)) );

end

c = cat( 1 , c0 , c(ii(s_ii),:) , c(jj(s_jj),:) );
l = cat( 1 , l0 , l(ii(s_ii),:) , l(jj(s_jj),:) );


if Nout == 3
   t = strhcat(l,'',1);
end

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function str = ger2eng(str)

% GER2ENG   Converts German MonthNotation to English
%
% EnglishString = GER2ENG( GermanString );
%
%  mär --> mar
%  mai --> may
%  okt --> oct
%  dez --> dec
%
% In other case you get problems to use DATEVEC, DATENUM
%

if nargin < 0
   str = '';
   return
end

if ~( iscellstr(str) | ischar(str) );
   error('Input must be a CharArray or CellStringArray.');
end


is_char = ischar(str);

if is_char
   str = cellstr(str);
end

rep = { ...
  'mär'  'mar'
  'mai'  'may'
  'okt'  'oct'
  'dez'  'dec'  };


for ii = 1 : prod(size(str))

    s = str{ii};

    if ~isempty(s)  

        s = lower(s); 

        for jj = 1 : size(rep,1)

            kk = findstr( s , rep{jj,1} );

            if ~isempty(kk)

               n = size(rep{jj,1},2);

               for ll = kk(:)'

                   ind = ll + (1:n) - 1;     

                    nn = find( ~( double( s(ind) ) == double( str{ii}(ind) ) ) );

                    str{ii}(ind) = rep{jj,2};

                    str{ii}(ind(nn)) = upper(str{ii}(ind(nn)));

               end
               % ll == kk(findstr)
            end
            % ~isempty(findstr) 
        end
        % jj = rep
    end
    % ~isempty(str{ii})
end
%ii

if is_char
   str = char(str);
end


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
%   default: N = NStrings;  NewLine = char(10);
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
 n = [];
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

if isempty(n)
   n = size(str,1) + 1;
end

str( : , 2 ) = { del };

str(n:n:size(str,1),2) = { nl };

str(    size(str,1),2) = { '' };

str = permute( str , [ 2  1 ] );

str = cat(2,str{:});


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,bb] = loadfile(file,varargin);

% LOADFILE  Load binary Data from File, using FREAD
%
% [ Msg , V ] = LOADFILE( FileName , MaxSize , Precision )
%
%   MaxSize     MaximumFileSize [Bytes],   default: 2097152 == 2MBytes
%   Precision   Precision for using FREAD, default: 'char'
%
%  For more Informations  type: >> help fread
%
%  In case of Precision 'char', a CharacterString will returned,
%   if all Bytes are valid Characters for conversion, if some
%   Characters are invalid, Msg is not empty.
%

msg = '';
bb  = [];

msg0 = 'LOADFILE: ';

%------------------------------------------------
% Defaults

m = 2*2^20; % 2 MBytes;  % MaximumFileSize [Bytes]
p = 'char';

Nin = nargin;

%************************************************
% Check File

if Nin == 0
   msg = [msg0 'Input File is missing.'];
   return
end

if isempty(file)
   return
end

if ~chkstr(file,1)
   msg = [msg0 'Input File must be a String.'];
   return
end

%------------------------------------------------
% Get MaxSize and Precision

for ii = 1 : Nin-1

    v = varargin{ii};

    if chkstr(v,1)

       p = v;

    elseif ( isnumeric(v)  &  ( prod(size(v)) == 1 ) )

       if ( isfinite(v)  &  ( v >= 0 ) )
          m = v;
       end

    end

end

%************************************************
% Open and Read File

  fid = fopen(file,'r');

  if fid == -1  
     msg = [ msg0  'Can''t open File.' ];
     return
  end

 %----------------------------------------------
 % Check Size of File

  d = dir(file);

  if isempty(d)
   d = dir( which(file) );  % which(file) gives the full Name
                            %  [ PathName FileName ]
  end

  if d.bytes > m
    msg = [ msg0 'File too large, Limit = '  ...
            sprintf('%.0f Bytes',m) '.' ];
    fclose(fid);
    return
  end


 %----------------------------------------------

  try
     bb = fread(fid,p);
  catch
     msg = [ msg0 'Error call FREAD.' char(10) lasterr ];
  end

  fclose(fid);

 %----------------------------------------------
 % Precision: 'char'  ==>  Transform to String

  if isempty(msg) & isequal( p , 'char' )
    
    % Check Characters

    if all( ( bb ==  9 ) |  ...
            ( bb == 10 ) |  ...
            ( bb == 13 ) |  ...
            (  28 <= bb  &   bb <= 126 ) | ...
            ( 160 <= bb  &   bb <= 255 )        );

       bb = char(bb(:)');

    else

       msg = [ msg0 'Invalid Characters in File.' ];

    end


  end

%***********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = chkstr(str,opt)


% CHKSTR  Checks Input for String
%
%  ok = chkstr(str,Option)
%
%  Option ~= 0 ==> only true for nonempty Strings
%
%   default:   Option == 0  (true for empty Strings)
%

 
if nargin < 2
   opt = 0;
end

ok = ( strcmp( class(str) , 'char' )      & ...
       ( prod(size(str)) == size(str,2) ) & ...
       ( isequal(opt,0) | ~isempty(str) )         );


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,out] = get_ref(str,s1,s2,c1,c2)

% GET_REF  Returns a Reference from a String
%
% [Msg,Reference] = GET_REF( String , S1 , S2 , C1 , C2 );
%
%  S1:  StartMarker  Characters which marks the Begin of Reference
%  S2:    EndMarker  Characters which marks the End   of Reference
%
%  C1:  Characters in Reference to remove,
%        default: [ CR NL TAB ]
%  C2:  Characters closed to C1 to remove,
%        default: [ Space '-' ]
%
%  Reference = { FullString  ReferenceString }  2-Column Cell-Array
%
%

Msg = '';

out = cell(0,2);


nl = char(10);

Msg0 = 'GET_REF: ';

nm0 = size(Msg0,2);

nl0 = char([ 10 32*ones(1,nm0+0) ]);

%---------------------------------------------------------

if ~( ischar(str) &  ( prod(size(str)) == size(str,2) ) )
  Msg = '1. Input must be a String.';
end
 
if ~( ischar(s1) &  ( prod(size(s1)) == size(s1,2) ) )
  Msg = [ Msg nl0(1:(end*(~isempty(Msg)))) ...
           'S2 must be a String.' ];
end
 
if ~( ischar(s2) &  ( prod(size(s2)) == size(s2,2) ) )
  Msg = [ Msg nl0(1:(end*(~isempty(Msg)))) ...
           'S2 must be a String.' ];
end

%---------------------------------------------------------

if nargin < 4

  c1 = [ 13  10  9 ];

else

  if ischar(c1)
    c1 = double(c1);
  end
  ok = isnumeric(c1);
  if ok & ~isempty(c1)
    c1 = c1(:)';
    ok = all( ( mod(c1,1) == 0 )  & ...
              ( c1 >= 0 ) & isfinite(c1)  );
  end
  if ~ok
      Msg = [ Msg nl0(1:(end*(~isempty(Msg)))) ...
              'Input C1 must be a String or ASCII-Codes.'];
  end

end   

%---------------------------------------------------------

if nargin < 5

  c2 = [ double('-')  32  ];

else

  if ischar(c2)
    c2 = double(c2);
  end
  ok = isnumeric(c2);
  if ok & ~isempty(c2)
    c2 = c2(:)';
    ok = all( ( mod(c2,1) == 0 )  & ...
              ( c2 >= 0 ) & isfinite(c2)  );
  end
  if ~ok
      Msg = [ Msg nl0(1:(end*(~isempty(Msg)))) ...
              'Input C1 must be a String or ASCII-Codes.'];
  end

end   

%---------------------------------------------------------

if ~isempty(Msg)
  Msg = [ Msg0  Msg ];
end

if isempty(str)  |  isempty(s1)  |  isempty(s2)  |  ...
   ~isempty(Msg)
 
  return

end


%**********************************************************

n = size(str,2);


%--------------------------------
% Start of Reference at End of s1

i1 = findstr(str,s1);

if isempty(i1)
  return
end

i1 = i1 + size(s1,2) - 1;
i1( find( i1 > n ) ) = [];

if isempty(i1)
  return
end
 
%--------------------------------
% End of Reference at Begin of s2

i2 = findstr(str,s2);

if isempty(i2)
  return
end

i2( find( i2 < i1(1) ) ) = [];

if isempty(i2)
  return
end

%----------------------------

ok = zeros(1,n);  % Start/End
ii = zeros(1,n);  % Index

ok(i1) = 1;  % Start of Reference
ok(i2) = 2;  % End   of Reference

ii(i1) = i1;
ii(i2) = i2;

jj = find(ok);

ok = ok(jj);
ii = ii(jj);


%----------------------------
% Following Start-End

ok = find( abs(diff(ok)) == 1 );

if isempty(ok)
   return
end

%----------------------------
% Start's

ok = ok(1:2:end);

n = prod(size(ok));

out      = cell(n,2);
out(:,1) = { [ s1  s2 ] };
out(:,2) = { '' };

for jj = 1 : n

  if ( ii(ok(jj)+1) - ii(ok(jj)) ) > 1
  % Reference not empty

    % String in Reference
 
    sr = str( ii(ok(jj))+1 : ii(ok(jj)+1)-1 );


    % Full String

    out{jj,1} = [ s1  sr  s2 ];


    % Remove bad Characters from Reference

    sr = rmblank( sr , 2 );

    b1 = zeros( 1 , size(sr,2) );
    for r1 = c1
      b1 = ( b1  |  ( double(sr) == r1 ) );
    end

    b2 = zeros( 1 , size(sr,2) );
    for r2 = c2
      b2 = ( b2  |  ( double(sr) == r2 ) );
    end
  

    i0 = find( diff(cat(2,0,( b1 | b2 )  )) ==  1 );  % Start of Group
    i1 = find( diff(cat(2,  ( b1 | b2 ),0)) == -1 );  % End   of Group
 
     b = 2 * b1 + 1 * b2;
    cb = cumsum(b,2);

    lg = i1 - i0 + 1;              % Length of Group
    sg = cb(i1) - cb(i0) + b(i0);  % Sum    of Group         

    ib = find( sg > lg );          % Groups incl. b1

    lg = lg(ib);
    i0 = i0(ib);

    % Build IndexVector to remove

     b = grp2ind(i0,lg);

     sr(b) = [];
 
     out{jj,2} = sr;

  end
  % Reference not empty

end
% jj

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,typ,link,name] = islink(name);

% ISLINK  Checks if argument is a symbolic Link
%
% [ IsLink , Type , Value , Target ] = ISLINK( Name )
%
%  IsLink:  1 if Name is a symbolic link by UNIX: test -L,
%           0 if Name is not a symboliclink or on non-UNIX systems
%          [] if Name is empty
%
%  Type     2 if Name or the target of the link is a file,
%           7 if Name or the taget of the link is a directory
%           0 if Name is not a File or Directory
%    
%  Value:  Value of symbolic Link (Derefer) by UNIX: readlink -m
%
%  Target: Target of symbolic Link by UNIX: ls -ld; or origin File or Directory
%
%
% uses UNIX commands: test, readlink, ls
%

Nout = nargout;

ok     = NaN;
typ    = 0;
link   = '';

if ischar(name) & isempty(name)
   ok = []; 
   return
elseif ~( ischar(link) &  ( size(typ,2) == prod(size(typ)) ) )
   error('Input must be a String.');
end

%------------------------------------------------------
% Check if Directory or File exist

typ = 2 * ( exist(name,'file') == 2 ) + ...
      7 * ( exist(name,'dir')  == 7 );

ok = ~( typ == 0 );

if ok
   link = name;
end

ok = ( ok & isunix );

if ~ok
    return
end


%------------------------------------------------------
% Remove FileSeparator from End in case of Directory !!!

name = name( 1 : end-( strcmp( link(end) , filesep ) & ( typ == 7 ) ) );

%------------------------------------------------------
% Check for Link, use UNIX: test

[s,w] = unix(sprintf('test -L "%s"',name));

ok = ( s == 0 );

if ~ok | ( Nout < 3 )
    return
end

%------------------------------------------------------
% Derefer Link, use UNIX: readlink

[s,link] = unix(sprintf('readlink -n -m "%s"',name));

if ~( s == 0 )
    link = '';
end

if ( Nout < 4 )
    return
end

%------------------------------------------------------
% Get Target of Link, use UNIX: ls

pre = '->';

[s,l] = unix(sprintf('ls -ld "%s"',name));

ii = findstr( l , pre );

if isempty(ii)
   return
end

ii = max(ii) + size(pre,2);

name = rmblank(l(ii:end),2);

