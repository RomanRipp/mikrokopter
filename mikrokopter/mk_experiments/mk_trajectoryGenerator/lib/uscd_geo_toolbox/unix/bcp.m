function bcp(mode)

% BCP   Wizard for BACKUP
%
% BCP(ArchiveMode)
%
% ArchiveMode = 'scr' | 'img' | 'dat' | 'arc' | 'matlab' | 'home' | 'all'
%
% ArchiveMode can be a CellStringArray of Modes
%
%  'scr' includes 'matlab'
%
%  'home'   ==  '$HOME'        , ':scr' , '-nl' 
%  'matlab' ==  '$HOME/matlab' , ':scr' , '-nl' 
%
%  'all' == { 'scr' 'img' 'dat' 'arc' 'matlab' 'home' };
%
 
%*******************************************************************
% Settings

Source = '/d1/project/matlab';
Destin = '/d0/backup/project/';   % FileSep at END !!!

ArchiveName  = 'project_matlab';

HomeSource   = getenv('HOME');

MatlabSource = fullfile(HomeSource,'matlab'); % Special for ArchiveMode 'src' (Script)

%*******************************************************************


if nargin == 0
   mode = 'scr';
end

if ischar(mode)
   mode = cellstr(mode);
   mode = mode(:)';
end

if any(strcmp(mode,'all'))
   mode = { 'scr' 'img' 'dat' 'arc' 'matlab' 'home' };
elseif any(strcmp(mode,'scr')) & ~any(strcmp(mode,'matlab'))
   mode = cat(2,mode,{'matlab'});
end

option = { '-mbz'  '*'  '-nl' '|database'  '-M1' };

for mm = mode

   switch mm{1}

     case 'home'

       opt = ':scr';
    
       backup(HomeSource,Destin,option{:},opt);

     case 'matlab'

       opt = ':scr';

       backup(MatlabSource,Destin,option{:},opt);

     otherwise

       opt = cat(2,':',mm{1});

       backup(Source,cat(2,Destin,ArchiveName),option{:},opt);

   end

end

