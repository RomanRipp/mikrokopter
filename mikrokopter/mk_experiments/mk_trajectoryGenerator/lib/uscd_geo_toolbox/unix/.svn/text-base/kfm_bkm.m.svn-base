function kfm_bkm(c,pfad);

% KFM_BKM  creates KFM-Bookmarks from BookmarkStructure
%
% KFM_BKM( S , Directory )
%
% default: Directory = '$HOME/.kde/share/apps/kfm/bookmarks'
%
% S is a BookmarkStructure, returned by BKMSTRUCT
%
% S = isdir:   0  |  1          % true for Folder
%      link:  URL | SubStruct
%      type:  LinkType { 'folder' | 'file' | 'http' | 'ftp' }
%      name:  LinkName
% 
% For Folders SubDirectories will created, for URLs KDELNK-Files.
%

Nin = nargin;

if Nin < 1
   error('Structure is missing.');
end

if Nin < 2
   pfad = fullfile(getenv('HOME'),'.kde/share/apps/kfm/bookmarks');
end

fs = filesep;

pfad = cat( 2 , pfad , fs(1:(end*(~strcmp(pfad(end),fs)))) );

%*****************************************************************
% Settings

icon = struct( 'ftp'  , { 'mini-run.xpm' } , ...
               'http' , { 'mini-gv.xpm'  } , ...
               'file' , { 'mini-doc.xpm' } , ...
               'none' , { 'html.xpm'     }       );

fi = fieldnames(icon);

def = { ...
'# KDE Config File'
'[KDE Desktop Entry]'
'URL=@LINK'
'Icon=@ICON'
'MiniIcon=@ICON'
'Name=@NAME'
'Comment=@NAME'
'Type=Link'       };

def = sprintf('%s\n',def{:});

ext = '.kdelnk';

%*****************************************************************

nc = prod(size(c));

nm = cell(nc,1);

for ii = 1 : nc

    name = c(ii).name;

    if ~strcmp(c(ii).type,'file')
        n = min(size(name,2),48);
        name = name(1:n);
    end

    name = strrep(name,'/','%2f');

    if ii > 1
       nm0 = name;
         z = 0;
       while any(strcmp(name,nm(1:ii-1)))
             name = sprintf('%s_%.0f',nm0,z);
             z = z+1;
       end
    end

    nm{ii} = name;

    %--------------------------------------------------------
    if c(ii).isdir

       command = sprintf('mkdir "%s%s"',pfad,name);

       fprintf(1,'Create Folder: "%s%s"',pfad,name);

       [s,w] = unix(command);

       if ~( s == 0 )
           fprintf(1,' %s\n','error');
       else
           fprintf(1,'\n');
           if isstruct(c(ii).link)
              kfm_bkm(c(ii).link,[pfad name]);
           end
       end


    %--------------------------------------------------------
    else

       if any(strcmp(c(ii).type,fi))
          ic = fi{find(strcmp(c(ii).type,fi))};
       else
          ic = 'none';
       end

       str = strrep(def,'@LINK',c(ii).link);
       str = strrep(str,'@NAME',c(ii).name);
       str = strrep(str,'@ICON',getfield(icon,ic));

       lnk = [ pfad  name ext ];

       fprintf(1,'Create   Link: "%s"',lnk);

       fid = fopen(lnk,'wt+');
       if ( fid == -1 )
           fprintf(1,' %s\n','error');
       else
           fprintf(1,'\n');
           fprintf(fid,'%s',str);
           fclose(fid);
       end

    end

end