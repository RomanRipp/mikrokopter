% File- and DirectoryHandling # Contents of "toolbox/diskfun"
%   
% Create backups:              BACKUP
% Search for Files:            WHICHFILE
% Update M-Files with:         CP2FIND, UPDATE
% Find Strings in Files:       MGREP, BGREP
% Get contents of MAT-Files:   WHOSFILE, LOOK_MAT
% Translate Strings in Files:  TRANSLATE
% Create "Contents.m"-Files:   MKCONT
%   
%-------------------------------------------------------------------------------------
%   
%   backup        - Creates a Backup-Archive of an Directory from selected Extensions
%   bgrep         - Search in binary Files for specific string
%   bkmstruct     - Converts NetscapeBookmarkFile to Structure
%   check4double  - Checks Directory recursively for duplicate Files
%   cp2find       - Copy Files recursive to Destination
%   cp_list       - Copy Files from ListFile using CP2FIND
%   dircont       - Returns Contents of Directory
%   dirhist       - Returns Directory-History
%   dirinfo       - Returns Info about Contents of Directory
%   dirlist       - Returns String-List from Output of DIRINFO
%   enc_dos       - Display's bad Characters from DOS-ASCII-Files
%   filelist      - Creates FileList
%   islink        - Checks if argument is a symbolic Link
%   loadfile      - Load binary Data from File, using FREAD
%   look_mat      - Display Informations of a MAT-File
%   lowername     - Renames Files to LowerCase-Names
%   mgrep         - Search in Files for specific string
%   mkcont        - Creates "Contents.m"- File of Directory
%   mkpcode       - Creates pre-parsed pseudo-code file
%   numfiles      - Numerate Files by ascending Order of Date
%   relpath       - Translates absolute Links into relative
%   remdir        - Removes Directory
%   translate     - Translates String in Files to another
%   update        - Update Files in Directory  by Suffix
%   whichfile     - Find Files in Directory or Matlab's SearchPath by Suffix or exact
%   whosfile      - Returns Contents of MAT-File as CellString
%   
