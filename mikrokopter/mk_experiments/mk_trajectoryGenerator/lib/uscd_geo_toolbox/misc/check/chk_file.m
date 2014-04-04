function [Msg,name,file] = chk_name(name,varargin);

% CHK_NAME   Checks FileName of valid Characters, Convention and Directory
%
% [ Msg , Name , FullName ] = CHK_NAME( Name , '|path' , '.ext' , '-pre' , Check , auto )
%
%    Check = 'new'    Warning, if File or Directory exists 
%            'none'   no check for existing File/Directory
%            'file'   Warning, if File not exist 
%            'dir'    Warning, if Directory not exist 
%
% 

Msg = '';

nl = char(10);

file = '';

%-------------------------------------------------------------

if nargin < 2
  pfad = '';
end

if nargin < 3
  ext = '';
end

if nargin < 4
  pre = '';
end

if nargin < 5
 chk = 'none';
end

%-------------------------------------------------------------

if ~( ischar(pfad)  &  ( prod(size(pfad)) == size(pfad,2) )  )
   Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
           'Input Pfad must be a String.'      ];
end

if ~( ischar(name)  &  ( prod(size(name)) == size(name,2) )  )
   Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
           'Input Name must be a String.'      ];
end

if ~( ischar(ext)  &  ( prod(size(ext)) == size(ext,2) )  )
   Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
           'Input Extension must be a String.'      ];
end

if ~( ischar(pre)  &  ( prod(size(pre)) == size(pre,2) )  )
   Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
           'Input Prefix must be a String.'      ];
end

if ~( ischar(chk)  &  ( prod(size(chk)) == size(chk,2) )  )
   Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
           'Input Check must be a String.'      ];
end

if ~isempty(Msg)
   return
end

%-------------------------------------------------------------

name = rmblank(name,2);

  fs = filesep;

if ~isempty(pfad)
    pfad = pfad( 1 : (end-strcmp(pfad(end),fs)) );
end


if ~isempty(name)

  jj = find( name == fs );

  if ~isempty(jj)

      jj = max(jj);

      if isempty(pfad)
         pfad = name( 1 : (jj-1) );
      else
         if ~strcmp(pfad,name(1:(jj-1)));
             Msg = sprintf('Directory "%s" is not conform with "%s".',pfad,name(1:(jj-1)));
             return
         end
      end
      
      name = name( (jj+1) : end );

   end

end


if ~isempty(name)

  %---------------------------------------------------
  % Check for valid Characters
  % 0 .. 9  |  A .. Z   |  a .. z  |  .  |   _  |  FileSeparator

  name = double(name); 
    
    fs = double(fs);

  ii = find( ~( (  48 <= name  &  name <= 57  )  | ...
                (  65 <= name  &  name <= 90  )  | ...
                (  97 <= name  &  name <= 122 )  | ...
                  name == 46   |  name == 95     |  name == 126  | ...
                ( name == fs   &  strcmp(chk,'dir') )   ) );

  name = char(name);

  if ~isempty(ii)
     Msg = 'Invalid Character''s in Name.';
     return
  end

  %---------------------------------------------------
  % Append pre

  if ~isempty(pre)

    s2 = size(name,2);

    sp = size(pre ,2);

    if ~strcmp(lower(pre),lower(name(1:min([sp s2]))))
       name = cat( 2 , pre , name );
    end

  end

  %---------------------------------------------------
  % Append ext

  if ~isempty(ext)


    f  = '.';

    ext = cat( 2 , f(1:(end*(~strcmp(ext(1),f)))) , ext );

    s2 = size(name,2);

    se = size(ext,2);

    if ~strcmp(lower(ext),lower(name(s2-min([se s2])+1:s2)))
       name = cat( 2 , name , ext( (1+strcmp(name(s2),f)) : se ) );
    end
    
  end
  
end 


%---------------------------------------------------
% Append pfad

if ~isempty(pfad)

   pfad = cat( 2 , pfad , fs );

   if isempty(name)
      file = pfad;
   else
      file = cat( 2 , pfad , name );
   end

else

   file = name;

end


%-----------------------------------------------------
% Check

isfile = ( exist(file,'file') == 2 );
isdir  = ( exist(file,'dir')  == 7 );

chk = lower(chk);

switch chk

   %---------------------------------------------------
   case 'new'

     if isfile
        Msg = [ 'File allready exist: ' file ];
     elseif isdir
        Msg = [ 'Directory allready exist: ' file ];
     end
   
   %---------------------------------------------------
   case 'file'

     if ~isfile
        Msg = [ 'File does not exist: ' file ];
     end

     if isdir
        Msg = [ 'Directory exist: ' file ];
     end

   %---------------------------------------------------
   case 'dir' 

     if ~isdir
        Msg = [ 'Directory does not exist: ' file ];
     end

     if isfile
        Msg = [ 'File exist: ' file ];
     end

end
