function out = mgrep(varargin)

% MGREP  Search in Files for specific string
% 
% mgrep(String,Pfad,'-r','-l','-v','-d','-n','-K#','.ext1','.ext2', ... )
% mgrep(String,File,'-r','-l','-v','-d','-n','-K#','.ext1','.ext2', ... )
%
%  '-r'  ==  recurse
%  '-l'  ==  follow DirectoryLinks
%  '-v'  ==  verbose, list Directories
%  '-d'  ==  no Display
%  '-n'  ==  no IntroDisplay
%  '-K#' ==  Max. Size of Files in KBytes 
%  '-M#' ==  Max. Size of Files in MBytes 
%
%  defaults:   no recurse 
%              no follow DirectoryLinks
%              Display
%              IntroDisplay
%             '.ext' == '.m'
%             '-K#'  == '-K1024'  % 1-MBytes MaxSize
%             '-M#'  == '-M1'     % 1-MBytes MaxSize
%
%  note:  give '.'  to grep all Extensions
%    
%  out = { FileName  Date  Bytes  Number LineNumber  String }
%


out = cell(0,6);

%****************************************************************************
% Get and Check Inputs

try
   [msg,str,Pfad,ext,recurse,follink,verbose,nodispl,nointro,kb] = checkin(varargin{:});
catch
   msg = sprintf('Error call MGREP:CHECKIN\n%s',lasterr);
end
 
if ~isempty(msg)
   error(msg);
end

%****************************************************************************

fs = filesep;
nl = char(10);
wc = '*';

 
if ~nointro & ~nodispl

  txt = 'MGREP: ';

  if ~recurse
    txt = [ txt 'no ' ];
  end

    txt = [ txt 'recurse, extension: ' ];
  
  if isempty(ext{1})
    txt =  [ txt '''*'', ' ];
  else
    txt = [ txt  ''''  strhcat(ext,''', ''')  ''', ' ];
  end   

    txt = [ txt ' MaxSize: ' sprintf('%g ',kb) 'KBytes' ];

    fprintf(1,'\n%s\n\n',txt);

end

verbose = ( verbose & ~nodispl );

if verbose
   fprintf(1,'%s',Pfad);
end

d = dir(Pfad);

if isempty(d)
   if ( exist(Pfad,'file') == 2 )
      Pfad = which(Pfad);
         d = dir(Pfad);
   end
end

if isempty(d)
   if verbose
      fprintf(1,' ... can''t read File or Directory\n');
   end
   return
end

%---------------------------------------------------------
% Check for File

IsFile = ( ( size(d,1) == 1 )  &  ( d(1).isdir == 0 ) );

if IsFile

   if isempty(fileparts(d(1).name))
      d(1).name = fullfile(fileparts(Pfad),d(1).name);
   end

   ext = { char(ones(1,0)) };
   wc  = ext{1};

else

    Pfad = cat( 2 , Pfad , fs(1:(end*(~strcmp(Pfad(end),fs)))) );
  
    is_dir = find(cat(1,d.isdir));

    dd = { d(is_dir).name }; % Store Diretories here

end

if verbose
   fprintf(1,'\n')
end

ns = size(str,2);

bytes = kb * 1024;

for ee = ext(:)'

  if ~IsFile
      d = dir(cat(2,Pfad,wc,ee{1}));
  end

  is_file = find(~cat(1,d.isdir));

  for ii = is_file(:)'

     if d(ii).bytes <= bytes

        name = d(ii).name;

        file = [ Pfad((1:(end*(~IsFile))))  name ];

        fid = fopen(file,'r');

        if ( fid == -1 ) & verbose

          fprintf(1,'   Error read %s\n',file);

        elseif ~( fid == -1 )

          bb = fread(fid,'char');
          fclose(fid);

          if ~any(bb==0)

            nn = cat(1,0,find(bb==10),size(bb,1)+1);
            bb = char(bb(:)');

            if ns <= size(bb,2)

              ff = findstr(bb,str);
              for jj = ff(:)'

                 i0 = sum( nn <= jj+0    );
                 i1 = sum( nn <= jj+ns-1 ) + 1;

                 ll = bb(nn(i0)+1:nn(i1)-1);

                out = cat(1,out,{ file d(ii).date d(ii).bytes jj i0  ll });

                if ~nodispl
                   fprintf(1,'%s, %5.0f: %s\n',file,i0,ll)
                end

              end
              % jj
            end
            % ns ok
 
         elseif verbose

                fprintf(1,'   Invalid Character in %s\n',file);
 
         end
          % bb ok
        end 
        % fid == -1
     end
     % bytes ok
  end
  % ii
end
% ee

if ~recurse
    return
end

opt = { '-l' '-r' '-n' sprintf('-K%.5f',kb) '-d' '-v' };


if isempty(ext{1})
   ext = {'.'};
end


n = size(opt,2);

ind = cat( 2 , ( 2-double(follink) : n-2 ) , ...
                (n-1) * ( 1 : double(nodispl) ) , ...
                 n    * ( 1 : double(verbose) )          );

opt = opt(ind);


for d = dd(:)'

    if ~any(strcmp(d{1},{ '.'  '..' }))

      pf = cat( 2 , Pfad , d{1} );

      ok = follink;

      if ~ok
         l  =  ls('-ld',['"' pf '"']);  % Without FileSep at End !!!
         ok = ~isequal(l(1),'l');
      end

      if ok
         out1 = mgrep(str,pf,opt{:},ext{:});
          out = cat(1,out,out1);
      end

    end

end

   

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,str,Pfad,ext,recurse,follink,verbose,nodispl,nointro,kb] = ...
            checkin(str,Pfad,varargin);


% CHECKIN  Check of Inputs 


msg = '';
nl  = char(10);

fs  = filesep;

Nin = nargin;


%*********************************************************
% Defaults

recurse = 0;
follink = 0;
verbose = 0;
nodispl = 0;
nointro = 0;

kb      = 1024;

ext  = cell(1,0);

ext0 = '.m';

%*********************************************************
% Get and Check Inputs

%---------------------------------------------------------
% Check Pfad and String

if Nin < 2
   Pfad = '';
end

if Nin < 1
   str = '';
   msg = 'Input String is missing.';
   return
end

if isempty(Pfad)
   Pfad = cd;
end

if isempty(Pfad)  |  any(strcmp(Pfad,{ '.'  '*' }))
   Pfad = cd;
end


if ~( ischar(Pfad)  &  ( size(Pfad,2) == prod(size(Pfad)) ) & ~isempty(Pfad) )
  msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                  'Pfad must be a String.');
end

if ~( ischar(str) &  ( size(str,2) == prod(size(str)) ) & ~isempty(str))
  msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                 'String must be a NonEmpty String.');
end



%---------------------------------------------------------
% Get Options from VarArg


if Nin <  3

  ext = { ext0 };

else

 
  VarIn = varargin;
  VarIn = VarIn(:);

  Nin = size(VarIn,1);

  msgK = '';
  msgM = '';

  for ii = 1 : Nin

    if ischar(VarIn{ii})  &  ~isempty(VarIn{ii})  &  ...
       ( size(VarIn{ii},2) == prod(size(VarIn{ii})) )

       %--------------------------------------------------------
       recurse = ( recurse | strcmp(VarIn{ii},'-r') );
       follink = ( follink | strcmp(VarIn{ii},'-l') );
       verbose = ( verbose | strcmp(VarIn{ii},'-v') );
       nodispl = ( nodispl | strcmp(VarIn{ii},'-d') );
       nointro = ( nointro | strcmp(VarIn{ii},'-n') );

       %--------------------------------------------------------
       if strcmp( VarIn{ii}(1) , '.' );
          ext = cat(2,ext,VarIn(ii));
       end

       %--------------------------------------------------------
       if ( size(VarIn{ii},2) > 2 )

         switch VarIn{ii}(1:2)

          %---------------------------------------
          case '-K'

            eval(['kb = '  VarIn{ii}(3:end) ';'],'msgK = lasterr;'); 

          %---------------------------------------
          case '-M'

            eval(['kb = '  VarIn{ii}(3:end) '*1024;'],'msgM = lasterr;'); 


         end
         % switch

       end
       % if

    end 
    % if char

  end
  % ii

  %---------------------------------------------------------------------------
  if ~isempty(msgK)
      msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                      'Invalid Value for Option ''-K''.' , nl ,  msgK )
  end

  %---------------------------------------------------------------------------
  if ~isempty(msgM)
      msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                      'Invalid Value for Option ''-M''.' , nl ,  msgM )
  end

  %---------------------------------------------------------------------------
  if isempty(ext)

     ext = { ext0 };

  elseif any( strcmp( ext , '.' ) )

     ext = { char(ones(1,0)) };

  end

end
% if Nin


%---------------------------------------------------------

follink = ( follink | (~isunix) );

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


