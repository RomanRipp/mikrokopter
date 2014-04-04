function [msg,varargout] =  loadasc(file,key,kmark,cmark)

% LOADASC  Load multiple numeric or character Arrays from an Ascii-File
%
% !!! OLD VERSION of LOAD_ASC; please use the NEW version !!!
%
% ******************************************************************
%
% [MSG,V] = LOADASC( FILENAME, KEYWORD, KEYMARKER, COMMENTMARKER)
%
% Load Data, and Commentarys of Variables, specified by KEYWORD's, 
%  from the ASCII-File FILENAME, created by WRT_ASC.
%
% The KeyWord's are marked in the File with the     KEYMARKER (default: '#').
%  A KeyWord ends normally with NL, COMMENTMARKER or a additional Marker,
%  specified as second Value for KEYMARKER in a 2-Element CellString
% The Commnet's are marked in the File with the COMMENTMARKER (default: '%').
%
% The Inputs  KEYMARKER, COMMENTMARKER are optional.
%
% DataTypes:  FILENAME   1 by n  CharArray (String)
%             KEYWORD    N by m  CharArray or CellArray of Strings (m=1)
%             KEYMARKER  1 by k  Char or String
%         COMMENTMARKER  1 by l  Char or String
%
%
% ******************************************************************
% OUTPUT
%
% MSG are the ErrorMessages. If all ok, MSG == ''
%
% 1. Single OutPut **************************************************
%
%  [MSG,V] =  LOADASC( ... )  returns an [ N by 2 ]-CellArray for N KeyWords.
%
%  The 1. Column contains an [ M by 3 ]-CellArray for each KeyWord,
%    where M > 1, if the KeyWord was multiple found in the File. 
%
%    The 1. Column in this CellArray contains the Data of the Variables.
%     If the KeyWord was multiple found, the Data are stored in a CellArray.
%    The 2. Column contains the Comments for the KeyWord as 1 String, 
%     if the KeyWord was multiple found, the Strings are stored in a CellArray
%    The 3. Column contains ErrorMessages in reading the Data
%
%  The 2. Column of the [ N by 2 ]-CellArray contains an ErrorMessage,
%   if the KeyWord wasn't found in the File
%
%
%  If there are invalid InputArguments, an error opening the File or
%    the File didn't contains any KeyWord's,  V is the ErrorMessage. 
%
%
% 2. Multiple OutPut **************************************************
%
% [MSG, V1, V2, ... VN ]    = LOADASC( ... ) returns the Data of the N KeyWords
%  in the N Variables. 
% [MSG, V1, V2, ... VN, V ] = LOADASC( ... ) returns the Data of the N KeyWords
%  in the N Variables and additional the Data, Comments, ErrorMessages in 
%   the Variable V, described above.
%
%  If an KeyWord was multiple found in the File,
%    the Data are stored in CellArrays.
%
%  If the Number of OutPutArguments didn't correspondend with the Number of
%   KeyWord's, like described above, an ErrorMessage in V1 is returned.
%
%
% 3. Special OutPut **************************************************
%
%  [MSG, KeyWords , Comments ]=  LOADASC( FILENAME,  '@LOOK'       , ...
%                                          KEYMARKER, COMMENTMARKER       )
%
%  returns an Overview about ALL Keywords and Comments, containing in the File.
%   (KEYMARKER and COMMENTMARKER are optional)
%
% 
% Note: LOADASC didn't produce ERROR's !
%
%
%
% Example:
%  
%  cmap  = jet(12);            % 12 by 3 RGB-Matrice
%  cdat  = rand(5,6);
%  txt   = ['Lon';'Lat'];
%  keys  = {'ColorMap','CData %Random','Text'};
%  form  = ['%6.2f';'%8.3f'];
%  comm1 = ['Jet         ';'use COLORMAP'];
%  comm3 = 'Labels for Axis'; 
%  
%  % Write the Data into ASCII-File
%  msg = wrt_asc('test.asc',keys,{cmap;cdat;txt},form,{ comm1 ; '' ;comm3}) 
%
%
%   keys = {'ColorMap','CData','Text'};       
%
%  % Load the Data from ASCII-File
%            V  = load_asc('test.asc',keys); 
%   % produce an 3 by 2 CellArray, contains Data, Comments, Messages
%
%  [msg,V1,V2,V3,V] = load_asc('test.asc',keys); 
%   % Store the Data into V1 .. V3
%   %  CellArray VC contains Comments and Messages
%
%
%
% See also:  WRT_ASC
%






Nin  = nargin;
Nout = nargout - 1;  


varargout = cell(1,Nout);


nl = char(10);


msg = '';

msg0 = 'LOAD_ASC: ';


if Nin < 2
 msg = [ msg0 ' Not enough InputArguments.' ];
 return
end

if Nin < 3
 kmark = '#';
end
if Nin < 4
 cmark = '%';
end


if ~( ischar(file) & ~isempty(file) & ( prod(size(file)) == size(file,2) ) )
 msg = [ msg0 ' FILENAME must be a String.'  ];
 return
end 

if ischar(kmark)
  kmark = cellstr(kmark);
end
 
if ~iscellstr(kmark)
 msg = [ msg0 ' KEYMARKER must be a String or 2-Element CellString.'  ];
 return
end

if ~ischar(cmark)  |  isempty(cmark)
 cmark = '%';
else
 cmark = rmblank(cmark,2);
end


kend = '';
if prod(size(kmark)) >= 2
 kend = kmark{2};
end

kmark = rmblank(kmark{1},2);


if ~( ischar(key)  | iscellstr(key)  )
 msg = [ msg0 ' KEYWORD must be a CharArray  or  CellArray of Strings '  ];
 return
end 



% Look for '@LOOK'  request ==> OVERVIEW

is_overview = 0;

if ischar(key)
 is_overview = strcmp(key,'@LOOK');
end


%************************************************************

fid = fopen(file,'r');

if ( fid == -1 )  

 msg = [ msg0  'Can''t open File '  file  '.'  ];
 
 return

end



% Check Size of File (max 1 MB)

 dd = dir( file );

 if isempty(dd)
  dd = dir( which(file) );  % which(file) gives the full Name
                            %  [ PathName FileName ]
 end

 if dd.bytes > 1*2^20
  msg = [ msg0 'File '  file  '  too large.' ];
  fclose(fid);
  return
 end

 

 str = fread(fid,'char');

 fclose(fid);

 str = double(str(:)');

 ok = all( ( str ==  9 ) |  ...
           ( str == 10 ) |  ...
           ( str == 13 ) |  ...
           (  28 <= str  &   str <= 126 ) | ...
           ( 160 <= str  &   str <= 255 )        );

 if ~ok
  msg = [ msg0 'Invalid Characters in File '  file  '.' ];
  return
 end
 

 str = ['      '   char(str)  '          '];

 str = strrep(str,char([13 10]),char(10)); 




 Nstr = length(str);

 Nmark = length(kmark);

 starts = sort( findstr( str , kmark ) );

 % Starts marks the Sequences between the KeyWordMarker
 


 if isempty(starts)

  msg = [ msg0 'No KeyMarker found in File '   file   '.'  ];
  return

 end



 ww = warnstat; warning('off');


  starts = [ starts  Nstr+1 ];



if is_overview

%*********************************************************
% OVERVIEW
%---------------------------------------------------------

  % Search all POSSIBLE Keys
  
   glob = cell( length(starts)-1 , 2 );  % { Key  Comment }

   for ii = 1 : length(starts)-1

    ind = [ starts(ii) : starts(ii+1)-1 ];

    % Search for Commentarys
      cc = findstr(str(ind),cmark);                      % CommentMarker in ind
      nn = [ find( double(str(ind))==13 | ...
                   double(str(ind))==10 ) , length(ind) ];  % NewLines

      ke = findstr(str(ind),kend);         % MarkerEnd in ind


    % Extract Keyword

     fkey = char( str(ind([ 1+length(kmark) : min([ nn(1) cc ke ])-1 ])));

 
    % Cut left,right Blanks

    fkey = rmblank(fkey,2);

     glob{ii,1} = fkey; 

    % Extract Commentary
      cind = 0*ind;                                      % CommentIndex in ind
      if ~isempty(cc)
        for ff = 1:length(cc)
          nf = min(find(nn>cc(ff)));  % Following NewLine 
          cind(cc(ff)+length(cmark):nn(nf)) = 1+0*(cc(ff)+1:nn(nf));
        end
        nn_ind = ind(nn);
        ind = ind(find(cind));
        if any( ind(length(ind)) == nn_ind )
         % Remove NewLine at End
          ind(length(ind)) = [];
        end
        glob{ii,2} = char( str(ind) );
      else
        glob{ii,2} = '';
      end

   end
   % ii

   warning(ww);

   for ii = 1 : min(Nout,2)
    varargout{ii} = glob(:,ii);
   end
   
   return

end



%*********************************************************
% Read KeyWords and Values
%---------------------------------------------------------


  if ischar(key) 
     key  = cellstr(key );
  end 

  key  = key(:);

  % Founded KeyWords in the Sequences, specified by starts
  % [ Start  End  KeyWordIndex ]

  is_key = NaN * ones( size(key,1)*size(starts,2) , 3 );

      zz = 0;   % Counter for is_key

  % Messages for KeyWords
  kmsg = cell(size(key,1),1);

  nseq = zeros(size(key,1),1);  % Number of Sequences, matching the Keyword

  is_eval = zeros(size(key,1),1);  % If KeyWord is not a Variable in WorkSpace,
                                   %  it could be evaluate later

   % This are some Variables, we need later, 
   %  KeyWords with this Names couldn't be evaluate
   start = []; ind = []; ok = []; v = []; 
   s_i = [];  cind = []; nn = []; ke = []; jj = []; ff = []; 
   ok0 = []; ok1 = []; ok2 = [];  
   k1  = [];    lw = []; v1 = []; v2 = [];

  % Define ValueCellArray for KeyWords
  val = cell(size(key,1),2);   % 2. Column: KeyWordMessages (kmsg)

  val(:,2) = { '' };

  val(:,1) = {{ []  []  '' }};


  % Running through the KeyWords

  for kk = 1:size(key,1)

    key{kk} = rmblank(key{kk}(1,:),2);

    if isempty(key{kk})

      val{kk,2} = [ val{kk,2} nl(1:(end*(~isempty(val{kk,2}))))  ...
                    int2str(kk) '. KeyWord is empty.'  ];

    else
 
     ok = 0;  % Ok for KeyWord found in Sequences

     % Search through the Sequences for the KeyWord
     for ii = 1:size(starts,2)-1

       seq_str = str(starts(1,ii)+Nmark:starts(1,ii+1)-1);

            ff = findstr( seq_str , key{kk} );


       if ~isempty(ff)
       % Found KeyWord in the Sequence

         ff = ff(1);

         % Check for KeyWord at Begin
         if ff > 1
            ok1 = isempty(rmblank(seq_str(1:ff-1)));
         else
            ok1 = 1;
         end

         % Sequence without KeyWord
         seq_str = seq_str( ff+size(key{kk},2) : end );

       
         % End of KeyWord
         ke = size(seq_str,2);
         k1 = ke;
         for jj = { kend char(10) cmark }
            nn =  findstr(seq_str,jj{1});
            if ~isempty(nn)
               ke = min( ke , ( nn(1) - 1 + ...
                   size(kend,2) * strcmp(jj{1},kend) + strcmp(jj{1},char(10)) ) );
               k1 = min( k1 , nn(1)-1 );
            end
         end

         % Check empty End
         if k1 > 0 
            ok2 = isempty(rmblank(seq_str(1:k1)));
         else
            ok2 = 1;
         end


        if ok1 & ok2
        % Found the SINGLE!!! KeyWord at the BEGIN!!! of the Sequence 

         ok = 1; 
         zz = zz+1;

         is_key(zz,:) = cat(2,starts(1,ii)+Nmark+(ff-1)+size(key{kk},2)+ke , ...
                              starts(1,ii+1)-1 , kk );

         nseq(kk) = nseq(kk)+1;   % One Sequence more for the KeyWord

        end
        % KeyWord at Begin
       end
       % Found KeyWord

     end
     % ii Search through Sequences

     if ~ok
       val{kk,2} = [ val{kk,2} nl(1:(end*(~isempty(val{kk,2}))))  ...
        [ 'Didn''t found ' int2str(kk) '. KeyWord "' kmark  key{kk}   '"' nl ...
          ' in File '   file   '.' nl ] ];
     end
    
     % The KeyWord could be evaluate later, 
     %  if it's not a Variable in WorkSpace now.
      is_eval(kk) = ( exist(key{kk}) ~= 1 );


   end
   % key{ii} empty


   % Define CellArray of Sequences for the KeyWord in val
   if nseq(kk)
    val{kk,1} = cell(nseq(kk),3);  % 1. Column: Variable
                                   % 2. Column: Comments
                                   % 3. Column: ReadErrors
   end

  end
  % kk  Search through KeyWords

  
  % Reset Counter for Sequences of the KeyWords
  nseq = 0*nseq;


  if zz
  % Found any KeyWords in Sequences 
 
    is_key = is_key(1:zz,:);
 
    [start , s_i] = sort( is_key(:,1) );

    start = [ start ; length(str)+1 ];


    % Running through the Founded Sequences
    for ii = 1 : length(start)-1
 
      kk = is_key(s_i(ii),3); % KeyWordIndex

      nseq(kk) = nseq(kk)+1;  % Counter +1

      ind = ( is_key(s_i(ii),1) : is_key(s_i(ii),2) ); % Define the Sequence

      % Search for Commentarys
      jj = findstr(str(ind),cmark);                      % CommentMarker in ind

      nn = [ find( double(str(ind))==13 | ...
                   double(str(ind))==10 ) , length(ind) ];  % NewLines

      ke = findstr(str(ind),cmark);                      % MarkerEnd in ind

      cind_c = 0*ind;                                      % CommentIndex in ind
      cind_r = 0*ind;                                      %  RemoveIndex in ind

      if ~isempty(jj)
        for ff = 1:length(jj)
          nf = min(find(nn>jj(ff)));  % Following NewLine 
          cind_c(jj(ff)+1:nn(nf)  ) = 1 + 0 * ( jj(ff)+1 : nn(nf)   );  
          cind_r(jj(ff)+1:nn(nf)-1) = 1 + 0 * ( jj(ff)+1 : nn(nf)-1 );  
          % Without NewLine    "-1" !!!
        end
        ccind = ind(find(cind_c));
        if any( ccind(length(ccind)) == ind(nn) )
         % Remove NewLine at End
          ccind(length(ccind)) = [];
        end
        val{kk,1}{nseq(kk),2} = char( str(ccind) );
        ind([jj find(cind_r)]) = [];   % Remove Commentarys and Marker
      end            

      % We will get some Problems:
      %
      % eval([' v = [' str '];' ]), % with str = "'XT' 0 360"
      % gives:  v = 'XT h'
      %
      % Therefore check with Warning (if invalid CharConversation) 
      %  AND with Datatypes using CellArray.

      % Reset WarningMessage

      lastwarn('');

       ok1 = 1;
       % Try Numeric- or CharArray
       eval([ 'v1 = [ ' ...
             char(str(ind))    ' ]; ' ] , 'ok1 = 0; v1 = ''''; ')
       ok1 = ok1*isempty(lastwarn);

       ok2 = 1;
       % Try CellArray
       eval([ 'v2 = { ' ...
            char(str(ind))    ' }; ' ] , 'ok2 = 0; v2 = {}; ')

       is_numeric = 0;
       for jj = 1 : size(v2,1)*size(v2,2)
        is_numeric = ( is_numeric | isnumeric(v2{jj}) );
       end

       ok1 = ok1 * (~( is_numeric  &  ischar(v1)  &  ok2 ));
      
       if     ok1,  v = v1; 
       elseif ok2,  v = v2; end


      if  ok1 | ok2  

        val{kk,1}{nseq(kk),1} = v;

        if is_eval(kk)
        % Evaluate the Variable in KeyWord, 
        %  to use as Variable for following Sequences
          eval( [ key{kk}   ' = v;'] , 'ok=0;' )
        end

      else
       val{kk,1}{nseq(kk),3} = [                           ...
          [ 'Error reading KeyWord "'  kmark  key{kk} '" :'   ...
            nl lasterr ]  ];
      end
      % ok

    end
    % ii Running through Sequences


  end
  % zz  found any KeyWords in the Sequences


  warning(ww);


   if Nout == 1

     varargout{1} = val;

   elseif Nout > 1  

     for ii = 1 : min(size(key,1),Nout)

      if size( val{ii,1}(:,1) , 1 ) == 1
       varargout{ii} = val{ii,1}{1};
      else
       varargout{ii} = val{ii,1}(:,1);
      end
    
     end
  
     if Nout >=  size(key,1)+1
       varargout{size(key,1)+1} = val;
     end

   end


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


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


%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ww = warnstat

% WARNSTAT  Returns global WarningStatus
%
%  WARNSTAT returns the Status of WARNING
%
% Matlab R<13   WARNING
% Matlab R>12   WARNING for Identifier ALL
%

ww = warning;

if isstruct(ww)   % New Matlab R>12 Syntax
   try
      id = strcmp({ww.identifier},'all');
      if any(id)
         id = find(id);
         ww = ww(id(1)).state;
      else
         ww = '';
      end
   catch
      ww = '';
   end
elseif ~chkstr(ww)
   ww = '';
end

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

