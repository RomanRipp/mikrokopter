function msg = wrt_asc(file,key,val,form,comm,kmark,cmark)

% WRT_ASC  Write multiple numeric or character Arrays to one Ascii-File
%
% WRT_ASC( FILENAME, KEYWORD,   VALUE, FORMAT , COMMENT, ...
%                    KEYMARKER, COMMENTMARKER)
%
% Write Data from the CellArray VALUE into File, specified by
%   FILENAME. The Data could be NumericArrays, CharArrays or 
%    CellArrays, contains  NumericArrays and/or  Char(Cell)Arrays.
% By default, a new File is created, an existing File with the same Name
%   will be overwritten. If the FILENAME begins with an '@',
%   the Data will be append on an existing File.
%
% The Data will be seperated by a KEYWORD, like a VariableName, 
%  the KeyWord's are marked by the KEYMARKER (default: '#').
%  Use some UPPERCASE Letters in the KeyWords!
%  Use LOAD_ASC to load the Data specified by the KEYWORD's.
%
% COMMENT's for the Data are marked by the CommentMarker (default: '%').
%  If a KeyWord-String contains a CommentMarker, the second part from
%   this KeyWord will be interpreted as Commentary in LOAD_ASC.
%  Multiple Comments for one DataSet are given by a CharArray in the CellArray,
%    where the Rows of the CharArray are the different Comments.
%
% Note: If You give any Carr.Return (char(13)) or NewLine (char(10)) Character
%        in a CommentString, a CommentMarker has to follow!
%       Don't use a [']-Character as Marker, this Character is for Strings! 
% 
% FORMAT gives the Format for numeric Data, using by FPRINTF (default: '%g'). 
%  
% The Inputs FORMAT, COMMENT, KEYMARKER, COMMENTMARKER are optional.
%
% MESSAGE = WRT_ASC( ... ) returns the ErrorMessages.
%
% Note: WRT_ASC didn't produce ERROR's !
%
%
% DataTypes:  FILENAME   1 by n  CharArray (String)
%             KEYWORD    N by 1  CharArray or CellArray of Strings 
%             VALUE      M by 1  CellArray, M <= N 
%                                 contains CharArrays or NumericArrays
%             FORMAT     L by 1  CharArray or CellArray of Strings    L <= N
%             COMMENT    K by 1  CharArray or CellArray of CharArrays K <= N
%             KEYMARKER  1 by m  Char or String
%         COMMENTMARKER  1 by k  Char or String
%
%
% see also:  LOAD_ASC
%
%----------------------------------------------------------------------------
% Example:
% 
%% DATA 
%  cmap  = jet(12);            % 12 by 3 RGB-Matrice
%  cdat  = rand(5,6);
%  init  = { ['XX';'YY']    [ 20 150 ; -30 30 ]  0  ;
%             'Z'           [ 0 5000 ]           1  ;
%            {'Temp' 'Sal'} [ 0 35   ]           2    };
%  txt   = { 'Latitude' 'Longitude' };    
%
%% KEYWORDS
%  keys  = {'ColorMap';'CData %Random';'VInit';'Text' };
%   % the 2. KeyWord is "CData" , "Random" is a Commentary
%% FORMAT
%  form  = ['%6.2f';'%8.3f'];
%% COMMENT
%  comm1 = ['Jet         ';'use COLORMAP'];
%  comm3 = 'Labels for Axis'; 
%
%  % Write the Data into ASCII-File  
%  msg = wrt_asc('test.asc',keys,{cmap;cdat;init;txt}, ...
%                     form,{ comm1 ; '' ; '' ;comm3} ) 
%
%
%  key = {'ColorMap','CData','Text'};       
%
%  % Load the Data from ASCII-File
%  [msg,V]  = load_asc('test.asc',key); 
%
%   % produce an 3 by 4 CellArray, contains KeyWord, Value, Comments, Class
%
%  [msg,V,v1,v2,v3] = load_asc('test.asc',key); 
%
%   % Returns the Values in v1 .. v3
%



Nin  = nargin;
Nout = nargout;

nl = char(10);
 if ~isunix
  nl = [ char(13) nl ];
 end


msg = '';

msg0 = ' WRT_ASC:';

if Nin < 3
  msg = [ msg0 ' Not enough InputArguments.' ];
  return
end

if ~ischar(file)
 msg = [ msg0 ' FILENAME must be a String.' ];
  return
end

if Nin < 4
 form = '%g';
end
if Nin < 5
 comm = {''};
end
if Nin < 6
 kmark = '#';
end
if Nin < 7
 cmark = '%';
end


if iscellstr(kmark)
 kmark = kmark{1};
end

if ~ischar(kmark)  |  isempty(kmark)
   kmark = '#';
else
   kmark = rmblank(kmark,2);
end


if ~ischar(cmark)   |  isempty(cmark)
 cmark = '%';
else
 cmark = rmblank(cmark,2);
end



if ~( ischar(form) | iscellstr(form) )
 msg = [ msg nl(1:(end*(~isempty(msg))))  ...
         msg0 ' FORMAT must be a CharArray  or  CellArray of Strings ' ];
end 
if ~( ischar(key)  | iscellstr(key)  )
 msg = [ msg nl(1:(end*(~isempty(msg))))  ...
         msg0 ' KEYWORD must be a CharArray  or  CellArray of Strings ' ];
end 
if ~( ischar(comm) | iscellstr(comm) )
 msg = [ msg nl(1:(end*(~isempty(msg))))  ...
         msg0 ' COMMENT must be a CharArray  or  CellArray of Strings ' ];
end



% Form CellArrays
if ischar(form), form = cellstr(form); end ,  form = form(:);
if ischar(key) , key  = cellstr(key ); end ,  key  = key(:) ; 
if ischar(comm), comm = cellstr(comm); end ,  comm = comm(:);

if iscell   (val ), val  = val(:) ; end

if size(key,1) < size(val,1)
 msg = [ msg nl(1:(end*(~isempty(msg))))  ...
         msg0 ' KEYWORD''s must match the Size of VALUE''s ' ];
end

if ~isempty(msg)
 if Nout == 0, msg = ''; end
 return
end



N = size(key,1);


if size(form,1) == 1  
    
  if isempty(form{1}), form{1} = '%g'; end

  form = cellstr( char( ones(N,1) * form{1} ) );

end


% Fill Formats
form( size(form,1)+1 : N ) = ...
     cellstr( char( ones(N-size(form,1),1) * '%g' ) );

% Fill Values
val( size(val,1)+1 : N ) = cell( N-size(val,1) , 1 );

% Fill Comments
comm( size(comm,1)+1 : N ) = cell( N-size(comm,1) , 1 );



kmark = duplchar(kmark,'%\');        % DUPLCHAR, defined at the End of File
cmark = duplchar(cmark,'%\');



if file(1) == '@'
 file = file(2:length(file));
 mode = 'a';
else
 mode = 'wt';
end


fid  = fopen(file,mode);

if fid ~= -1
 
 for ii = 1:N

  if isempty(key{ii})
   key{ii} = sprintf('VARIABLE%4.4d',ii);
  else
   key{ii} = key{ii}(1,:);
  end

  % Cut left,right Blanks

  key{ii} = rmblank(key{ii},2);
 
   if isempty(key{ii})
    key{ii} = sprintf('VARIABLE%4.4d',ii);
   end

   fprintf(fid,[ kmark ' ' duplchar(key{ii},'%\')  nl ]); 

  if ~isempty( comm{ii} );
   for jj = 1:size(comm{ii},1)
    fprintf(fid,[ cmark ' '  duplchar(comm{ii}(jj,:),'%\')  nl ]);
   end
  end

   fprintf(fid, nl );

   var = val{ii};    
   


   if ischar(var)
   % Write CharArrays
   %-----------------------------------------------------------

    % Write Strings
    if isempty(var)
      fprintf(fid,[ ' '''' ' nl ]); 
    else
      for jj = 1:size(var,1)
        fprintf(fid,[ ' ''' duplchar(var(jj,:),'%\''') '''' nl ]);
      end     
    end

   elseif iscellstr(var)
   % Write CellArray of Char
   %------------------------------------------------------------

    var = var(:);

     for mm = 1:size(var,1)
        vv = var{mm};       % CharArray
        if isempty(vv) 
           fprintf(fid,[ ' '''' ' nl ]); 
        elseif size(vv,1) == 1
           fprintf(fid,[ ' ''' duplchar(vv,'%\''') '''  ' ]);
        else
           fprintf(fid,' [')
          for jj = 1:size(vv,1)
           dd = char( ';' + ( ']' - ';' ) * ( jj == size(vv,1) ) );
           fprintf(fid,[ ' ''' duplchar(vv(jj,:),'%\''') ''' ' dd  ]); 
          end 
        end      
        fprintf(fid,nl);
     end    
     % mm 


   elseif isnumeric(var)
   % Write Numerics, using Format
   %----------------------------------------------------------- 

    formm = form{ii};

    % Fill Format for Columns of Variable to write
    kk = find( formm == '%' );
    if length(kk) == 1 
      formm = [32 formm]'*ones(1,size(var,2)); 
      formm = char(reshape(formm,1,prod(size(formm))));
    end

    fprintf(fid,[formm nl],var');



   elseif iscell(var)   
   % WriteCellArrays
   %-----------------------------------------------------------


     for mm = 1 : size(var,1)

      for nn = 1 : size(var,2)

        vv = var{mm,nn};
         

        if ischar(vv)
         
          if isempty(vv)
            fprintf(fid,[ ' '''' ' ]);
          elseif size(vv,1) == 1 
            fprintf(fid,[ ' ''' duplchar(vv,'%\') ''' ' ]);
          else
            fprintf(fid,' [');
            for jj = 1:size(vv,1)
             dd = char( ';' + ( ']' - ';' ) * ( jj == size(vv,1) ) );
             fprintf(fid,[ ' ''' duplchar(vv(jj,:),'%\') ''' ' dd  ]); 
            end 
          end
       
        elseif iscellstr(vv)     
          
          fprintf(fid,' {');
          for jj = 1:size(vv,1)
            for ll = 1:size(vv,2)
              fprintf(fid,[ ' ''' duplchar(vv{jj,ll},'%\') ''' ' ]);
            end
            dd = char( ';' + ( '}' - ';' ) * ( jj == size(vv,1) ) );
            fprintf(fid, dd );
          end

        elseif isnumeric(vv)

           formm = form{ii};

           % Fill Format for Columns of Variable to write
           kk = find( formm == '%' );
           if length(kk) == 1 
            formm = [32 formm]'*ones(1,size(vv,2)); 
            formm = char(reshape(formm,1,prod(size(formm))));
           end


           si = size(vv);

           if max(si) == 1 
             fprintf(fid,formm,vv);
           elseif  ~isempty(vv)
            fprintf(fid,' [');
            for jj = 1:si(1)
              fprintf(fid,formm,vv(jj,:));
              dd = char( ';' + ( ']' - ';' ) * ( jj == size(vv,1) ) );
              fprintf(fid,dd);
            end
           end  
            
        end               
        % ischar | iscellstr | isnumeric (vv)

        fprintf(fid,'  ');

      end
      % nn

        fprintf(fid,nl);

     end
     % mm


   end 
   % ischar | iscellstr | isnumeric | iscell (var)  


  fprintf(fid,[nl nl]);


 end
 % ii

 fclose(fid);

else

 msg = [ msg0 ' Error open File: '  file  '.'  ];

end



%****************************************************************

function nstr = duplchar(str,str0)

% NEWSTR = DUPLCHAR(STR,CHAR), duplicates CHARs in String STR
%   for using FPRINTF, SPRINTF  (like '%' and '\')
%
% nstr = duplchar('%test\01','%')  gives nstr == '%%test\01'
% nstr = duplchar('%test\01','%\') gives nstr == '%%test\\01'
%

str  = str (:)';
str0 = str0(:);

for jj = 1:size(str0,1);

  ind = ones(1,size(str,2));

  ii = findstr(str,str0(jj));

  ind(ii) = 1+ind(ii);

  nstr = str0(jj) + zeros( 1 , size(str,2)+length(ii) );

  nstr( cumsum(ind) ) = str;

  nstr = char(nstr);

   str = nstr;
end


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


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
